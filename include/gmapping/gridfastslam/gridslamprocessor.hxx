
#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x==FP_NAN)
#endif

/**Just scan match every single particle.
If the scan matching fails, the particle gets a default likelihood.*/
inline void GridSlamProcessor::scanMatch(const double *plainReading) {

    if (!m_doImprovePose)
        return;

    // sample a new pose from each scan in the reference
    double sumScore = 0;

    if (m_outputStream.is_open()) {
        m_outputStream << "OPT_SCORE " << m_particles.size() << " ";
        //m_outputStream << std::setiosflags(std::ios::fixed) << std::setprecision(6);
    }

    for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
        OrientedPoint corrected;
        double score = 0, l, s;

        score = m_matcher.optimize(corrected, it->map, it->pose, plainReading);
        if(m_outputStream.is_open())
            m_outputStream << score << " " << (score > m_minimumScore) << " ";

        //    it->pose=corrected;
        if (score > m_minimumScore) {
            it->pose = corrected;
        } else {
            if (m_infoStream) {
                m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l << std::endl;
                m_infoStream << "lp:" << m_lastPartPose.x << " " << m_lastPartPose.y << " " << m_lastPartPose.theta
                             << std::endl;
                m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " " << m_odoPose.theta << std::endl;
            }
        }

        sumScore += score;
    }
    if (m_outputStream.is_open())
        m_outputStream << std::endl;
    if (m_infoStream)
        m_infoStream << "Average Scan Matching Score=" << sumScore / m_particles.size() << std::endl;
}

void GridSlamProcessor::weightParticles(const double *plainReading){
    if (m_particleWeighting == ScanMatcher::ParticleWeighting::ClosestMeanHitLikelihood){
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {

            double l, s;

            m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);

            it->weight += l;
            it->weightSum += l;
        }
    }
    else if(m_particleWeighting == ScanMatcher::ParticleWeighting::MeasurementLikelihood){
        int num_particles = m_particles.size();
        int initial_beam_skip = m_matcher.getinitialBeamsSkip();
        int likelihood_skip = m_matcher.getlikelihoodSkip();
        int skip = 0;

        const double *angles = m_matcher.laserAngles();
        double max_range = std::max(m_matcher.getusableRange(), m_matcher.getlaserMaxRange());
        // Get the maximum number of cells that a scan could have
        double max_beam_cells = max_range / m_delta;

        // Using 180 beams with 20.0m max_range and a map resolution of 0.05m as reference with
        // 18 batches, i.e.: batch size of 10 beams.
        double num_batches = max_beam_cells * (180 / 10) / (20.0 / 0.05);
        int batch_beams = m_beams / num_batches;
        // Compute the number of beams a batch should process before normalizing the weights.
        batch_beams = std::max(1, batch_beams);

        // initialize particle weights with 1
        for (ParticleVector::iterator it=m_particles.begin(); it < m_particles.end(); it++) {
            it->weight = 0;
            it->weightSum = 0;
        }

        std::vector<double> weights(num_particles, 0.0);

        // For each batch of beams
        for (int b = initial_beam_skip; b < m_beams; b += batch_beams){

            std::vector<double> log_likelihoods(num_particles, 0.0);
            // For each beam in the batch
            for (int i = 0; i < batch_beams && b + i < m_beams; i++) {
                skip++;
                skip = skip > likelihood_skip ? 0 : skip;
                if (skip) continue;

                double reading_range = plainReading[b + i];
                double reading_bearing = angles[b + i];

                // For each particle
                for (int part = 0; part < num_particles; part++) {

                    if (log_likelihoods[part] == -std::numeric_limits<double>::max())
                        continue;

                    double tmp_likelihood = m_matcher.measurementLikelihood(m_particles[part].pose,
                                                                            reading_range, reading_bearing,
                                                                            m_particles[part].map);

                    if (isnan(tmp_likelihood))
                        log_likelihoods[part] = -std::numeric_limits<double>::max();
                    else
                        log_likelihoods[part] += tmp_likelihood;
                }
            }

            softmax_normalize(log_likelihoods);

            for (int j = 0; j < num_particles; j++)
                weights[j] += log(log_likelihoods[j]);

            softmax_normalize(weights);

            for (int j = 0; j < num_particles; j++)
                weights[j] = log(weights[j]);

        }
        for (int j = 0; j < num_particles; j++)
            m_particles[j].weight = weights[j];
    }

    for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
        //set up the selective copy of the active area
        //by detaching the areas that will be updated
        m_matcher.invalidateActiveArea();
        m_matcher.computeActiveArea(it->map, it->pose, plainReading);
    }
}

inline void GridSlamProcessor::normalize() {

    m_weights.clear();
    for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
        m_weights.push_back(it->weight);
    }
    double gain = 1;

    if (m_particleWeighting == ScanMatcher::ParticleWeighting::ClosestMeanHitLikelihood)
        gain = 1. / (m_obsSigmaGain * m_particles.size());

    m_neff = softmax_normalize(m_weights, gain);

}

inline double GridSlamProcessor::softmax_normalize(std::vector<double> &values, double gain){
    // Perform a Softmax over the particle weights to normalize them an make them add up to one.
    double v_max = -std::numeric_limits<double>::max();

    // Find max value
    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
        v_max = *it > v_max ? *it : v_max;
    }

    double v_acc = 0;

    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
        // Exponentiate all weights minus the maximum for numerical stability.
        *it = exp(gain * (*it - v_max));
    }

    return linear_normalize(values);
}

inline double GridSlamProcessor::linear_normalize(std::vector<double> &values){

    double v_acc = 0, v_cnt = 0;

    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
        v_acc += *it;
        v_cnt ++;
    }

    if (v_acc != 0)
        v_acc = 1 / v_acc;
    if (v_cnt != 0)
        v_cnt = 1 / v_cnt;

    double neff = 0;
    // Divide all weights by the sum of all values.
    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
        if (v_acc == 0)
            *it = v_cnt;
        else
            *it *= v_acc;
        neff += (*it) * (*it);
    }

    return 1 / neff;
}

inline bool GridSlamProcessor::resample(const double *plainReading, int adaptSize, const RangeReading *reading) {

    bool hasResampled = false;

    TNodeVector oldGeneration;
    for (unsigned int i = 0; i < m_particles.size(); i++) {
        oldGeneration.push_back(m_particles[i].node);
    }

    if (m_neff < m_resampleThreshold * m_particles.size()) {

        if (m_infoStream)
            m_infoStream << "*************RESAMPLE***************" << std::endl;

        uniform_resampler<double, double> resampler;
        m_indexes = resampler.resampleIndexes(m_weights, adaptSize);

        if (m_outputStream.is_open()) {
            m_outputStream << "RESAMPLE " << m_indexes.size() << " ";
            for (std::vector<unsigned int>::const_iterator it = m_indexes.begin(); it != m_indexes.end(); it++) {
                m_outputStream << *it << " ";
            }
            m_outputStream << std::endl;
        }

        onResampleUpdate();
        //BEGIN: BUILDING TREE
        ParticleVector temp;
        unsigned int j = 0;
        std::vector<unsigned int> deletedParticles;        //this is for deleteing the particles which have been resampled away.

        //		cerr << "Existing Nodes:" ;
        for (unsigned int i = 0; i < m_indexes.size(); i++) {
            //			cerr << " " << m_indexes[i];
            while (j < m_indexes[i]) {
                deletedParticles.push_back(j);
                j++;
            }
            if (j == m_indexes[i])
                j++;
            Particle &p = m_particles[m_indexes[i]];
            TNode *node = 0;
            TNode *oldNode = oldGeneration[m_indexes[i]];
            //			cerr << i << "->" << m_indexes[i] << "B("<<oldNode->childs <<") ";
            node = new TNode(p.pose, 0, oldNode, 0);
            //node->reading=0;
            node->reading = reading;
            //			cerr << "A("<<node->parent->childs <<") " <<endl;

            temp.push_back(p);
            temp.back().node = node;
            temp.back().previousIndex = m_indexes[i];
        }
        while (j < m_indexes.size()) {
            deletedParticles.push_back(j);
            j++;
        }
        //		cerr << endl;
        std::cerr << "Deleting Nodes:";
        for (unsigned int i = 0; i < deletedParticles.size(); i++) {
            std::cerr << " " << deletedParticles[i];
            delete m_particles[deletedParticles[i]].node;
            m_particles[deletedParticles[i]].node = 0;
        }
        std::cerr << " Done" << std::endl;

        //END: BUILDING TREE
        std::cerr << "Deleting old particles...";
        m_particles.clear();
        std::cerr << "Done" << std::endl;
        std::cerr << "Copying Particles";

        if (m_doMapUpdate)
                     std::cerr << "and  Registering  scans...";

        for (ParticleVector::iterator it = temp.begin(); it != temp.end(); it++) {
            it->setWeight(0);
            if (m_doMapUpdate) {
                m_matcher.invalidateActiveArea();
                m_matcher.registerScan(it->map, it->pose, plainReading);
            }
            m_particles.push_back(*it);
        }
        std::cerr << " Done" << std::endl;
        hasResampled = true;
    } else {
        int index = 0;
        std::cerr << "Registering Scans:";
        TNodeVector::iterator node_it = oldGeneration.begin();
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
            //create a new node in the particle tree and add it to the old tree
            //BEGIN: BUILDING TREE
            TNode *node = 0;
            node = new TNode(it->pose, 0.0, *node_it, 0);

            //node->reading=0;
            node->reading = reading;
            it->node = node;

            //END: BUILDING TREE
            if (m_doMapUpdate) {
                m_matcher.invalidateActiveArea();
                m_matcher.registerScan(it->map, it->pose, plainReading);
            }
            it->previousIndex = index;
            index++;
            node_it++;

        }
        std::cerr << "Done" << std::endl;

    }
    //END: BUILDING TREE

    return hasResampled;
}
