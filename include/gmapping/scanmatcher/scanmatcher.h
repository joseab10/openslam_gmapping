#ifndef SCANMATCHER_H
#define SCANMATCHER_H

#include "gmapping/scanmatcher/icp.h"
#include "gmapping/scanmatcher/smmap.h"
#include "gmapping/scanmatcher/gridlinetraversal.h"
#include <gmapping/utils/macro_params.h>
#include <gmapping/utils/stat.h>
#include <iostream>
#include <gmapping/utils/gvalues.h>

#define LASER_MAXBEAMS 2048

namespace GMapping {

    /*
    enum MapModel {
        ReflectionModel,
        ExpDecayModel
    };

    enum ParticleWeighting {
        ClosestMeanHitLikelihood,
        ForwardSensorModel,
        MeasurementLikelihood
    };
    */

    class ScanMatcher {
    public:
        typedef Covariance3 CovarianceMatrix;

        ScanMatcher();

        ~ScanMatcher();

        double icpOptimize(OrientedPoint &pnew, const ScanMatcherMap &map, const OrientedPoint &p,
                           const double *readings) const;

        double
        optimize(OrientedPoint &pnew, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        double optimize(OrientedPoint &mean, CovarianceMatrix &cov, const ScanMatcherMap &map, const OrientedPoint &p,
                        const double *readings) const;

        double registerScan(ScanMatcherMap &map, const OrientedPoint &p, const double *readings);

        void setLaserParameters
                (unsigned int beams, double *angles, const OrientedPoint &lpose);

        enum ParticleWeighting {
            ClosestMeanHitLikelihood,
            ForwardSensorModel,
            MeasurementLikelihood
        };

        void setMatchingParameters
                (double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations,
                 double likelihoodSigma = 1, unsigned int likelihoodSkip = 0,
                 ScanMatcherMap::MapModel mapModel = ScanMatcherMap::MapModel::ReflectionModel,
                 ParticleWeighting particleWeighting = ParticleWeighting::ClosestMeanHitLikelihood);

        void invalidateActiveArea();

        void computeActiveArea(ScanMatcherMap &map, const OrientedPoint &p, const double *readings);

        inline double
        icpStep(OrientedPoint &pret, const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;


        inline double score(const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const;

        inline unsigned int likelihoodAndScore(double &s, double &l, const ScanMatcherMap &map, const OrientedPoint &p,
                                               const double *readings) const;

        double likelihood(double &lmax, OrientedPoint &mean, CovarianceMatrix &cov, const ScanMatcherMap &map,
                          const OrientedPoint &p, const double *readings);

        double likelihood(double &_lmax, OrientedPoint &_mean, CovarianceMatrix &_cov, const ScanMatcherMap &map,
                          const OrientedPoint &p, Gaussian3 &odometry, const double *readings, double gain = 180.);

        inline const double *laserAngles() const { return m_laserAngles; }

        inline unsigned int laserBeams() const { return m_laserBeams; }

        static const double nullLikelihood;
    protected:
        //state of the matcher
        bool m_activeAreaComputed;

        /**laser parameters*/
        unsigned int m_laserBeams;
        double m_laserAngles[LASER_MAXBEAMS];
        //OrientedPoint m_laserPose;

        double computeCellR(const ScanMatcherMap &map, Point beamStart, Point beamEnd, IntPoint cell) const;

    PARAM_SET_GET(OrientedPoint, laserPose, protected, public, public)

    PARAM_SET_GET(double, laserMaxRange, protected, public, public)
        /**scan_matcher parameters*/
    PARAM_SET_GET(double, usableRange, protected, public, public)

    PARAM_SET_GET(double, gaussianSigma, protected, public, public)

    PARAM_SET_GET(double, likelihoodSigma, protected, public, public)

    PARAM_SET_GET(int, kernelSize, protected, public, public)

    PARAM_SET_GET(double, optAngularDelta, protected, public, public)

    PARAM_SET_GET(double, optLinearDelta, protected, public, public)

    PARAM_SET_GET(unsigned int, optRecursiveIterations, protected, public, public)

    PARAM_SET_GET(unsigned int, likelihoodSkip, protected, public, public)

    PARAM_SET_GET(double, llsamplerange, protected, public, public)

    PARAM_SET_GET(double, llsamplestep, protected, public, public)

    PARAM_SET_GET(double, lasamplerange, protected, public, public)

    PARAM_SET_GET(double, lasamplestep, protected, public, public)

    PARAM_SET_GET(bool, generateMap, protected, public, public)

    PARAM_SET_GET(ScanMatcherMap::MapModel, mapModel, protected, public, public)

    PARAM_SET_GET(ParticleWeighting, particleWeighting, protected, public, public)

    PARAM_SET_GET(double, enlargeStep, protected, public, public)

    PARAM_SET_GET(double, fullnessThreshold, protected, public, public)

    PARAM_SET_GET(double, angularOdometryReliability, protected, public, public)

    PARAM_SET_GET(double, linearOdometryReliability, protected, public, public)

    PARAM_SET_GET(double, freeCellRatio, protected, public, public)

    PARAM_SET_GET(unsigned int, initialBeamsSkip, protected, public, public)

        // allocate this large array only once
        IntPoint *m_linePoints;
    };

    inline double ScanMatcher::icpStep(OrientedPoint &pret, const ScanMatcherMap &map, const OrientedPoint &p,
                                       const double *readings) const {
        const double *angle = m_laserAngles + m_initialBeamsSkip;
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;
        unsigned int skip = 0;
        double freeDelta = map.getDelta() * m_freeCellRatio;
        std::list<PointPair> pairs;

        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++) {
            skip++;
            skip = skip > m_likelihoodSkip ? 0 : skip;
            if (*r > m_usableRange || *r == 0.0) continue;
            if (skip) continue;
            Point phit = lp;
            phit.x += *r * cos(lp.theta + *angle);
            phit.y += *r * sin(lp.theta + *angle);
            IntPoint iphit = map.world2map(phit);
            Point pfree = lp;
            pfree.x += (*r - map.getDelta() * freeDelta) * cos(lp.theta + *angle);
            pfree.y += (*r - map.getDelta() * freeDelta) * sin(lp.theta + *angle);
            pfree = pfree - phit;
            IntPoint ipfree = map.world2map(pfree);
            bool found = false;
            Point bestMu(0., 0.);
            Point bestCell(0., 0.);
            for (int xx = -m_kernelSize; xx <= m_kernelSize; xx++)
                for (int yy = -m_kernelSize; yy <= m_kernelSize; yy++) {
                    IntPoint pr = iphit + IntPoint(xx, yy);
                    IntPoint pf = pr + ipfree;
                    //AccessibilityState s=map.storage().cellState(pr);
                    //if (s&Inside && s&Allocated){
                    const PointAccumulator &cell = map.cell(pr);
                    const PointAccumulator &fcell = map.cell(pf);
                    if (((double) cell) > m_fullnessThreshold && ((double) fcell) < m_fullnessThreshold) {
                        Point mu = phit - cell.mean();
                        if (!found) {
                            bestMu = mu;
                            bestCell = cell.mean();
                            found = true;
                        } else if ((mu * mu) < (bestMu * bestMu)) {
                            bestMu = mu;
                            bestCell = cell.mean();
                        }

                    }
                    //}
                }
            if (found) {
                pairs.push_back(std::make_pair(phit, bestCell));
                //std::cerr << "(" << phit.x-bestCell.x << "," << phit.y-bestCell.y << ") ";
            }
            //std::cerr << std::endl;
        }

        OrientedPoint result(0, 0, 0);
        //double icpError=icpNonlinearStep(result,pairs);
        std::cerr << "result(" << pairs.size() << ")=" << result.x << " " << result.y << " " << result.theta
                  << std::endl;
        pret.x = p.x + result.x;
        pret.y = p.y + result.y;
        pret.theta = p.theta + result.theta;
        pret.theta = atan2(sin(pret.theta), cos(pret.theta));
        return score(map, p, readings);
    }

    inline double ScanMatcher::score(const ScanMatcherMap &map, const OrientedPoint &p, const double *readings) const {
        double s = 0;
        const double *angle = m_laserAngles + m_initialBeamsSkip;
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;
        unsigned int skip = 0;
        double freeDelta = map.getDelta() * m_freeCellRatio;
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++) {
            skip++;
            skip = skip > m_likelihoodSkip ? 0 : skip;
            if (skip || *r > m_usableRange || *r == 0.0) continue;
            Point phit = lp;
            phit.x += *r * cos(lp.theta + *angle);
            phit.y += *r * sin(lp.theta + *angle);
            IntPoint iphit = map.world2map(phit);
            Point pfree = lp;
            pfree.x += (*r - map.getDelta() * freeDelta) * cos(lp.theta + *angle);
            pfree.y += (*r - map.getDelta() * freeDelta) * sin(lp.theta + *angle);
            pfree = pfree - phit;
            IntPoint ipfree = map.world2map(pfree);
            bool found = false;
            Point bestMu(0., 0.);
            for (int xx = -m_kernelSize; xx <= m_kernelSize; xx++)
                for (int yy = -m_kernelSize; yy <= m_kernelSize; yy++) {
                    IntPoint pr = iphit + IntPoint(xx, yy);
                    IntPoint pf = pr + ipfree;
                    //AccessibilityState s=map.storage().cellState(pr);
                    //if (s&Inside && s&Allocated){
                    const PointAccumulator &cell = map.cell(pr);
                    const PointAccumulator &fcell = map.cell(pf);
                    if (((double) cell) > m_fullnessThreshold && ((double) fcell) < m_fullnessThreshold) {
                        Point mu = phit - cell.mean();
                        if (!found) {
                            bestMu = mu;
                            found = true;
                        } else
                            bestMu = (mu * mu) < (bestMu * bestMu) ? mu : bestMu;
                    }
                    //}
                }
            if (found)
                s += exp(-1. / m_gaussianSigma * bestMu * bestMu);
        }
        return s;
    }

    inline unsigned int
    ScanMatcher::likelihoodAndScore(double &s, double &l, const ScanMatcherMap &map, const OrientedPoint &p,
                                    const double *readings) const {
        using namespace std;
        l = 0;
        s = 0;
        const double *angle = m_laserAngles + m_initialBeamsSkip;
        OrientedPoint lp = p;
        lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
        lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
        lp.theta += m_laserPose.theta;
        IntPoint ilp = map.world2map(lp);

        double noHit = nullLikelihood / (m_likelihoodSigma);
        unsigned int skip = 0;
        unsigned int c = 0;
        double freeDelta = map.getDelta() * m_freeCellRatio;
        for (const double *r = readings + m_initialBeamsSkip; r < readings + m_laserBeams; r++, angle++) {
            skip++;
            skip = skip > m_likelihoodSkip ? 0 : skip;
            //if (*r > m_usableRange) continue;
            double range = *r;
            bool out_of_range = false;

            if (range > m_usableRange) {
                if (m_particleWeighting == ClosestMeanHitLikelihood)
                    continue;
                else {
                    range = m_usableRange;
                    out_of_range = true;
                }
            }

            if (skip) continue;
            Point phit = lp;
            phit.x += range * cos(lp.theta + *angle);
            phit.y += range * sin(lp.theta + *angle);
            IntPoint iphit = map.world2map(phit);

            /*
             * Original Particle weight computation method from openslam_gmapping.
             * It computes a point pfree which is (roughly) one grid cell before the hit point in the
             * direction of the beam.
             * Then, for a given window around the hit point (+/- kernelSize in x and y), it checks if each cell
             * is occupied, and the cell in pfree relative to it is free according to a threshold value.
             * If it is the case, then the difference between the actual hitpoint and the mean of all hits of said cell
             * is computed. Finally, the minimum of all distances from hits and means is used to compute a
             * log likelihood from a gaussian, and added to the particle as it's weight.
             */
            if (m_particleWeighting == ClosestMeanHitLikelihood) {


                Point pfree = lp;
                pfree.x += (range - freeDelta) * cos(lp.theta + *angle);
                pfree.y += (range - freeDelta) * sin(lp.theta + *angle);
                pfree = pfree - phit;
                IntPoint ipfree = map.world2map(pfree);
                bool found = false;
                Point bestMu(0., 0.);
                for (int xx = -m_kernelSize; xx <= m_kernelSize; xx++)
                    for (int yy = -m_kernelSize; yy <= m_kernelSize; yy++) {
                        IntPoint pr = iphit + IntPoint(xx, yy);
                        IntPoint pf = pr + ipfree;
                        //AccessibilityState s=map.storage().cellState(pr);
                        //if (s&Inside && s&Allocated){
                        const PointAccumulator &cell = map.cell(pr);
                        const PointAccumulator &fcell = map.cell(pf);
                        if (((double) cell) > m_fullnessThreshold && ((double) fcell) < m_fullnessThreshold) {
                            Point mu = phit - cell.mean();
                            if (!found) {
                                bestMu = mu;
                                found = true;
                            } else
                                bestMu = (mu * mu) < (bestMu * bestMu) ? mu : bestMu;
                        }
                        //}
                    }
                if (found) {
                    s += exp(-1. / m_gaussianSigma * bestMu * bestMu);
                    c++;
                }
                if (!skip) {
                    double f = (-1. / m_likelihoodSigma) * (bestMu * bestMu);
                    l += (found) ? f : noHit;
                }
            }
            else {
                /*
                 * For all other methods, we will consider all cells crossed by the beam and not just the endpoint.
                 */
                GridLineTraversalLine line;
                line.points = m_linePoints;
                GridLineTraversal::gridLine(ilp, iphit, &line);

                double alpha_prior = map.getAlpha();
                double beta_prior  = map.getBeta();

                int i = 0;
                // For all the cells that the beam travelled through (misses)
                for (i = 0; i < line.num_points - 1 + out_of_range; i++) {
                    const IntPoint i_miss_cell = line.points[i];
                    const PointAccumulator &miss_cell = map.cell(i_miss_cell);
                    int Hi = miss_cell.n - miss_cell.n_inc;

                    if (m_mapModel == ScanMatcherMap::MapModel::ReflectionModel) {
                        int Mi = miss_cell.visits - Hi - miss_cell.visits_inc; // Subtract 1 to get the previous number of misses

                        double denominator = Hi + alpha_prior + Mi + beta_prior;

                        if (denominator == 0.0)
                            l += noHit;
                        else {
                            if (m_particleWeighting == MeasurementLikelihood)
                                // p = (Mi + beta) / denominator
                                l += log(Mi + beta_prior) - log(denominator);
                            else if (m_particleWeighting == ForwardSensorModel)
                                // p = 1 - (Hi + alpha) / denominator
                                l += log(1 - (Hi + alpha_prior) / denominator);
                        }
                    }
                    else if (m_mapModel == ScanMatcherMap::MapModel::ExpDecayModel){
                        double ri = computeCellR(map, lp, phit, i_miss_cell);
                        double Ri = miss_cell.R - miss_cell.R_inc; // subtract ri to get the previous Ri

                        if (m_particleWeighting == MeasurementLikelihood) {
                            double denominator = Ri + beta_prior + ri;
                            if (denominator == 0.0 || Ri + beta_prior == 0.0)
                                l += noHit;
                            else
                                // p = pow(((Ri + beta) / denominator), (Hi + alpha))
                                l += (Hi + alpha_prior) * (log(Ri + beta_prior) - log(denominator));
                        }
                        else if (m_particleWeighting == ForwardSensorModel) {
                            if (Ri == 0.0)
                                l += noHit;
                            else {
                                //double lambda = (Hi + alpha) / denominator;
                                double lambda = Hi / Ri;
                                // p = exp(-(lambda * ri))
                                l += -(lambda * ri);
                            }
                        }
                    }
                }
                // For the endpoint cell
                if (! out_of_range){
                    IntPoint i_hit_cell = line.points[line.num_points -1];
                    const PointAccumulator &hit_cell = map.cell(i_hit_cell);

                    int Hi = hit_cell.n - hit_cell.n_inc; // Subtract 1 to get the previous number of hits

                    if (m_mapModel == ScanMatcherMap::MapModel::ReflectionModel) {
                        int Mi = hit_cell.visits - Hi - hit_cell.visits_inc; // Subtract 1 to get the previous number of misses

                        if (m_particleWeighting == MeasurementLikelihood) {
                            double denominator = Hi + alpha_prior + Mi + beta_prior;
                            if (denominator == 0)
                                l += noHit;
                            else
                                // p = (Hi + alpha) / denominator
                                l += log(Hi + alpha_prior) - log(denominator);
                        }
                        else if (m_particleWeighting == ForwardSensorModel) {
                            double denominator = Hi + Mi;
                            if (denominator == 0)
                                l += noHit;
                            else
                                // p = Hi / denominator
                                l += log(Hi) - log(denominator);
                        }
                    }
                    else if (m_mapModel == ScanMatcherMap::MapModel::ExpDecayModel){
                        double ri = computeCellR(map, lp, phit, i_hit_cell);
                        double Ri = hit_cell.R - hit_cell.R_inc; // subtract ri to get the previous Ri

                        if (m_particleWeighting == MeasurementLikelihood) {
                            double denominator = Ri + beta_prior + ri;
                            if (denominator == 0.0)// || Hi + alpha == 0 || Ri + beta == 0)
                                l += noHit;
                            else {
                                double n1 = Ri + beta_prior;
                                double n2 = Hi + alpha_prior;
                                // p = pow(((Ri + beta) / denominator), (Hi + alpha)) * ((Hi + alpha) / denominator)
                                l += log(n2) + (n2 * (log(n1) - log(denominator))) - log(denominator);
                            }
                        }
                        else if (m_particleWeighting == ForwardSensorModel){
                            double denominator = Ri;// + beta;

                            if (!Hi || denominator == 0.0)
                                l += noHit;
                            else {
                                //double lambda = (Hi + alpha) / denominator;
                                double lambda = Hi / denominator;
                                // p = lambda * exp(-(lambda * ri))
                                l += log(lambda) - (lambda * ri);
                            }
                        }

                    }
                }
                PointAccumulator cell = map.cell(iphit);
                if (! out_of_range && cell.n){
                    Point mu = phit - cell.mean();
                    s += exp(-1. / m_gaussianSigma * mu * mu);
                    c++;
                }

            }
        }
        return c;
    }

};

#endif
