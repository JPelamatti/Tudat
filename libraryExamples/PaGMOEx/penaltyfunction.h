#ifndef PENALTYFUNCTION_H
#define PENALTYFUNCTION_H

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Eigen/Core>
#include <cmath>


double PenaltyFunction(std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > stateVectors);

#endif // PENALTYFUNCTION_H
