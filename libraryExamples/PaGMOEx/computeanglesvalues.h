#ifndef COMPUTEANGLESVALUES_H
#define COMPUTEANGLESVALUES_H

#include <iostream>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"


// External libraries: Eigen
#include <Eigen/Core>


Eigen::Vector2d computeAngles(double t,std::vector<std::vector<double>> coeff, double period);


#endif // COMPUTEANGLESVALUES_H
