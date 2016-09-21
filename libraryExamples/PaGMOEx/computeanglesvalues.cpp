#include <iostream>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"


// External libraries: Eigen
#include <Eigen/Core>

// Header
#include "/home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/computeanglesvalues.h"


// Function that computes the values of the attitude angles at a given time given the polynomial coefficients
// Works for an integration period of a year


// MUST CHECK EQUATIONS

Eigen::Vector2d computeAngles(double t,std::vector<std::vector<double>> coeff, double period)
{
    std::vector<double> coneAngleCoeffs, clockAngleCoeffs;

    double coneAngleValue;
    double clockAngleValue;
    int numOfRev,n;
    double deltaT;

    numOfRev = t/period;
    deltaT = numOfRev*period;

    coneAngleCoeffs = coeff.at(0);
    clockAngleCoeffs = coeff.at(1);

    if (t <= period*365/4)
    {
      n = 0;
      coneAngleValue = coneAngleCoeffs[0+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
              + coneAngleCoeffs[1+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
              + coneAngleCoeffs[2+n]*(t-deltaT+coneAngleCoeffs[4+n]) + coneAngleCoeffs[3+n];

      clockAngleValue = clockAngleCoeffs[0+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
              + clockAngleCoeffs[1+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
              + clockAngleCoeffs[2+n]*(t-deltaT+clockAngleCoeffs[4+n]) + clockAngleCoeffs[3+n];
    }

    else if (t > period*365/4 && t<= 2.0*period*365/4)
    {
        n = 5;
        coneAngleValue = coneAngleCoeffs[0+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
                + coneAngleCoeffs[1+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
                + coneAngleCoeffs[2+n]*(t-deltaT+coneAngleCoeffs[4+n]) + coneAngleCoeffs[3+n];

        clockAngleValue = clockAngleCoeffs[0+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
                + clockAngleCoeffs[1+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
                + clockAngleCoeffs[2+n]*(t-deltaT+clockAngleCoeffs[4+n]) + clockAngleCoeffs[3+n];
    }

    else if (t > 2.0*period*365/4 && t<= 3.0*period*365/4)
    {
        n = 10;
        coneAngleValue = coneAngleCoeffs[0+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
                + coneAngleCoeffs[1+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
                + coneAngleCoeffs[2+n]*(t-deltaT+coneAngleCoeffs[4+n]) + coneAngleCoeffs[3+n];

        clockAngleValue = clockAngleCoeffs[0+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
                + clockAngleCoeffs[1+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
                + clockAngleCoeffs[2+n]*(t-deltaT+clockAngleCoeffs[4+n]) + clockAngleCoeffs[3+n];
    }

    else if (t > 3.0*period*365/4 && t<= 4.0*period*365/4)
    {
        n = 15;
        coneAngleValue = coneAngleCoeffs[0+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
                + coneAngleCoeffs[1+n]*(t-deltaT+coneAngleCoeffs[4+n])*(t-deltaT+coneAngleCoeffs[4+n])
                + coneAngleCoeffs[2+n]*(t-deltaT+coneAngleCoeffs[4+n]) + coneAngleCoeffs[3+n];

        clockAngleValue = clockAngleCoeffs[0+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
                + clockAngleCoeffs[1+n]*(t-deltaT+clockAngleCoeffs[4+n])*(t-deltaT+clockAngleCoeffs[4+n])
                + clockAngleCoeffs[2+n]*(t-deltaT+clockAngleCoeffs[4+n]) + clockAngleCoeffs[3+n];
    }

    else
    {
        std::cout<<"error";
    }

    //angles[0] = coneAngleValue;
    //angles[1] = clockAngleValue;
    Eigen::Vector2d angles (coneAngleValue,clockAngleValue);

      return angles;
}
