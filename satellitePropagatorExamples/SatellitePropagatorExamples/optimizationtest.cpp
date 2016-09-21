
#include <limits>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"



// External libraries: Eigen
#include <Eigen/Core>


// Simulation
#include <SatellitePropagatorExamples/singlePerturbedSatellitePropagator.cpp>

// Header
#include <SatellitePropagatorExamples/optimizationtest.h>

namespace pagmo {
namespace problem {

SolarSailFormationFlying::SolarSailFormationFlying(
    const std::vector< std::vector< double > > problemBounds ) :
    base( problemBounds[ 0 ], problemBounds[ 1 ], 0, 1 ),
    problemBounds_( problemBounds )
{ }

//! Clone method.
base_ptr SolarSailFormationFlying::clone( ) const {
        return base_ptr( new SolarSailFormationFlying( *this ) );
}

//! Descriptive name of the problem
std::string SolarSailFormationFlying::get_name() const {
    return "Formation flying solar sailing optimization problem";
}

//! Implementation of the objective function.
void SolarSailFormationFlying::objfun_impl( fitness_vector &f, const decision_vector &xv ) const
{

    // Set Keplerian elements for Athos.
    basic_mathematics::Vector6d athosInitialStateInKeplerianElements;
    athosInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    athosInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    athosInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    athosInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    athosInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    athosInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Set Keplerian elements for Porthos.
    basic_mathematics::Vector6d porthosInitialStateInKeplerianElements;
    porthosInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    porthosInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    porthosInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    porthosInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    porthosInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    porthosInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Set Keplerian elements for Aramis.
    basic_mathematics::Vector6d aramisInitialStateInKeplerianElements;
    aramisInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    aramisInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    aramisInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    aramisInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    aramisInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    aramisInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );



    double mu_E = physical_constants::EARTH_GRAVITATIONAL_PARAMETER;
    double period = std::sqrt(7500.0E3*7500.0E3*7500.0E3*4*mathematical_constants::PI*mathematical_constants::PI/mu_E);

    //Definition of the attitude angles function
    std::vector<std::vector<double>> coefficients = polynomialCoefficients(&xv,period);

    boost::function<std::vector<double> (double t, std::vector<std::vector<double>>& coeff, double period)> angleFunction;

    struct fun {
      std::vector<double> operator()(double t,std::vector<std::vector<double>>& coeff, double period) const
      {
        std::vector<double> coneAngleCoeffs, clockAngleCoeffs;
        std::vector<double> angles(2);

        double coneAngleValue;
        double clockAngleValue;
        int numOfRev,n;
        double deltaT;

        numOfRev = t/period;
        deltaT = numOfRev*period;



        coneAngleCoeffs = coeff[0];
        clockAnglCoeffs = coeff[1];

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

        angles[0] = coneAngleValue;
        angles[1] = clockAngleValue;
          return angles;
      }
      };

    angleFunction = fun();





    f[0] = FFSSSimulation(athosInitialStateInKeplerianElements,porthosInitialStateInKeplerianElements,
            aramisInitialStateInKeplerianElements, angleFunction );

}


} // namespace problem
} //namespace pagmo

BOOST_CLASS_EXPORT_IMPLEMENT( pagmo::problem::SolarSailFormationFlying );
