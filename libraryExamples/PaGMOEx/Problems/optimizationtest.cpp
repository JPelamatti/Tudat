
#include <limits>
#include <iostream>


#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"


#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"



// External libraries: Eigen
#include <Eigen/Core>

// Header
#include "/home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/Problems/optimizationtest.h"

// Simulation
#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/singlePerturbedSatellitePropagator.h>

// Additional code
#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/angleCoefficientsComputation.h>
#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/computeanglesvalues.h>


namespace pagmo {
namespace problem {

SolarSailFormationFlying::SolarSailFormationFlying(
    const std::vector< std::vector< double > > problemBounds ) :
    base( problemBounds[ 0 ], problemBounds[ 1 ], 0, 1 ),
    problemBounds_( problemBounds )
{ }

//! Clone method.
base_ptr SolarSailFormationFlying::clone( ) const
{
        return base_ptr( new SolarSailFormationFlying( *this ) );
}

//! Descriptive name of the problem
std::string SolarSailFormationFlying::get_name() const {
    return "Formation flying solar sailing optimization problem";
}


//! Implementation of the objective function.
void SolarSailFormationFlying::objfun_impl( fitness_vector &f, const decision_vector &xv ) const
{

    unsigned int numberOfParameters = 24;
    std::vector<double> given_coeffs(numberOfParameters);

    for (int i = 0; i < numberOfParameters; i ++)
    {
        given_coeffs[i] = xv[i];
    }

    // Set Keplerian elements for Athos.
    tudat::basic_mathematics::Vector6d athosInitialStateInKeplerianElements;
    athosInitialStateInKeplerianElements( tudat::orbital_element_conversions::semiMajorAxisIndex ) = 41732.0E3;
    athosInitialStateInKeplerianElements( tudat::orbital_element_conversions::eccentricityIndex ) = 0;
    athosInitialStateInKeplerianElements( tudat::orbital_element_conversions::inclinationIndex ) = tudat::unit_conversions::convertDegreesToRadians( 0.1679 );
    athosInitialStateInKeplerianElements( tudat::orbital_element_conversions::argumentOfPeriapsisIndex )
            = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
    athosInitialStateInKeplerianElements( tudat::orbital_element_conversions::longitudeOfAscendingNodeIndex )
            = tudat::unit_conversions::convertDegreesToRadians( 49.1696);
    athosInitialStateInKeplerianElements( tudat::orbital_element_conversions::trueAnomalyIndex ) = tudat::unit_conversions::convertDegreesToRadians( 265.0217 );

    // Set Keplerian elements for Porthos.
    tudat::basic_mathematics::Vector6d porthosInitialStateInKeplerianElements;
    porthosInitialStateInKeplerianElements( tudat::orbital_element_conversions::semiMajorAxisIndex ) = 41732.0E3;
    porthosInitialStateInKeplerianElements( tudat::orbital_element_conversions::eccentricityIndex ) = 0;
    porthosInitialStateInKeplerianElements( tudat::orbital_element_conversions::inclinationIndex ) = tudat::unit_conversions::convertDegreesToRadians( 0.1839  );
    porthosInitialStateInKeplerianElements( tudat::orbital_element_conversions::argumentOfPeriapsisIndex )
            = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
    porthosInitialStateInKeplerianElements( tudat::orbital_element_conversions::longitudeOfAscendingNodeIndex )
            = tudat::unit_conversions::convertDegreesToRadians( 340.4865 );
    porthosInitialStateInKeplerianElements( tudat::orbital_element_conversions::trueAnomalyIndex ) = tudat::unit_conversions::convertDegreesToRadians( 334.2449 );


    // Set Keplerian elements for Aramis. NOT USED FOR NOW
    tudat::basic_mathematics::Vector6d aramisInitialStateInKeplerianElements;
    aramisInitialStateInKeplerianElements( tudat::orbital_element_conversions::semiMajorAxisIndex ) = 7500.0E3;
    aramisInitialStateInKeplerianElements( tudat::orbital_element_conversions::eccentricityIndex ) = 0.1;
    aramisInitialStateInKeplerianElements( tudat::orbital_element_conversions::inclinationIndex ) = tudat::unit_conversions::convertDegreesToRadians( 85.3 );
    aramisInitialStateInKeplerianElements( tudat::orbital_element_conversions::argumentOfPeriapsisIndex )
            = tudat::unit_conversions::convertDegreesToRadians( 235.7 );
    aramisInitialStateInKeplerianElements( tudat::orbital_element_conversions::longitudeOfAscendingNodeIndex )
            = tudat::unit_conversions::convertDegreesToRadians( 23.4 );
    aramisInitialStateInKeplerianElements( tudat::orbital_element_conversions::trueAnomalyIndex ) = tudat::unit_conversions::convertDegreesToRadians( 139.87 );



    double mu_E = tudat::celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER;
    double period = std::sqrt(7500.0E3*7500.0E3*7500.0E3*4*tudat::mathematical_constants::PI*tudat::mathematical_constants::PI/mu_E);

    //Definition of the attitude angles function
    std::vector<double> athosGivenCoeffs ;
    std::vector<double> porthosGivenCoeffs;

    for(int k = 0; k < given_coeffs.size()/2; k++)
    {
        athosGivenCoeffs.push_back(given_coeffs[k]);
        porthosGivenCoeffs.push_back(given_coeffs[k+given_coeffs.size()/2]);
    }

    std::vector<std::vector<double>> coefficientsAthos = polynomialCoefficients( athosGivenCoeffs,period);
    std::vector<std::vector<double>> coefficientsPorthos = polynomialCoefficients( porthosGivenCoeffs,period);

    //boost::function<std::vector<double> (double t)> angleFunction;
    boost::function< Eigen::Vector2d( double t) >angleFunctionAthos;
    boost::function< Eigen::Vector2d( double t) >angleFunctionPorthos;


    // in theory chaning the output type of compute angles should be enough
    angleFunctionAthos = boost::bind( &computeAngles, _1, coefficientsAthos, period );
    angleFunctionAthos = boost::bind( &computeAngles, _1, coefficientsPorthos, period );



    // Computation of the fitness function value
    f[0] = FFSSSimulation(athosInitialStateInKeplerianElements,porthosInitialStateInKeplerianElements,
            aramisInitialStateInKeplerianElements, angleFunctionAthos, angleFunctionPorthos );

}


} // namespace problem
} //namespace pagmo

BOOST_CLASS_EXPORT_IMPLEMENT( pagmo::problem::SolarSailFormationFlying );
