/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// External libraries: Boost
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

// Tudat library
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include </home/yeargh/tudatBundle/tudat/Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createBodies.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>

#include <fstream>

#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/penaltyfunction.h>
#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/singlePerturbedSatellitePropagator.h>


// C++ Standard library
#include <iostream>

// External libraries: Eigen
#include <Eigen/Core>

//#include "SatellitePropagatorExamples/applicationOutput.h"




double FFSSSimulation(tudat::basic_mathematics::Vector6d athosInitialStateInKeplerianElements,
                      tudat::basic_mathematics::Vector6d porthosInitialStateInKeplerianElements,
        tudat::basic_mathematics::Vector6d aramisInitialStateInKeplerianElements,
                      boost::function< Eigen::Vector2d( double ) >angleFunctionAthos,
                      boost::function< Eigen::Vector2d( double ) >angleFunctionPorthos)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace ephemerides;
    using namespace interpolators;
    //using namespace basic_astrodynamics;
    //using namespace input_output;



    //IMPORTING DATA


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT       //////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 365.0/2*tudat::physical_constants::JULIAN_DAY;
    const double fixedStepSize = 60.0;

    //const unsigned int numberOfTimeSteps = simulationEndEpoch/fixedStepSize;


    // Set number of satellites in constellation.
    const unsigned int numberOfSatellites = 2;
    const unsigned int sizeOfState = 6;


    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLES           /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create spacecraft objects.
    bodyMap[ "Athos" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Athos" ]->setConstantBodyMass( 3.0 );
    bodyMap[ "Porthos" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Porthos" ]->setConstantBodyMass( 3.0 );
    /*bodyMap[ "Aramis" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Aramis" ]->setConstantBodyMass( 3.0 );*/

    // SRP acceleration settings
    boost::shared_ptr< RadiationPressureInterfaceSettings > athosRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
            "Sun", 25.0, 2.0, boost::assign::list_of( "Earth" )( "Moon" ) );
    bodyMap[ "Athos" ]->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
            athosRadiationPressureSettings, "Athos", bodyMap ) );

   boost::shared_ptr< RadiationPressureInterfaceSettings > porthosRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
            "Sun", 25.0, 1.0, boost::assign::list_of( "Earth" )( "Moon" ) );
    bodyMap[ "Porthos" ]->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
            porthosRadiationPressureSettings, "Porthos", bodyMap ) );

    /*boost::shared_ptr< RadiationPressureInterfaceSettings > aramisRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
            "Sun", 4.0, 1.0, boost::assign::list_of( "Earth" )( "Moon" ) );
    bodyMap[ "Aramis" ]->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
            aramisRadiationPressureSettings, "Aramis", bodyMap ) );*/

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Definition of intial conditions
    Eigen::MatrixXd initialConditions( sizeOfState, numberOfSatellites );


    //THIS MIGHT NOT BE NECESSARY
    std::cout << bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Satellites initial conditions
    initialConditions.col(0) = convertKeplerianToCartesianElements(
                athosInitialStateInKeplerianElements, earthGravitationalParameter );


    initialConditions.col(1) = convertKeplerianToCartesianElements(
                porthosInitialStateInKeplerianElements, earthGravitationalParameter );


    initialConditions.col(2) = convertKeplerianToCartesianElements(
                aramisInitialStateInKeplerianElements, earthGravitationalParameter );

    // Definition of the system initial state
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( sizeOfState * numberOfSatellites );
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        systemInitialState.segment( i * 6, 6 ) = initialConditions.col( i );
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAthos;

    accelerationsOfAthos[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );

    accelerationsOfAthos[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAthos[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
   accelerationsOfAthos[ "Sun" ].push_back( boost::make_shared< FlatPlateRadiationPressureAccelerationSettings >(
                                             angleFunctionAthos ) );



    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfPorthos;

    accelerationsOfPorthos[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );

    accelerationsOfPorthos[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfPorthos[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfPorthos[ "Sun" ].push_back( boost::make_shared< FlatPlateRadiationPressureAccelerationSettings >(
                                              angleFunctionPorthos) );



    /*std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAramis;

    accelerationsOfAramis[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );

    accelerationsOfAramis[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAramis[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );*/
   // accelerationsOfAramis[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
         //                                         basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationMap[  "Athos" ] = accelerationsOfAthos;
    accelerationMap[  "Porthos" ] = accelerationsOfPorthos;
    //accelerationMap[  "Aramis" ] = accelerationsOfAramis;

    bodiesToPropagate.push_back( "Athos" );
    bodiesToPropagate.push_back( "Porthos" );
    //bodiesToPropagate.push_back( "Aramis" );


    centralBodies.push_back( "Earth" );
    centralBodies.push_back( "Earth" );
    //centralBodies.push_back( "Earth" );


    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

   dependentVariables.push_back(
               boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
                 "Athos", "Earth" ) );
  //  dependentVariables.push_back(
  //              boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
 //                   "Athos", "Porthos" ) );
//    dependentVariables.push_back(
 //               boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
//                    "Athos", "Aramis" ) );
  dependentVariables.push_back(
               boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
                   "Porthos", "Earth" ) );

 //   dependentVariables.push_back(
//                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
 //                   spherical_harmonic_gravity, "Athos", "Earth", 1 ) );
  //  dependentVariables.push_back(
 //               boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
 //                   spherical_harmonic_gravity, "Athos", "Earth", 0 ) );
    dependentVariables.push_back(
              boost::make_shared< SingleDependentVariableSaveSettings >(
                   total_acceleration_dependent_variable, "Athos" ) );

    dependentVariables.push_back(
             boost::make_shared< SingleDependentVariableSaveSettings >(
                   total_acceleration_dependent_variable, "Porthos" ) );

       dependentVariables.push_back(
                  boost::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable,
                      "Athos", "Sun" ) );

  // dependentVariables.push_back(
         //        boost::make_shared< SingleAccelerationDependentVariableSaveSettings
        //      >( single_plate_radiation_pressure,
           //        "Athos", "Sun",0 ) );

  // dependentVariables.push_back(
    //             boost::make_shared< SingleAccelerationDependentVariableSaveSettings
     //           >( cannon_ball_radiation_pressure,
       //             "Porthos", "Sun",0 ) );



    // Create acceleration models and propagation settings.

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              simulationEndEpoch, cowell,
              boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );



    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    /*std::map< double, Eigen::VectorXd > dependentVariableSoution =
            dynamicsSimulator.getDependentVariableHistory( );*/


    // Computes the value of the penalty function
    double pentaltyValue = PenaltyFunction(numericalSolution);
    double FitnessFunctionValue = 0; // Must be set proportional to the sum of the second time
                                    //derivatives of the attitude angles


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        OUTPUT WRITING              ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Write perturbed satellite propagation history to file.
    /*
    input_output::writeDataMapToTextFile( dependentVariableSoution,
                                          "dependentVariables1year.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( numericalSolution,
                                          "stateVectors1year.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );*/



    return (pentaltyValue + FitnessFunctionValue);
}

