#ifndef SINGLEPERTURBEDSATELLITEPROPAGATOR_H
#define SINGLEPERTURBEDSATELLITEPROPAGATOR_H



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
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createBodies.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>

#include <Tudat/InputOutput/matrixTextFileReader.h>

#include <fstream>


// C++ Standard library
#include <iostream>

// External libraries: Eigen
#include <Eigen/Core>



double FFSSSimulation(tudat::basic_mathematics::Vector6d athosInitialStateInKeplerianElements,
                      tudat::basic_mathematics::Vector6d porthosInitialStateInKeplerianElements,
        tudat::basic_mathematics::Vector6d aramisInitialStateInKeplerianElements,
                      boost::function< Eigen::Vector2d( double ) >angleFunctionAthos,
                      boost::function< Eigen::Vector2d( double ) >angleFunctionPorthos);


#endif // SINGLEPERTURBEDSATELLITEPROPAGATOR_H
