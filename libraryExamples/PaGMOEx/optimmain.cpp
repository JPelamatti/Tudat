/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>
#include <cmath>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <Eigen/Core>

#include <pagmo/src/pagmo.h>
#include <pagmo/src/rng.h>

#include "/home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/Problems/optimizationtest.h"
#include "/home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/setbounds.h"


using namespace pagmo;
using boost::format;

//! Execute main
int main( )
{
    // Set the PRNG seed, such that results are reproducable
    int seed = 123456;
    pagmo::rng_generator::set_seed( seed );

    // If we have archipelagos, we also set the seed there
    // arch.set_seeds(sim_id);
    // Similarly set the seed for any other PRNG we might be using:
    // srand( seed );

    // We have two decision variables each with a lower and upper
    // bound, create a vector of vectors that will contain these.
    std::vector< std::vector< double > > bounds = SetBounds();

    // Define the problem
    problem::SolarSailFormationFlying prob( bounds );

    // Create a population (8 is minimum for jDE)
    population pop( prob, 8 );

    // Select the self-adaptive differential evolution algorithm.
    // One generation per evolution step.
    algorithm::jde algo( 1 );

    unsigned int i = 0;

    // For n generation optimise the population
    for( ; i < 25; i++ )
    {
        algo.evolve( pop );
        int c = pop.get_best_idx( );
        decision_vector cx = pop.get_individual( c ).cur_x;
        fitness_vector  cf = pop.get_individual( c ).cur_f;
        std::cout << "GEN=" << i << " ID=" << c << " DV=" << cf[ 0 ]
                  << "m/s DEP=" << cx[ 0 ] << "JD TOF=" << cx[ 1 ] << "d" << std::endl;
    }

    // NOT SURE ABOUT THE RETURN
    return 0;
}
