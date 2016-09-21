#include <iostream>
#include <vector>

std::vector< std::vector< double > > SetBounds()
{

    unsigned int numberOfParameters = 24; // TO BE ADAPTED TO MY CASE, originally 2

    std::vector< std::vector< double > > bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );

    // Bounds for the optimized parameters, 2xn vector
    for (int i = 0;  i < numberOfParameters; i ++)
    {
        bounds[0][i] = -3.0;
        bounds[1][i] = 3.0;
    }

    return bounds;

}
