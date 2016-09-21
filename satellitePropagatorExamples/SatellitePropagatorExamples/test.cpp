#include <Eigen/Core>
#include <vector>


int main()
{

    std::map< double, std::vector< double> > stateVectors;
    std::vector<double> vec1 = {1,2};
    std::vector<double> vec2 = {3,4};
    std::vector<double> vec3 = {5,6};


    stateVectors[1.1].push_back(vec1);
    stateVectors[3,2].push_back(vec2);
    stateVectors[5,7].push_back(vec3);


    for (i = stateVectors.begin(); i != stateVectors.end(); i++)
    {
        std::cout << stateVectors[i];
    }

}
