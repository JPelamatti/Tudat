#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Eigen/Core>
#include <cmath>

#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/penaltyfunction.h>

// Function that computes the value of the optimization penalty function given the state vector history

double PenaltyFunction(std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > stateVectors)
{

    Eigen::Matrix< double, Eigen::Dynamic, 3 > positionAthos;
    Eigen::Matrix< double, Eigen::Dynamic, 3 > positionPorthos;
    Eigen::Matrix< double, Eigen::Dynamic, 3 > positionAramis;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > vectors;

    std::vector<double> xFC,yFC,zFC,rFC;
    std::vector<double> inclinationFC;
    std::vector<double>  dAthosPorthos,dAthosAramis,dPorthosAramis;
    std::vector<double> timeVector;

    double a,b,c,d;
    double inclinationFCMax,dAthosPorthosSum,dAthosAramisMax,dPorthosAramisMax;
    double dAthosPorthosNom, dAthosAramisNom, dPorthosAramisNom;


    //Penalty function weights, to be tuned!!!
    a = 360/2/3.14; // approx rad to degrees conversion
    b = 1/1e5;
    c = 1;
    d = 1;

    //Nominal distances, to be tuned (or make a function to define them f(t)
    dAthosPorthosNom = 4e5; // [m]
    dAthosAramisNom = 1;
    dPorthosAramisNom = 1;


    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 >>::iterator i;
    unsigned int k;

    for(i = stateVectors.begin(), k = 0; i != stateVectors.end(); i = i++, k++)
    {
        timeVector.push_back(i->first);
        vectors = i->second;
        positionAthos.row(k) = vectors.block(0,0,3,1);
        positionPorthos.row(k) = vectors.block(6,0,3,1);
        //positionAramis.row(k) = vectors.block(12,0,3,1);
    }

    for(k = 0; k < timeVector.size(); k++ )
    {
    xFC.push_back((positionAthos(k,0)+positionPorthos(k,0))/2);
    yFC.push_back((positionAthos(k,1)+positionPorthos(k,1))/2);
    zFC.push_back((positionAthos(k,2)+positionPorthos(k,2))/2);

    rFC.push_back(std::sqrt(xFC[k]*xFC[k]+yFC[k]*yFC[k]+zFC[k]*zFC[k]));

    inclinationFC.push_back(std::asin(zFC[k]/rFC[k]));

    dAthosPorthos.push_back(std::sqrt((positionAthos(k,0)-positionPorthos(k,0))*(positionAthos(k,0)-positionPorthos(k,0))
                                   +(positionAthos(k,1)-positionPorthos(k,1))*(positionAthos(k,1)-positionPorthos(k,1))
                                   +(positionAthos(k,2)-positionPorthos(k,2))*(positionAthos(k,2)-positionPorthos(k,2))));

    /*dAthosAramis.push_back(((positionAthos(k,0)-positionAramis(k,0))*(positionAthos(k,0)-positionAramis(k,0))
                                   +(positionAthos(k,1)-positionAramis(k,1))*(positionAthos(k,1)-positionAramis(k,1))
                                   +(positionAthos(k,2)-positionAramis(k,2))*(positionAthos(k,2)-positionAramis(k,2))));

    dPorthosAramis.push_back(((positionAramis(k,0)-positionPorthos(k,0))*(positionAramis(k,0)-positionPorthos(k,0))
                                   +(positionAramis(k,1)-positionPorthos(k,1))*(positionAramis(k,1)-positionPorthos(k,1))
                                   +(positionAramis(k,2)-positionPorthos(k,2))*(positionAramis(k,2)-positionPorthos(k,2))));*/
    }

    inclinationFCMax = *max_element(std::begin(inclinationFC), std::end(inclinationFC));

    for( k = 0, dAthosPorthosSum = 0; k < dAthosPorthos.size(); k++)
    {
        dAthosPorthosSum = dAthosPorthosSum + std::abs(dAthosPorthos[k] - dAthosPorthosNom);
    }
    //dAthosAramisMax = *max_element(std::begin(dAthosAramis), std::end(dAthosAramis));
    //dPorthosAramisMax = *max_element(std::begin(dPorthosAramis), std::end(dPorthosAramis));



    return (a*std::abs(inclinationFCMax) + b*dAthosPorthosSum);
}
