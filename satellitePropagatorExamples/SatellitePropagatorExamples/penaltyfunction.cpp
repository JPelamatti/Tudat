#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Eigen/Core>
#include <cmath>


double PenaltyFunction(std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > stateVectors)
{

    Eigen::Matrix< double, Eigen::Dynamic, 3 > positionAthos;
    Eigen::Matrix< double, Eigen::Dynamic, 3 > positionPorthos;
    Eigen::Matrix< double, Eigen::Dynamic, 3 > positionAramis;

    Eigen::VectorXd <double> xFC,yFC,zFC,rFC;
    Eigen::VectorXd <double> inclinationFC;
    Eigen::VectorXd <double> dAthosPorthos;
    Eigen::VectorXd <double> dAthosAramis;
    Eigen::VectorXd <double> dPorthosAramis;

    double a,b,c,d;
    double inclinationFCMax,dAthosPorthosMax,dAthosAramisMax,dPorthosAramisMax;
    double dAthosPorthosNom, dAthosAramisNom, dPorthosAramisNom;


    //Penalty function weights, to be tuned!!!
    a = 1;
    b = 1;
    c = 1;
    d = 1;

    //Nominal distances, to be tuned (or make a function to define them f(t)
    dAthosPorthosNom = 1;
    dAthosAramisNom = 1;
    dPorthosAramisNom = 1;

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 >>::iterator i;
    int k;

    for(i = stateVectors.begin(), k = 0; i != stateVectors.end(); i = i++, k++)
    {
        timeVector = i->first;
        vectors = i->second;
        positionAthos.row(k) = vectors[i].segment(0,3);
        positionPorthos.row(k) = vectors[i].segment(6,3);
        positionAramis.row(k) = vectors[i].segment(12,3);
    }

    for(k = 0; k < timeVector.size(); k++ )
    {
    xFC(k) = (positionAthos(k,0)+positionPorthos(k,0)+positionAramis(k,0))/3;
    yFC(k) = (positionAthos(k,1)+positionPorthos(k,1)+positionAramis(k,1))/3;
    zFC(k) = (positionAthos(k,2)+positionPorthos(k,2)+positionAramis(k,2))/3;

    rFC(k) = cmath::sqrt(xFc(k)*xFC(k)+yFC*yFC(k)+zFC(k)*zFC(k));

    inclinationFC(k) = cmath::asin(zFC(k)/rFC(k));

    dAthosPorthos(k) = cmath::sqrt((positionAthos(k,0)-positionPorthos(k,0))*(positionAthos(k,0)-positionPorthos(k,0))
                                   +(positionAthos(k,1)-positionPorthos(k,1))*(positionAthos(k,1)-positionPorthos(k,1))
                                   +(positionAthos(k,2)-positionPorthos(k,2))*(positionAthos(k,2)-positionPorthos(k,2)));

    dAthosAramis(k) = cmath::sqrt((positionAthos(k,0)-positionAramis(k,0))*(positionAthos(k,0)-positionAramis(k,0))
                                   +(positionAthos(k,1)-positionAramis(k,1))*(positionAthos(k,1)-positionAramis(k,1))
                                   +(positionAthos(k,2)-positionAramis(k,2))*(positionAthos(k,2)-positionAramis(k,2)));

    dPorthosAramis(k) = cmath::sqrt((positionAramis(k,0)-positionPorthos(k,0))*(positionAramis(k,0)-positionPorthos(k,0))
                                   +(positionAramis(k,1)-positionPorthos(k,1))*(positionAramis(k,1)-positionPorthos(k,1))
                                   +(positionAramis(k,2)-positionPorthos(k,2))*(positionAramis(k,2)-positionPorthos(k,2)));
    }

    inclinationFCMax = *max_element(std::begin(inclinationFC), std::end(inclinationFC));
    dAthosPorthosMax = *max_element(std::begin(dAthosPorthos), std::end(dAthosPorthos));
    dAthosAramisMax = *max_element(std::begin(dAthosAramis), std::end(dAthosAramis));
    dPorthosAramisMax = *max_element(std::begin(dPorthosAramis), std::end(dPorthosAramis));

    pentaltyFunValue = a*std::abs(inclinationFCMax) + b*std::abs(dAthosPorthosMax-dAthosPorthosNom)
                        + c*std::abs(dAthosAramisMax-dAthosAramisNom) + d*std::abs(dPorthosAramisMax-dPorthosAramisNom);

    return pentaltyFunValue;

}
