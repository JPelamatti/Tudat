#include <iostream>

#include </home/yeargh/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/angleCoefficientsComputation.h>
#include <vector>



// Function that computes the missing polynomial coefficients using the continuity constraints
std::vector<std::vector<double>> polynomialCoefficients ( std::vector<double> xv, double deltaT)
{
    double a1,b1,c1,d1,e1;
    double a2,b2,c2,d2,e2;
    double a3,b3,c3,d3,e3;
    double a4,b4,c4,d4,e4;
    double ka1,kc1;
    double ka2,kc2;
    double ka3,kc3;
    double ka4,kc4;
    std::vector< std::vector< double > > coeffs( 2, std::vector< double >( 20, 0.0 ) );


    b1 = xv[0];
    d1 = xv[1];
    e1 = xv[2];
    e2 = xv[3];
    e3 = xv[4];
    e4 = xv[5];

    //CONE ANGLE COMPUTATION

    //TO BE CHECKED

    ka1 = 2*deltaT/3/(e1*e1-(e1+deltaT)*(e1+deltaT));
    ka2 = 2*deltaT/3/(e2*e2-(e2+deltaT)*(e2+deltaT));
    ka3 = 2*deltaT/3/(e3*e3-(e3+deltaT)*(e3+deltaT));
    ka4 = 2*deltaT/3/(e4*e4-(e4+deltaT)*(e4+deltaT));

    kc1 = 2/3*(e1*e1*e1-(e1+deltaT)*(e1+deltaT)*(e1+deltaT))/(e1*e1-(e1+deltaT)*(e1+deltaT))+(e1*e1-(e1+deltaT)*(e1+deltaT))/deltaT;
    kc2 = 2/3*(e2*e2*e2-(e2+deltaT)*(e2+deltaT)*(e2+deltaT))/(e2*e2-(e2+deltaT)*(e2+deltaT))+(e2*e2-(e2+deltaT)*(e2+deltaT))/deltaT;
    kc3 = 2/3*(e3*e3*e3-(e3+deltaT)*(e3+deltaT)*(e3+deltaT))/(e3*e3-(e3+deltaT)*(e3+deltaT))+(e3*e3-(e3+deltaT)*(e3+deltaT))/deltaT;
    kc4 = 2/3*(e4*e4*e4-(e4+deltaT)*(e4+deltaT)*(e4+deltaT))/(e4*e4-(e4+deltaT)*(e4+deltaT))+(e4*e4-(e4+deltaT)*(e4+deltaT))/deltaT;

    a1 = 2*b1*deltaT/3/(e1*e1-(e1+deltaT)*(e1*deltaT));
    c1 = (a1*(e1*e1*e1-(e1+deltaT)*(e1+deltaT)*(e1+deltaT))+b1*(e1*e1-(e1+deltaT)*(e1+deltaT)))/deltaT;

    b2 = b1*(3*ka1*(e1+deltaT)*(e1+deltaT)+2*(e1+deltaT)+kc1)/(3*ka2*e2*e2+2*e2+kc2);
    a2 = 2*b2*deltaT/3/(e2*e2-(e2+deltaT)*(e2*deltaT));
    c2 = (a2*(e2*e2*e2-(e2+deltaT)*(e2+deltaT)*(e2+deltaT))+b2*(e2*e2-(e2+deltaT)*(e2+deltaT)))/deltaT;

    b3 = b2*(3*ka2*(e2+deltaT)*(e2+deltaT)+2*(e2+deltaT)+kc2)/(3*ka3*e3*e3+2*e3+kc3);
    a3 = 2*b3*deltaT/3/(e3*e1-(e3+deltaT)*(e3*deltaT));
    c3 = (a3*(e3*e3*e3-(e3+deltaT)*(e3+deltaT)*(e3+deltaT))+b3*(e3*e3-(e3+deltaT)*(e3+deltaT)))/deltaT;

    b4 = b3*(3*ka3*(e3+deltaT)*(e3+deltaT)+2*(e3+deltaT)+kc3)/(3*ka4*e4*e4+2*e4+kc4);
    a4 = 2*b4*deltaT/3/(e4*e4-(e4+deltaT)*(e4*deltaT));
    c4 = (a4*(e4*e4*e4-(e4+deltaT)*(e4+deltaT)*(e4+deltaT))+b4*(e4*e4-(e4+deltaT)*(e4+deltaT)))/deltaT;

    d2 = a1*(e1+deltaT)*(e1+deltaT)*(e1+deltaT)+b1*(e1+deltaT)*(e1+deltaT)+c1*(e1+deltaT)+d1-a2*e2*e2*e2-b2*e2*e2-c2*e2;
    d3 = a2*(e2+deltaT)*(e2+deltaT)*(e2+deltaT)+b2*(e2+deltaT)*(e2+deltaT)+c2*(e2+deltaT)+d2-a3*e3*e3*e3-b3*e3*e3-c3*e3;
    d4 = a3*(e3+deltaT)*(e3+deltaT)*(e3+deltaT)+b3*(e3+deltaT)*(e3+deltaT)+c3*(e3+deltaT)+d3-a4*e4*e4*e4-b4*e4*e4-c4*e4;

    coeffs[0][0] = a1;
    coeffs[0][1] = b1;
    coeffs[0][2] = c1;
    coeffs[0][3] = d1;
    coeffs[0][4] = e1;
    coeffs[0][5] = a2;
    coeffs[0][6] = b2;
    coeffs[0][7] = c2;
    coeffs[0][8] = d2;
    coeffs[0][9] = e2;
    coeffs[0][10] = a3;
    coeffs[0][11] = b3;
    coeffs[0][12] = c3;
    coeffs[0][13] = d3;
    coeffs[0][14] = e3;
    coeffs[0][15] = a4;
    coeffs[0][16] = b4;
    coeffs[0][17] = c4;
    coeffs[0][18] = d4;
    coeffs[0][19] = e4;




   // CLOCK ANGLE COMPUTATION

    b1 = xv[6];
    d1 = xv[7];
    e1 = xv[8];
    e2 = xv[9];
    e3 = xv[10];
    e4 = xv[11];


    //TO BE CHECKED

    ka1 = 2*deltaT/3/(e1*e1-(e1+deltaT)*(e1+deltaT));
    ka2 = 2*deltaT/3/(e2*e2-(e2+deltaT)*(e2+deltaT));
    ka3 = 2*deltaT/3/(e3*e3-(e3+deltaT)*(e3+deltaT));
    ka1 = 2*deltaT/3/(e4*e4-(e4+deltaT)*(e4+deltaT));

    kc1 = 2/3*(e1*e1*e1-(e1+deltaT)*(e1+deltaT)*(e1+deltaT))/(e1*e1-(e1+deltaT)*(e1+deltaT))+(e1*e1-(e1+deltaT)*(e1+deltaT))/deltaT;
    kc2 = 2/3*(e2*e2*e2-(e2+deltaT)*(e2+deltaT)*(e2+deltaT))/(e2*e2-(e2+deltaT)*(e2+deltaT))+(e2*e2-(e2+deltaT)*(e2+deltaT))/deltaT;
    kc3 = 2/3*(e3*e3*e3-(e3+deltaT)*(e3+deltaT)*(e3+deltaT))/(e3*e3-(e3+deltaT)*(e3+deltaT))+(e3*e3-(e3+deltaT)*(e3+deltaT))/deltaT;
    kc1 = 2/3*(e4*e4*e4-(e4+deltaT)*(e4+deltaT)*(e4+deltaT))/(e4*e4-(e4+deltaT)*(e4+deltaT))+(e4*e4-(e4+deltaT)*(e4+deltaT))/deltaT;

    a1 = 2*b1*deltaT/3/(e1*e1-(e1+deltaT)*(e1*deltaT));
    c1 = (a1*(e1*e1*e1-(e1+deltaT)*(e1+deltaT)*(e1+deltaT))+b1*(e1*e1-(e1+deltaT)*(e1+deltaT)))/deltaT;

    b2 = b1*(3*ka1*(e1+deltaT)*(e1+deltaT)+2*(e1+deltaT)+kc1)/(3*ka2*e2*e2+2*e2+kc2);
    a2 = 2*b2*deltaT/3/(e2*e2-(e2+deltaT)*(e2*deltaT));
    c2 = (a2*(e2*e2*e2-(e2+deltaT)*(e2+deltaT)*(e2+deltaT))+b2*(e2*e2-(e2+deltaT)*(e2+deltaT)))/deltaT;

    b3 = b2*(3*ka2*(e2+deltaT)*(e2+deltaT)+2*(e2+deltaT)+kc2)/(3*ka3*e3*e3+2*e3+kc3);
    a3 = 2*b3*deltaT/3/(e3*e1-(e3+deltaT)*(e3*deltaT));
    c3 = (a3*(e3*e3*e3-(e3+deltaT)*(e3+deltaT)*(e3+deltaT))+b3*(e3*e3-(e3+deltaT)*(e3+deltaT)))/deltaT;

    b4 = b3*(3*ka3*(e3+deltaT)*(e3+deltaT)+2*(e3+deltaT)+kc3)/(3*ka4*e4*e4+2*e4+kc4);
    a4 = 2*b4*deltaT/3/(e4*e4-(e4+deltaT)*(e4*deltaT));
    c4 = (a4*(e4*e4*e4-(e4+deltaT)*(e4+deltaT)*(e4+deltaT))+b4*(e4*e4-(e4+deltaT)*(e4+deltaT)))/deltaT;

    d2 = a1*(e1+deltaT)*(e1+deltaT)*(e1+deltaT)+b1*(e1+deltaT)*(e1+deltaT)+c1*(e1+deltaT)+d1-a2*e2*e2*e2-b2*e2*e2-c2*e2;
    d3 = a2*(e2+deltaT)*(e2+deltaT)*(e2+deltaT)+b2*(e2+deltaT)*(e2+deltaT)+c2*(e2+deltaT)+d2-a3*e3*e3*e3-b3*e3*e3-c3*e3;
    d4 = a3*(e3+deltaT)*(e3+deltaT)*(e3+deltaT)+b3*(e3+deltaT)*(e3+deltaT)+c3*(e3+deltaT)+d3-a4*e4*e4*e4-b4*e4*e4-c4*e4;

    coeffs[1][0] = a1;
    coeffs[1][1] = b1;
    coeffs[1][2] = c1;
    coeffs[1][3] = d1;
    coeffs[1][4] = e1;
    coeffs[1][5] = a2;
    coeffs[1][6] = b2;
    coeffs[1][7] = c2;
    coeffs[1][8] = d2;
    coeffs[1][9] = e2;
    coeffs[1][10] = a3;
    coeffs[1][11] = b3;
    coeffs[1][12] = c3;
    coeffs[1][13] = d3;
    coeffs[1][14] = e3;
    coeffs[1][15] = a4;
    coeffs[1][16] = b4;
    coeffs[1][17] = c4;
    coeffs[1][18] = d4;
    coeffs[1][19] = e4;

return coeffs;
}
