#include"Master.hpp"
#include<iostream>
SecondOrderOde::SecondOrderOde(double coeffUxx, double coeffUx, double coeffU,double (*righthandside)(double),double xMinimum,double xMaximum)
{
    mCoeffofUxx=coeffUxx;
    mCoeffofUx=coeffUx;
    mCoeffofU=coeffU;
    mpRhsFunc=righthandside;
    mXmax=xMaximum;
    mXmin=xMinimum;        
}