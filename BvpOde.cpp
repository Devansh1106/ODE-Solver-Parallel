#include <string>
#include<fstream>
#include<cassert>
#include<iostream>
#include"Master.hpp"
using namespace std;
BvpOde::BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs,int numNodes)
{
    mpOde = pOde;
    mpBconds = pBcs;

    mNumNodes = numNodes;
    mpGrid = new FiniteDifferenceGrid (mNumNodes, pOde->mXmin, pOde->mXmax);
    mpRhsVec = new Vector(mNumNodes);
    mpLhsMat = new Matrix(mNumNodes, mNumNodes);

    mpLinearSystem = NULL;
}
BvpOde::~BvpOde()
{
    //deletes memory allocated by constructor
    delete mpRhsVec;
    delete mpLhsMat;
    delete mpGrid;
    delete mpLinearSystem;   
}
void BvpOde::Solve()
{
    PopulateMatrix();
    PopulateVector();
    ApplyBoundaryConditions();
    mpLinearSystem = new LinearSystem (*mpLhsMat, *mpRhsVec);
    mpLinearSystem->Mumps_input_gen();
}

void BvpOde::PopulateMatrix()
{
    for (int i=1; i < mNumNodes-1; i++)
    {
        //xm, x and xp are x(i-1) ,x(i) and x(i+1)
        double xm = mpGrid->mNodes[i-1].coordinate; 
        double x = mpGrid->mNodes[i].coordinate; 
        double xp = mpGrid->mNodes[i+1].coordinate;
        double alpha = 2.0/(xp-xm)/(x-xm); 
        double beta = -2.0/(xp-x)/(x-xm); 
        double gamma = 2.0/(xp-xm)/(xp-x);
        (*mpLhsMat)(i+1,i) = (mpOde->mCoeffofUxx)*alpha - (mpOde->mCoeffofUx)/(xp-xm);
        (*mpLhsMat)(i+1,i+1) = (mpOde->mCoeffofUxx)*beta + (mpOde->mCoeffofU);
        (*mpLhsMat)(i+1,i+2) = (mpOde->mCoeffofUxx)*gamma + (mpOde->mCoeffofUx)/(xp-xm);
    }
}

void BvpOde::PopulateVector()
{
    for (int i=0; i < mNumNodes-1; i++)
    {
        double x = mpGrid->mNodes[i].coordinate;
        (*mpRhsVec)(i+1) = mpOde->mpRhsFunc(x);
    }
}

void BvpOde::ApplyBoundaryConditions()
{
    bool left_bc_applied = false;
    bool right_bc_applied = false;

    if (mpBconds->mLhsBcIsDirichlet)
    {
        (*mpLhsMat)(1,1) = 1.0;
        (*mpRhsVec)(1) = mpBconds->mLhsBcValue;
        left_bc_applied = true;
    }
    if (mpBconds->mRhsBcIsDirichlet)
    {
        (*mpLhsMat)(mNumNodes,mNumNodes) = 1.0;
        (*mpRhsVec)(mNumNodes) = mpBconds->mRhsBcValue;
        right_bc_applied = true;
    }
    if (mpBconds->mLhsBcIsNeumann)
    {
        assert(left_bc_applied == false);
        double h = mpGrid->mNodes[1].coordinate - mpGrid->mNodes[0].coordinate;
        (*mpLhsMat)(1,1) = -1.0/h;
        (*mpLhsMat)(1,2) = 1.0/h;
        (*mpRhsVec)(1) = mpBconds->mLhsBcValue;
        left_bc_applied = true;
    }
    if (mpBconds->mRhsBcIsNeumann)
    {
        assert(right_bc_applied == false);
        double h = mpGrid->mNodes[mNumNodes-1].coordinate - mpGrid->mNodes[mNumNodes-2].coordinate;
        (*mpLhsMat)(mNumNodes,mNumNodes-1) = -1.0/h;
        (*mpLhsMat)(mNumNodes,mNumNodes) = 1.0/h;
        (*mpRhsVec)(mNumNodes) = mpBconds->mRhsBcValue;
        right_bc_applied = true;
    }
    
    //Check that boundary conditions have been applied on both boundaries
    assert(right_bc_applied);
    assert(left_bc_applied);
}