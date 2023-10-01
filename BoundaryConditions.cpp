#include"Master.hpp"
#include<cassert>
BoundaryConditions::BoundaryConditions()
{
    mLhsBcIsDirichlet = false;
    mRhsBcIsDirichlet = false;
    mLhsBcIsNeumann = false;
    mRhsBcIsNeumann = false;
}
void BoundaryConditions::SetLhsDirichletBc(double lhsvalue)
{
    assert(!mLhsBcIsNeumann);
    mLhsBcIsDirichlet = true;
    mLhsBcValue = lhsvalue;
}
void BoundaryConditions::SetRhsDirichletBc(double rhsvalue)
{
    assert(!mRhsBcIsNeumann);
    mRhsBcIsDirichlet = true;
    mRhsBcValue = rhsvalue;
}
void BoundaryConditions::SetLhsNeumannBc(double lhsDerivValue)
{
    assert(!mLhsBcIsDirichlet);
    mLhsBcIsNeumann = true;
    mLhsBcValue = lhsDerivValue;
}
void BoundaryConditions::SetRhsNeumannBc(double rhsDerivValue)
{
    assert(!mRhsBcIsDirichlet);
    mRhsBcIsNeumann = true;
    mRhsBcValue = rhsDerivValue;
}