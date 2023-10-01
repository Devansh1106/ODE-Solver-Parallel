#include<iostream>
#include<cassert>
#include"Master.hpp"

//Defining Constructor
Vector::Vector(int size)
{
    assert(size > 0);
    mSize=size;
    mData= new double [mSize];
    for (int i = 0; i < mSize; i++)
    {
        mData[i]=0.0;
    }
}

//to access private variables
int Vector::GetSize() const {
    return mSize;
}

//Overloading square braackets
//uses 'zero based' indexing
double Vector::operator[](int i)
{   assert(i > -1);
    assert(i < mSize);
    return mData[i];
}

//Overloading round brackets
//uses 'one-based' indexing
double& Vector::operator()(int i)
{
    assert(i > 0);
    assert(i <= mSize);
    return mData[i-1];
}