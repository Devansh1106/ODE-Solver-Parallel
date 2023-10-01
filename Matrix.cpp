#include"Master.hpp"
#include<cassert>
#include<cmath>

// Allocate memory for new matrix using double pointer
// and copy entries into the matrix
Matrix::Matrix (const Matrix& otherMatrix)
{
    mNumRows=otherMatrix.mNumRows;
    mNumcols=otherMatrix.mNumcols;          
    mData= new double *[mNumRows];         
    for (int i=0; i<mNumRows; i++)
    {
        mData[i]= new double [mNumcols];
    }
    for(int i=0; i<mNumRows; i++)
    {
        for(int j=0; j<mNumcols; j++)
        {
            mData[i][j]=otherMatrix.mData[i][j];
        }
    }
}

//Constructor for matrix of given row and columns
//Memory allocation for matrix and
//Initialization of entries to 0.0
Matrix::Matrix (int numRows, int numCols)
{
    assert(numRows > 0);
    assert(numCols > 0);
    mNumRows = numRows;
    mNumcols = numCols;
    mData = new double* [mNumRows];
    for (int i=0; i < mNumRows; i++)
    {
        mData[i] = new double [mNumcols];
    }
    for (int i=0; i < mNumRows; i++ )
    {
        for (int j=0; j<mNumcols; j++)
        {
            mData[i][j] = 0.0;
        }
    }
}

//overwritten destructor to free memory
Matrix::~Matrix()
{
    for (int i=0; i<mNumRows; i++)
    {
        delete [] mData[i];
    }
    delete [] mData;
}

//Method to get number of rows of matrix 
int Matrix::GetNumberOfRows() const
{
    return mNumRows;
}

//Method to get number of columns
int Matrix::GetNumberOfCols() const
{
    return mNumcols;
}

//Overloading round brackets
//Uses 1-based indexing
//and a check on validity of index
double& Matrix::operator()(int i, int j)
{
    assert(i > 0);
    assert(i < mNumRows+1);
    assert(j > 0);
    assert(j < mNumcols+1);
    return mData[i-1][j-1];
}