#include"Master.hpp"
#include<cassert>
#include<fstream>
#include<iostream>

using namespace std;

// copy matrix and vector so that by gaussian elimination original matrix and vector remains unchanged
LinearSystem::LinearSystem(const Matrix& A, const Vector& b )
{
    //Check for size compatibility of matrix and vector
    int local_size= A.GetNumberOfRows();
    assert(A.GetNumberOfCols() == local_size);
    assert(b.GetSize() == local_size);

    //Set variables for linear system
    mSize=local_size;
    mpA = new Matrix(A);
    mpb = new Vector(b);
}

// Destructor frees memory
LinearSystem::~LinearSystem()
{
    delete mpA;
    delete mpb;
}

//Solve the linear system using Gaussian elimination
//This will change the matrix mpA
void LinearSystem::Mumps_input_gen()
{
    //Introducing pointer for easy access to Matrix and vector
    Matrix rA = *mpA;
    Vector rb = *mpb;

    //Printing the matrix and rhs in form of irn and jcn form
    std::ofstream file_1("irn.txt");
    std::ofstream file_2("jcn.txt");
    std::ofstream file_3("rhs.txt");
    std::ofstream file_4("a.txt");
    //std::cout<<"HI\n";
    if (!file_1.is_open() && !file_2.is_open() && !file_3.is_open() && !file_4.is_open())
    {
        std::cout<<"File cannot be opened while creating irn.txt etc.!";
    }
    else
    {
        assert(file_1.is_open());
        assert(file_2.is_open());
        assert(file_3.is_open());
        assert(file_4.is_open());

        int N = rA.GetNumberOfRows();
        // std::cout<<"\n"<<N<<"\n";
        for(int i = 1; i < N+1; i++)
        {
            for(int j = 1; j < N+1; j++)
            {
                if(rA(i,j) != 0)
                {
                    file_1<<i<<"\n";
                    file_2<<j<<"\n";
                    double y = rA(i,j);
                    file_4<<y<<"\n";
                }
            }
            double z = rb[i-1];
            file_3<<z<<"\n";
        }
        //file_1.flush();
        file_1.close();

        //file_2.flush();
        file_2.close();

        //file_3.flush();
        file_3.close();

        //file_4.flush();
        file_4.close();
    }
    //std::cout<<"Ouptput done!\n";
    // //forward sweep of gaussian elimination (conversion to upper triangular form)
    // for (int k=0; k < mSize-1; k++)
    // {
    //     //see if pivoting is necessary
    //     double max=0.0;
    //     int row = -1;
    //     for(int i = k; i < mSize; i++)
    //     {
    //         if (fabs(rA(i+1,k+1)) > max)
    //         {
    //             row = i;
    //         }
    //     }
    //     assert (row > 0);

    //     // pivot if necessary
    //     if (row != k)
    //     {
    //         //swap matrix k+1 with row+1
    //         for (int i=0; i < mSize; i++)
    //         {
    //             double temp = rA(k+1,i+1);
    //             rA(k+1,i+1) = rA(row+1,i+1);
    //             rA(row+1,i+1) = temp;
    //         }
    //         //swap vector enteries k+1 with row+1
    //         double temp = rb(k+1);
    //         rb(k+1) = rb(row+1);
    //         rb(row+1) = temp;
    //     }

    //     //create 0 in the lower part of column k
    //     for (int i=k+1; i < mSize; i++)
    //     {
    //         m(i+1) = rA(i+1,k+1)/rA(k+1,k+1);
    //         for (int j=k; j < mSize; j++)
    //         {
    //             rA(i+1,j+1) -= rA(k+1,j+1)*m(i+1);
    //         }
    //         rb(i+1) -= rb(k+1)*m(i+1);
    //     }
    // }

    // //back substitution
    // for (int i=mSize-1; i > -1; i--)
    // {
    //     solution(i+1) = rb(i+1);
    //     for (int j=i+1; j < mSize; j++ )
    //     {
    //         solution(i+1) -= rA(i+1,j+1)*solution(j+1);
    //     }
    //     solution (i+1) /= rA(i+1,i+1);
    // }
    // return solution;
}