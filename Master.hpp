#ifndef MASTERHEADERDEF
#define MASTERHEADERDEF
#include<vector>
#include<string>

class BoundaryConditions {
    private:
        bool mLhsBcIsDirichlet;
        bool mRhsBcIsDirichlet;
        bool mLhsBcIsNeumann;
        bool mRhsBcIsNeumann;
        double mLhsBcValue;
        double mRhsBcValue;
    public:
        friend class BvpOde;
        BoundaryConditions();
        void SetLhsDirichletBc(double lhsValue);
        void SetRhsDirichletBc(double rhsValue);
        void SetLhsNeumannBc (double lhsDerivValue);
        void SetRhsNeumannBc (double rhsDerivValue);
};

class Node
{
    public:
        double coordinate;
};

class Vector{
    private:
        double* mData; //data stored in a vector
        int mSize; //size of a vector
    public:
    Vector(int size);
    int GetSize() const;
    double operator[] (int i);
    double& operator() (int i);    
};

class SecondOrderOde{
    friend class BvpOde;
public:
    double mCoeffofUxx;
    double mCoeffofUx;
    double mCoeffofU;
    double (*mpRhsFunc) (double x);
    double mXmin;
    double mXmax;
    SecondOrderOde()=default;
    SecondOrderOde(double coeffUxx, double coeffUx, double coeffU,double (*righthandside)(double),double xMinimum,double xMaximum);    
};

class Matrix
{
    private:
        double** mData; 
        int mNumRows,mNumcols;
    public:
        Matrix(const Matrix& otherMatrix);
        Matrix(int numRows,int numCols);
        ~Matrix();
        int GetNumberOfRows() const;
        int GetNumberOfCols() const;
        double& operator()(int i,int j); // 1-based indexing
};
class LinearSystem
{
    private:
        int mSize; // Size of the linear system
        Matrix* mpA; // matrix for linear system
        Vector* mpb; // vector for linear system(rhs vector)
        /*Copy constructor is private. Only allow constructor that specifies matrix and vector */
        LinearSystem (const LinearSystem& otherLinearSystem);
    public:
        LinearSystem(const Matrix& A, const Vector& b);
        virtual ~LinearSystem(); // destructor frees memory
        void Mumps_input_gen();
};

class FiniteDifferenceGrid{ 
    public:
        std::vector<Node> mNodes;
        FiniteDifferenceGrid(std::vector<Node>::size_type numNodes, double xMin, double xMax);
};

class BvpOde{
    private:
        FiniteDifferenceGrid* mpGrid;
        SecondOrderOde* mpOde;
        BoundaryConditions* mpBconds;
        Vector* mpRhsVec;
        Matrix* mpLhsMat;
        LinearSystem* mpLinearSystem;
        std::string mFilename;
        void PopulateMatrix();
        void PopulateVector();
        void ApplyBoundaryConditions();
        int mNumNodes;
    public:
        BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs,int numNodes);
        ~BvpOde();
        void SetFilename(const std::string& name){
            mFilename= name;
        }
        void Solve();
};
extern "C"
{
    void Solver(int argc, char** argv);
}
#endif
