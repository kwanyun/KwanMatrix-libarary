#pragma once

#ifndef MATRIXLIBRARY_API
#define MATRIXLIBRARY_API __declspec(dllexport)
#else
#define MATRIXLIBRARY_API __declspec(dllimport)
#endif


class MATRIXLIBRARY_API KwanMat {
private:
    float RecursiveDet(const float mat[25], const int n);
    void Cofactor(const float M[25], float t[16], int p, int q, int n);
    void ADJ(float M[25], float adj[25]);
public:
    float* matNums;
    unsigned int numRows, numCols;

    explicit KwanMat(const int rows, const int cols, float matrices[]);
    ~KwanMat();

    //uniform calculation of float
    void add(const float value);
    void sub(const float value);
    void mul(const float value);
    void div(const float value);

    //memberwise calculation
    KwanMat& operator+(const KwanMat& ref) const;
    KwanMat& operator-(const KwanMat& ref) const;
    KwanMat& operator*(const KwanMat& ref) const;

    //outer product
    KwanMat& MatMul(const KwanMat& ref);

    //transpose
    KwanMat& T();

    //for advanced calculations
    //Calcualte determinant with Laplace expansion for genenral size of square matrice for less than 6
    float DET();
    KwanMat& Inverse();

    //Deprecated
    int getRowSize();
    int getColumnSize();

};

MATRIXLIBRARY_API KwanMat & Identity(const unsigned int len);

MATRIXLIBRARY_API KwanMat & Zero(const unsigned int N, const unsigned int M);

MATRIXLIBRARY_API KwanMat & Rotation2D(const float theta);

MATRIXLIBRARY_API KwanMat & Rotation3D(const float theta,const char axis);

MATRIXLIBRARY_API bool IsSquare(const KwanMat & mat);

MATRIXLIBRARY_API bool IsSymmetric(const KwanMat & mat);

