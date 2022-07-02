#pragma once

#ifndef MATRIXLIBRARY_API
#define MATRIXLIBRARY_API __declspec(dllexport)
#else
#define MATRIXLIBRARY_API __declspec(dllimport)
#endif



extern "C" class MATRIXLIBRARY_API KwanMat {
private:

public:
    float* matNums;
    unsigned int numRows, numCols;

    KwanMat(const int rows, const int cols, float matrices[]);
    ~KwanMat();

    //uniform calculation of float
    void add(const float value);
    void sub(const float value);
    void mul(const float value);
    void div(const float value);

    //outer product
    KwanMat& MatMul(KwanMat& ref);

    //transpose
    KwanMat& T();

    //memberwise calculation
    KwanMat& operator+(const KwanMat& ref);
    KwanMat& operator-(const KwanMat& ref);
    KwanMat& operator*(const KwanMat& ref);

    //Deprecated
    int getRowSize();
    int getColumnSize();

};

extern "C" MATRIXLIBRARY_API KwanMat & Identity(const unsigned int len);

extern "C" MATRIXLIBRARY_API KwanMat & Zero(const unsigned int N, const unsigned int M);

