#pragma once

#ifndef MATRIXLIBRARY_API
#define MATRIXLIBRARY_API __declspec(dllexport)
#else
#define MATRIXLIBRARY_API __declspec(dllimport)
#endif

namespace KMath {

    class MATRIXLIBRARY_API Mat {
    private:
        unsigned int N, M;
    public:
        float** matNums;

        Mat(const int rows, const int cols, void* matrices) : N(rows), M(cols), matNums((float**)matrices) {};
        //~Mat();
        //uniform calculation of float
        void add(const float value);
        void sub(const float value);
        void mul(const float value);
        void div(const float value);

        //outer product
        Mat& MatMul(const Mat& ref);

        //memberwise calculation
        Mat& operator+(const Mat& ref);
        //Mat& operator-(const Mat& ref);
        //Mat& operator*(const Mat& ref);

        int getRowSize();
        int getColumnSize();

    };

    extern "C" MATRIXLIBRARY_API Mat & I(const unsigned int len);

    extern "C" MATRIXLIBRARY_API Mat & Zero(const unsigned int N, const unsigned int M);
}