#pragma once

#ifndef MATRIXLIBRARY_API
#define MATRIXLIBRARY_API __declspec(dllexport)
#else
#define MATRIXLIBRARY_API __declspec(dllimport)
#endif


namespace KMath {

    extern "C" class MATRIXLIBRARY_API Mat {
    private:
        
    public:
        float* matNums;
        unsigned int N, M;

        Mat(const int rows, const int cols, float matrices[]);
        ~Mat();

        //uniform calculation of float
        void add(const float value);
        void sub(const float value);
        void mul(const float value);
        void div(const float value);

        //outer product
        Mat& MatMul(Mat& ref);

        //transpose
        Mat& T();

        //memberwise calculation
        Mat& operator+(const Mat& ref);
        Mat& operator-(const Mat& ref);
        Mat& operator*(const Mat& ref);

        //Deprecated
        int getRowSize();
        int getColumnSize();

    };

    extern "C" MATRIXLIBRARY_API Mat & Identity(const unsigned int len);

    extern "C" MATRIXLIBRARY_API Mat & Zero(const unsigned int N, const unsigned int M);

    extern "C" MATRIXLIBRARY_API void CopyArr(float start[], float dest[], unsigned int size);

}