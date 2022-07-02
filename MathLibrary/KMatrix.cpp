// Kmatrhix.cpp : Defines the exported functions for the DLL.
// Originally made by Kwan Yun.
#include "pch.h"
#include "KMatrix.h"

namespace KMath {
	//initialize. get 1D array and
	Mat::Mat(const int rows, const int cols, float matrices[]) : N(rows), M(cols) 
	{
		matNums = new float[N * M];
		if (matrices == 0)
		{
			matNums = { 0 };
		}
		else
		{
			for (unsigned int i = 0; i < N * M; i++)
				matNums[i] = matrices[i];
		}
	};

	Mat::~Mat() 
	{
		delete[] matNums; 
	};


	void Mat::add(const float value)
	{
		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] += value;
		}
	}

	void Mat::sub(const float value)
	{
		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] -= value;
		}
	}

	void Mat::mul(const float value)
	{
		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] *= value;
		}
	}

	void Mat::div(const float value)
	{
		if (value == 0)
		{
			throw std::overflow_error("Divide by zero exception");
		}

		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] /= value;
		}
	}


	Mat& Mat::MatMul(Mat& ref)
	{
		if (M != ref.N)
		{
			throw std::out_of_range("Matrix multiplication should match the size");
		}

		Mat& theMat = Zero(N, ref.M);

		for (unsigned int row = 0; row < N; row++) {
			for (unsigned int col = 0; col < ref.M; col++) {
				// Multiply the row of A by the column of B to get the row, column of product.
				for (unsigned int inner = 0; inner < M; inner++) {
					theMat.matNums[row*ref.M+col] += matNums[row*M+inner] * ref.matNums[inner*ref.M+col];
				}
			}
		}
		return theMat;
	}
	
	Mat& Mat::T()
	{
		//copy temporally
		float* tempNums = new float[N*M];
		for (unsigned int i = 0; i < N * M; i++)
			tempNums[i] = matNums[i];
		for (unsigned int i = 0; i < N; i++)
		{
			for (unsigned int j = 0; j < M; j++)
			{
				matNums[i*M+j] = tempNums[j * M + i];
			}
		}
		delete[] tempNums;

		int tempN = N;
		N = M;
		M = tempN;

		return *this;
	}


	Mat& Mat::operator+(const Mat& ref)
	{
		if (N != ref.N || M != ref.M)
		{
			throw std::out_of_range("Size should be Same for + operation");
		}

		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] += ref.matNums[i];
		}
		return *this;
	}

	Mat& Mat::operator-(const Mat& ref)
	{
		if (N != ref.N || M != ref.M)
		{
			throw std::out_of_range("Size should be Same for + operation");
		}

		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] -= ref.matNums[i];
		}
		return *this;
	}

	Mat& Mat::operator*(const Mat& ref)
	{
		if (N != ref.N || M != ref.M)
		{
			throw std::out_of_range("Size should be Same for + operation");
		}

		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i] *= ref.matNums[i];
		}
		return *this;
	}


	int Mat::getRowSize()
	{
		return N;
	}

	int Mat::getColumnSize()
	{
		return M;
	}

	/// <summary>
	/// Functions
	/// </summary>
	Mat& Identity(const unsigned int len)
	{
		Mat* theMat = new Mat(len, len, 0);

		for (unsigned int i = 0; i < len * len; i++)
		{
			if (i / len == i % len)
				theMat->matNums[i] = 1;
		}
		return *theMat;
	}

	Mat& Zero(const unsigned int N, const unsigned int M)
	{
		Mat* theMat = new Mat(N, M, 0);
		return *theMat;
	}

	void CopyArr(float start[], float dest[],unsigned int size)
	{
		for (unsigned int i = 0; i < size; i++)
			dest[i] = start[i];
	}
}