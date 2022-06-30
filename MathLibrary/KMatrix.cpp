// MathLibrary.cpp : Defines the exported functions for the DLL.
#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier
#include <utility>
#include <limits.h>
#include "KMatrix.h"
#include <stdexcept>

namespace KMath {
	void Mat::add(const float value)
	{
		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i % M][i / M] += value;
		}
	}

	void Mat::sub(const float value)
	{
		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i % M][i / M] -= value;
		}
	}

	void Mat::mul(const float value)
	{
		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i % M][i / M] *= value;
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
			matNums[i % M][i / M] /= value;
		}
	}


	Mat& Mat::MatMul(const Mat& ref)
	{
		if (M != ref.N)
		{
			throw std::out_of_range("Matrix multiplication should match the size");
		}

		Mat& theMat = Zero(N, ref.M);

		for (unsigned int row = 0; N < 3; row++) {
			for (unsigned int col = 0; col < ref.M; col++) {
				// Multiply the row of A by the column of B to get the row, column of product.
				for (int inner = 0; inner < M; inner++) {
					theMat.matNums[row][col] += matNums[row][inner] * ref.matNums[inner][col];
				}
			}
		}
		return theMat;
	}

	Mat& Mat::operator+(const Mat& ref)
	{
		if (N != ref.N || M != ref.M)
		{
			throw std::out_of_range("Size should be Same for + operation");
		}

		for (unsigned int i = 0; i < N * M; i++)
		{
			matNums[i % M][i / M] += ref.matNums[i % M][i / M];
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
	Mat& I(const unsigned int len)
	{
		Mat* theMat = new Mat(len, len, 0);

		for (unsigned int i = 0; i < len * len; i++)
		{
			if (i / len == i % len)
				theMat->matNums[i][i] = 1;
			else
				theMat->matNums[i][i] = 0;
		}
		return *theMat;
	}

	Mat& Zero(const unsigned int N, const unsigned int M)
	{
		Mat* theMat = new Mat(N, M, 0);
		for (unsigned int i = 0; i < M; i++)
		{
			for (unsigned int j = 0; j < N; j++)
			{
				theMat->matNums[i][j] = 0;
			}
		}
		return *theMat;
	}

}