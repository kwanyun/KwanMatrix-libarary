// Kmatrhix.cpp : Defines the exported functions for the DLL.
// Originally made by Kwan Yun.
#include "pch.h"
#include "KMatrix.h"


//initialize. get 1D array and
KwanMat::KwanMat(const int rows, const int cols, float matrices[]) : numRows(rows), numCols(cols) 
{
	matNums = new float[numRows * numCols];
	if (matrices == 0)
	{
		for (unsigned int i = 0; i < numRows * numCols; i++)
			matNums[i] = 0;
	}
	else
	{
		for (unsigned int i = 0; i < numRows * numCols; i++)
			matNums[i] = matrices[i];
	}
};

KwanMat::~KwanMat() 
{
	delete[] matNums; 
};


void KwanMat::add(const float value)
{
	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] += value;
	}
}

void KwanMat::sub(const float value)
{
	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] -= value;
	}
}

void KwanMat::mul(const float value)
{
	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] *= value;
	}
}

void KwanMat::div(const float value)
{
	if (value == 0)
	{
		throw std::overflow_error("Divide by zero exception");
	}

	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] /= value;
	}
}


KwanMat& KwanMat::MatMul(KwanMat& ref)
{
	if (numCols != ref.numRows)
	{
		throw std::out_of_range("Matrix multiplication should match the size");
	}

	KwanMat& theMat = Zero(numRows, ref.numCols);

	for (unsigned int row = 0; row < numRows; row++) {
		for (unsigned int col = 0; col < ref.numCols; col++) {
			for (unsigned int inner = 0; inner < numCols; inner++) {
				theMat.matNums[row*ref.numCols+col] += matNums[row*numCols+inner] * ref.matNums[inner*ref.numCols+col];
			}
		}
	}
	return theMat;
}

KwanMat& KwanMat::T()
{
	//copy temporally
	float* tempNums = new float[numRows*numCols];
	for (unsigned int i = 0; i < numRows * numCols; i++)
		tempNums[i] = matNums[i];
	for (unsigned int i = 0; i < numRows; i++)
	{
		for (unsigned int j = 0; j < numCols; j++)
		{
			matNums[i*numCols+j] = tempNums[j * numCols + i];
		}
	}
	delete[] tempNums;

	int tempN = numRows;
	numRows = numCols;
	numCols = tempN;

	return *this;
}


KwanMat& KwanMat::operator+(const KwanMat& ref)
{
	if (numRows != ref.numRows || numCols != ref.numCols)
	{
		throw std::out_of_range("Size should be Same for + operation");
	}

	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] += ref.matNums[i];
	}
	return *this;
}

KwanMat& KwanMat::operator-(const KwanMat& ref)
{
	if (numRows != ref.numRows || numCols != ref.numCols)
	{
		throw std::out_of_range("Size should be Same for + operation");
	}

	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] -= ref.matNums[i];
	}
	return *this;
}

KwanMat& KwanMat::operator*(const KwanMat& ref)
{
	if (numRows != ref.numRows || numCols != ref.numCols)
	{
		throw std::out_of_range("Size should be Same for + operation");
	}

	for (unsigned int i = 0; i < numRows * numCols; i++)
	{
		matNums[i] *= ref.matNums[i];
	}
	return *this;
}


int KwanMat::getRowSize()
{
	return numRows;
}

int KwanMat::getColumnSize()
{
	return numCols;
}

/// <summary>
/// Functions
/// </summary>
KwanMat& Identity(const unsigned int len)
{
	KwanMat* theMat = new KwanMat(len, len, 0);

	for (unsigned int i = 0; i < len * len; i++)
	{
		if (i / len == i % len)
			theMat->matNums[i] = 1;
	}
	return *theMat;
}

KwanMat& Zero(const unsigned int N, const unsigned int M)
{
	KwanMat* theMat = new KwanMat(N, M, 0);
	return *theMat;
}
