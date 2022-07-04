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

float KwanMat::RecursiveDet(const float mat[25], const int n)
{
	float det = 0;
	float subMat[25];
	if (n == 2)
		return mat[0] * mat[3] - mat[1] * mat[2];
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					subMat[subi*(n-1)+subj] = mat[i*n+j];
					subj++;
				}
				subi++;
			}
			det += x % 2 ? -mat[x] * RecursiveDet(subMat, n - 1) : mat[x] * RecursiveDet(subMat, n - 1);
			//det = det + (pow(-1, x) * mat[x] * RecursiveDet(submatrix, n - 1));
		}
	}
	return det;
}

void KwanMat::Cofactor(const float M[25],float t[16],int p,int q, int n)
{
	int i = 0, j = 0;
	for (int r = 0; r < n; r++) {
		for (int c = 0; c < n; c++) //Copy only those elements which are not in given row r and column c
			if (r != p && c != q) {
				t[i*(n-1)+j++] = M[r*n+c]; //If row is filled increase r index and reset c index
				if (j == n - 1) {
					j = 0; i++;
				}
			}
	}
}

void KwanMat::ADJ(float M[25], float adj[25])
{
	int s = 1;
	float t[25] = {0};
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numRows; j++) {
			//To get cofactor of M[i][j]
			Cofactor(M, t, i, j, numRows);
			s = ((i + j) % 2 == 0) ? 1 : -1; //sign of adj[j][i] positive if sum of row and column indexes is even.
			adj[j*numRows+i] = (s) * (RecursiveDet(t, numRows-1)); //Interchange rows and columns to get the transpose of the cofactor matrix
		}
	}
}



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
	//for general purpose
	//float* tempNums = new float[numRows*numCols];

	float tempNums[1024];

	//copy temporally
	for (unsigned int i = 0; i < numRows * numCols; i++)
		tempNums[i] = matNums[i];
	for (unsigned int i = 0; i < numRows; i++)
	{
		for (unsigned int j = 0; j < numCols; j++)
		{
			matNums[i*numCols+j] = tempNums[j * numCols + i];
		}
	}
	//delete[] tempNums;

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


float KwanMat::DET()
{
	if (!IsSquare(*this))
		throw std::runtime_error("Determinant is 0 or matrix is not square");

	//Laplace expansion
	if (numRows == 1) return matNums[0];
	else return RecursiveDet(matNums, numRows);
}

KwanMat& KwanMat::Inverse()
{
	float det = KwanMat::DET();
	if (det == 0) {
		throw std::runtime_error("Determinant is 0");
	}
	float dInv = 1 / det;
	
	float adj[25]; 
	ADJ(matNums, adj);
	for (int i = 0; i < numRows* numRows; i++)
		matNums[i] = adj[i] * dInv;

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


KwanMat& Identity(const unsigned int len)
{
	KwanMat* theMat = new KwanMat(len, len, 0);

	for (unsigned int i = 0; i < len; i++)
	{
		theMat->matNums[i*len+i] = 1;
	}
	return *theMat;
}

KwanMat& Zero(const unsigned int N, const unsigned int M)
{
	KwanMat* theMat = new KwanMat(N, M, 0);
	return *theMat;
}

KwanMat& Rotation2D(const float theta)
{
	float c = static_cast<float>(cos(theta));
	float s = static_cast<float>(sin(theta));
	float rotArr[4] = { c,-s,s,c };
	KwanMat* rotMat = new KwanMat(2, 2, rotArr);
	return *rotMat;
}

KwanMat& Rotation3D(const float theta,const char axis)
{
	float c = static_cast<float>(cos(theta));
	float s = static_cast<float>(sin(theta));
	float rotArr[9];
	if (axis == 'x' || axis == '0')
	{
		float rotArr[9] = {1,0,0,0, c,-s,0,s,c };
	}
	else if (axis == 'y' || axis == '1')
	{
		float rotArr[9] = { c,0,s,0,1,0,-s,0,c};
	}
	else if (axis == 'z' || axis == '2')
	{
		float rotArr[9] = { c,-s,0,s,c,0,0,0,1 };
	}
	else
	{
		throw std::out_of_range("axis should be x,y or z");
	}
	KwanMat* rotMat = new KwanMat(2, 2, rotArr);
	return *rotMat;
}

bool IsSquare(const KwanMat& mat)
{
	if (mat.numCols != mat.numRows) return false;
	else return true;
}

bool IsSymmetric(const KwanMat& mat)
{
	if(!IsSquare(mat))
		return false;

	int len = mat.numRows;
	for (unsigned int i = 0; i <= len/2; i++)
	{
		for (unsigned int j = 0; j <= len/2; j++)
		{
			if (mat.matNums[i * len + j] != mat.matNums[j * len + i]) return false;
		}
	}
	return true;
}