
#matrix with 1d array based DLL.



##ex) transpose implementation with
<code>
>  for (unsigned int i = 0; i < numRows; i++)
>  {
>  	for (unsigned int j = 0; j < numCols; j++)
>  	{
>  		matNums[i*numCols+j] = tempNums[j * numCols + i];
>  	}
>  }
</code>


##ex) matrix multiplication with

<code>
>  for (unsigned int row = 0; row < numRows; row++) {
>  		for (unsigned int col = 0; col < ref.numCols; col++) {
>  			for (unsigned int inner = 0; inner < numCols; inner++) {
>  				theMat.matNums[row*ref.numCols+col] += matNums[row*numCols+inner] * ref.matNums[inner*ref.numCols+col];
>  			}
>  		}
>  }
</code>
