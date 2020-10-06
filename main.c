#include <stdio.h>
#include "geodesy.h"
#include "commons.h"
#define PI 3.1415926
int main()
{
	double a[2][3]={{1,2,3},{4,5,6}};
	double b[2][3]={{4,5,6},{1,2,3}};
	double c[2][3]={{4,5,6},{1,2,3}};
	struct Matrix *aMat=initMatrix(2,3,a);
	struct Matrix *bMat=initMatrix(2,3,b);
	struct Matrix *cMat=initMatrix(2,3,c);
	struct Matrix *add=matrixAddition(3,aMat,bMat,cMat);
	printMatrix(add);
	return 0;

}
