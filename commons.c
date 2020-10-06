#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "commons.h"
#define PI 3.1415926
#define true 1
#define false 0
typedef int bool;

/*
从角度换算为弧度
参数格式:[-]ddmmss.ssss
*/
double dmsToDegree(double dms)
{
	bool flag=false;
	if(dms<0)
	{
		flag=true;
		dms=-dms;
	}
	int d=(int)(dms / 10000);
	int m=(int)((dms - d*10000) / 100);
	double s=dms - 10000 * d - 100 * m;
	printf("%d %d %lf\n",d,m,s);
	double degree = (d + m/60.0 + s/3600)* PI / 180;

	return flag ? -degree : degree;
}

/*
计算两点间的平面距离
*/
float planeDist(struct Point *ptA,struct Point *ptB)
{
	float dist=sqrt((ptA->x-ptB->x)*(ptA->x-ptB->x)+(ptA->y-ptB->y)*(ptA->y-ptB->y));
	return dist;
}

/*
计算两条直线的交点
*/
struct Point *intersection(struct Line *aLine,struct Line *bLine)
{
	if(aLine->A==0 && bLine->A==0 && aLine->B!=bLine->B)
		return NULL;
	if(aLine->B==0 && bLine->B==0 && aLine->A!=bLine->A)
		return NULL;
	if(aLine->A*bLine->B == aLine->B*bLine->A)
		return NULL;
	struct Point *pt=(struct Point*)malloc(sizeof(struct Point));
	if(aLine->A==0 && bLine->B==0)
	{
		pt->x=-bLine->C / bLine->A;
		pt->y=-aLine->C / aLine->B;
		return pt;
	}
	if(aLine->B==0 && bLine->A==0)
	{
		pt->x=-aLine->C / aLine->A;
		pt->y=-bLine->C / bLine->B;
		return pt;
	}
	float x=(aLine->C * bLine->B - bLine->C * aLine->B) / (aLine->B * bLine->A - bLine->B * aLine->A);
	float y=(bLine->C * aLine->A - aLine->C * bLine->A) / (aLine->B * bLine->A - bLine->B * aLine->A);
	pt->x=x;
	pt->y=y;
	return pt;
}
/*
两点确定一条直线
*/
struct Line *makeLine(struct Point *ptA,struct Point *ptB)
{
	struct Line *line=(struct Line*)malloc(sizeof(struct Line));
	if(ptA->x == ptB->x && ptA->y == ptB->y)
		return NULL;
	else if(ptA->x == ptB->x)
	{
		line->A=1;
		line->B=0;
		line->C=-ptA->x;
	}
	else if(ptA->y == ptB->y)
	{
		line->A=0;
		line->B=1;
		line->C=-ptA->y;
	}
	else
	{
		line->A=(ptA->y - ptB->y) / (ptB->x - ptA->x);
		line->B=1;
		line->C=(ptB->x * ptA->y - ptA->x * ptB->y)/(ptA->x - ptB->x);
	}
	return line;
}
/*
计算两条直线的夹角，返回弧度
*/
float dimangular(struct Line *aLine,struct Line* bLine)
{
	if(aLine->B==0)
	{
		if(bLine->B==0) return 0;
		if(bLine->A==0) return PI/2;
		return PI/2 - atan(-bLine->A / bLine->B);
	}
	else
	{
		if(bLine->B==0) return PI/2-atan(-aLine->A / aLine->B);
		float k1=-aLine->A/aLine->B;
		float k2=-bLine->A/bLine->B;
		return atan(fabs(k1-k2)/(1+k1*k2));
	}
}

/*
计算点到直线的距离
*/
float lineDist(struct Line *line,struct Point *pt)
{
	return fabs(line->A*pt->x + line->B*pt->y + line->C)/sqrt(line->A*line->A + line->B*line->B);
}



//Matrix operation
//分配矩阵空间
struct Matrix *mallocMatrix(int rows,int columns)
{
	struct Matrix *mat=(struct Matrix*)malloc(sizeof(struct Matrix));
	int i;
	mat->rows = rows;
	mat->columns = columns;
	mat->data=(double**)malloc(mat->rows * sizeof(double*));
	for(i=0;i<rows;i++)
		mat->data[i]=(double*)malloc(sizeof(double) * mat->columns);
	return mat;
}

//initialize
struct Matrix *initMatrix(int rows,int columns,double a[rows][columns])
{
	int i,j;
	struct Matrix *mat=mallocMatrix(rows,columns);
	for(i=0;i<rows;i++)
		for(j=0;j<columns;j++)
			mat->data[i][j]=a[i][j];
	return mat;
}

//打印矩阵元素
void printMatrix(struct Matrix *mat)
{
	int i,j;
	if(mat == NULL)
		return;
	for(i=0;i<mat->rows;i++){
		for(j=0;j<mat->columns;j++){
			printf("%.7lf ",mat->data[i][j]);
		}
		printf("\n");
	}
}

//矩阵元素深复制
struct Matrix *deepCopy(struct Matrix *srcMat)
{
	struct Matrix *mat=mallocMatrix(srcMat->rows,srcMat->columns);
	for(int i=0;i<srcMat->rows;i++)
		for(int j=0;j<srcMat->columns;j++)
			mat->data[i][j]=srcMat->data[i][j];

	return mat;
}

//矩阵加矩阵
struct Matrix *matrixAddition0(struct Matrix *mat1,struct Matrix *mat2)
{
	int i,j;
	if(mat1->rows != mat2->rows || mat1->columns != mat2->columns)
		return NULL;
	struct Matrix *mat=mallocMatrix(mat1->rows,mat1->columns);
	for(i=0;i<mat->rows;i++)
		for(j=0;j<mat->columns;j++)
		{
			mat->data[i][j]=mat1->data[i][j]+mat2->data[i][j];
		}
	return mat;
}

//矩阵加矩阵
struct Matrix *matrixAddition(int n,...)
{
	va_list ap;
	struct Matrix *mat=NULL;
	struct Matrix *mat1,*mat2;
	struct Matrix *matArr[n-1];
	va_start(ap,n);
	mat1=va_arg(ap,struct Matrix*);
	mat2=va_arg(ap,struct Matrix*);
	matArr[0]=matrixAddition0(mat1,mat2);

	for(int i=2;i<n;i++)
	{
		mat1=matArr[i-2];
		mat2=va_arg(ap,struct Matrix*);
		matArr[i-1]=matrixAddition0(mat1,mat2);
	}

	va_end(ap);

	mat=deepCopy(matArr[n-2]);
	for(int i=0;i<n-1;i++)
		freeMatrix(matArr[i]);

	return mat;
}

//矩阵减矩阵
struct Matrix *matrixSubtration0(struct Matrix *mat1,struct Matrix *mat2)
{

	int i,j;
	if(mat1->rows != mat2->rows || mat1->columns != mat2->columns)
		return NULL;
	struct Matrix *mat=mallocMatrix(mat1->rows,mat1->columns);
	for(i=0;i<mat->rows;i++)
		for(j=0;j<mat->columns;j++)
		{
			mat->data[i][j]=mat1->data[i][j]-mat2->data[i][j];
		}
	return mat;
}

//矩阵减矩阵
struct Matrix *matrixSubtration(int n,...)
{
	va_list ap;
	struct Matrix *mat=NULL;
	struct Matrix *mat1,*mat2;
	struct Matrix *matArr[n-1];
	va_start(ap,n);
	mat1=va_arg(ap,struct Matrix*);
	mat2=va_arg(ap,struct Matrix*);
	matArr[0]=matrixSubtration0(mat1,mat2);

	for(int i=2;i<n;i++)
	{
		mat1=matArr[i-2];
		mat2=va_arg(ap,struct Matrix*);
		matArr[i-1]=matrixSubtration0(mat1,mat2);
	}

	va_end(ap);

	mat=deepCopy(matArr[n-2]);
	for(int i=0;i<n-1;i++)
		freeMatrix(matArr[i]);

	return mat;
}

//矩阵乘以常数
struct Matrix* matrixMultConst(struct Matrix *mat,double times)
{

	int i,j;
	struct Matrix *retMat=mallocMatrix(mat->rows,mat->columns);
	for(i=0;i<mat->rows;i++)
		for(j=0;j<mat->columns;j++)
		{
			retMat->data[i][j]=mat->data[i][j] * times;
		}
	return retMat;
}


//矩阵乘以矩阵
struct Matrix *matrixMultMatrix(int n,...)
{
	va_list ap;
	struct Matrix *mat=NULL;
	struct Matrix *mat1,*mat2;
	struct Matrix *matArr[n-1];
	va_start(ap,n);
	mat1=va_arg(ap,struct Matrix*);
	mat2=va_arg(ap,struct Matrix*);
	matArr[0]=matrixMultMatrix0(mat1,mat2);

	for(int i=2;i<n;i++)
	{
		mat1=matArr[i-2];
		mat2=va_arg(ap,struct Matrix*);
		matArr[i-1]=matrixMultMatrix0(mat1,mat2);
	}

	va_end(ap);

	mat=deepCopy(matArr[n-2]);
	for(int i=0;i<n-1;i++)
		freeMatrix(matArr[i]);

	return mat;
}

struct Matrix *matrixMultMatrix0(struct Matrix *mat1,struct Matrix *mat2)
{
	struct Matrix *mat=NULL;
	int i,j,k;
	double m=0;
	if(mat1->columns != mat2->rows)
		return NULL;
	mat=mallocMatrix(mat1->rows,mat2->columns);
	for(i=0;i<mat1->rows;i++){
		for(j=0;j<mat2->columns;j++){
			m=0;
			for(k=0;k<mat1->columns;k++){
			m=m+mat1->data[i][k]*mat2->data[k][j];
			}
			mat->data[i][j]=m;
		}
	}
	return mat;
}
//矩阵的两个行相互交换
void matrixRowChange(int row1,int row2,struct Matrix *mat1,struct Matrix *mat2)
{
	int i;
	int k=mat1->columns + mat2->columns;
	double *arr=(double*)malloc(sizeof(double)*k);
	//将row1中的数据拷贝
	for(i=0;i<mat1->columns;i++)
		arr[i]=mat1->data[row1][i];
	for(i=0;i<mat2->columns;i++)
		arr[mat1->columns + i]=mat2->data[row1][i];
	//将row2中的数据赋值给row1
	for(i=0;i<mat1->columns;i++)
		mat1->data[row1][i]=mat1->data[row2][i];
	for(i=0;i<mat2->columns;i++)
		mat2->data[row1][i]=mat2->data[row2][i];
	//将拷贝数据赋值给row2
	for(i=0;i<mat1->columns;i++)
		mat1->data[row2][i]=arr[i];
	for(i=0;i<mat2->columns;i++)
		mat2->data[row2][i]=arr[mat1->columns + i];
	//回收内存
	free(arr);
	arr=NULL;
}

//增广矩阵放大倍数
void scaleMatrix(struct Matrix *mat1,struct Matrix *mat2,int row,double scale)
{
	int i;
	for(i=0;i<mat1->columns;i++)
		mat1->data[row][i] = mat1->data[row][i] * scale;
	for(i=0;i<mat2->columns;i++)
		mat2->data[row][i] = mat2->data[row][i] * scale;
}

//矩阵两行相加,row2的值最终变化
void matrixRowAddition(struct Matrix *mat1,struct Matrix *mat2,int row1,int row2,double scale)
{
	int i;
	for(i=0;i<mat1->columns;i++)
		mat1->data[row2][i] = mat1->data[row2][i] + mat1->data[row1][i] * scale;
	for(i=0;i<mat2->columns;i++)
		mat2->data[row2][i] = mat2->data[row2][i] + mat2->data[row1][i] * scale;
}

//逆矩阵,使用增广矩阵求逆
struct Matrix *inverseMatrix(struct Matrix *mat)
{
	bool flag=false;
	int i,j,k;
	double scale;
	if(mat->rows != mat->columns)
		return NULL;
	struct Matrix *invMat=mallocMatrix(mat->rows,mat->rows);
	for(i=0;i<mat->rows;i++)
		for(int j=0;j<mat->rows;j++){
			if(i != j)
				invMat->data[i][j]=0;
			else
				invMat->data[i][j]=1;
		}

	for(i=0;i<mat->rows;i++){
	//主对角线元素为零，向下查找不为零的元素
		if(mat->data[i][i] == 0){
			flag=true;
			for(k=i+1;k<mat->rows;k++)
			{
				if(mat->data[k][i] != 0)
				{
					flag=false;
					matrixRowChange(i,k,mat,invMat);
					break;
				}
			}
		}

		//某一列在某一行下全为零
		if(flag) return NULL;

		for(j=0;j<mat->rows;j++)
		{
			if(i==j) continue;
			if(mat->data[j][i]==0) continue;
			scale = -mat->data[j][i] / mat->data[i][i];
			matrixRowAddition(mat,invMat,i,j,scale);
		}
	}
	for(i=0;i<mat->rows;i++)
	{
		scale= 1.0 / mat->data[i][i];
		scaleMatrix(mat,invMat,i,scale);
	}
	return invMat;
}

//转置矩阵
struct Matrix *transposition(struct Matrix *mat)
{
	int i,j;
	struct Matrix *trans=mallocMatrix(mat->columns,mat->rows);
	for(i=0;i<mat->rows;i++)
		for(j=0;j<mat->columns;j++)
			trans->data[j][i]=mat->data[i][j];

	return trans;
}

//释放矩阵内存
void freeMatrix(struct Matrix *mat)
{
	if(mat == NULL) return;
	for(int i=0;i<mat->rows;i++)
		free(mat->data[i]);
	free(mat->data);
	free(mat);
}
