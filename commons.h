#ifndef COMMONS_H_
#define COMMONS_H_

struct Point{
float x;
float y;
};

struct Line
{
	float A;
	float B;
	float C;
};//Ax+By+C=0

struct Matrix
{
	int rows;
	int columns;
    double **data;
};
double dmsToDegree(double dms);
float planeDist(struct Point *ptA,struct Point *ptB);
struct Point* intersection(struct Line *aLine,struct Line *bLine);
float dimangular(struct Line *aLine,struct Line *bLine);
struct Line* makeLine(struct Point *ptA,struct Point *ptB);
float lineDist(struct Line *line,struct Point *pt);

//矩阵运算
//分配矩阵空间
struct Matrix *mallocMatrix(int rows,int columns);

//initialize
struct Matrix *initMatrix(int rows,int columns,double a[rows][columns]);

//print matrix
void printMatrix(struct Matrix *mat);

//矩阵加矩阵
struct Matrix *matrixAddition(int n,...);

//矩阵减矩阵
struct Matrix *matrixSubtration(int n,...);

//矩阵乘以常数
struct Matrix* matrixMultConst(struct Matrix *mat,double times);

//矩阵乘以矩阵
struct Matrix *matrixMultMatrix(int n,...);
struct Matrix *matrixMultMatrix0(struct Matrix *mat1,struct Matrix *mat2);

//逆矩阵
struct Matrix *inverseMatrix(struct Matrix *mat);

//转置矩阵
struct Matrix *transposition(struct Matrix *mat);
//释放矩阵内存
void freeMatrix(struct Matrix *mat);
#endif
