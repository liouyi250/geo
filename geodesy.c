#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "geodesy.h"

#define PI 3.1415926

struct Ellipsoid ellips[3]={{6378245,6356863.0188,0.00669342162297},//bj54
			   {6378140,6356755.2882,0.00669438499959},//xa80
			   {6378137,6356752.3141,0.00669438002290}};//cgcs2000

struct Ellipsoid *ell = ellips ;


//法线长
double normalLine(double b)
{
	return ell->a / sqrt(1 - ell->e2 * sin(b) * sin(b));
}

//大地坐标转空间直角坐标
struct Cartisian *geoToCart(double b,double l,float h)
{
	double n=normalLine(b);
	double x=(n + h) * cos(b) * cos(l);
	double y=(n + h) * cos(b) * sin(l);
	double z=(n * (1 - ell->e2) + h) * sin(b);
	struct Cartisian *cartisian = (struct Cartisian*)malloc(sizeof(struct Cartisian));
	cartisian->x=x;
	cartisian->y=y;
	cartisian->z=z;
	return cartisian;
}

//空间直角坐标转大地坐标
struct Geodetic *cartToGeo(float x,float y,float z)
{
	double b0,b1;
	struct Geodetic *geo=(struct Geodetic*)malloc(sizeof(struct Geodetic));
	double l=atan(y / x);
	b0=z/sqrt(x * x + y * y);
	b1=(z + ell->a * ell->e2 * b0 / sqrt(1 + b0 * b0 - ell->e2 * b0 * b0))/sqrt(x*x + y*y);
	while(fabs(atan(b0) - atan(b1)) * 180 / PI * 3600 > 0.0001)
	{
		b0 = b1;
		b1=(z + ell->a * ell->e2 * b0 / sqrt(1 + b0 * b0 - ell->e2 * b0 * b0))/sqrt(x*x + y*y);
	}
	b1 = atan(b1);
	double n=normalLine(b1);
	float h=sqrt(x*x + y*y) / cos(b1) - n;

	geo->b=b1;
	geo->l=l;
	geo->h=h;
	return geo;
}

//子午线弧长(b1纬度低，b2纬度高，参数为弧度)
double meridianArcLen(double b1,double b2)
{
	double A=1 + 3 * ell->e2 / 4 + 45 * ell->e2 * ell->e2 / 64 + 175 * ell->e2 * ell->e2 * ell->e2 / 256 + 11025 * ell->e2 * ell->e2 * ell->e2 * ell->e2 / 16384;
	double B=3 * ell->e2 /4 + 15 * ell->e2 * ell->e2 / 16 + 525 * ell->e2 * ell->e2 * ell->e2 /512 + 2205 * pow(ell->e2,4) / 2048;
	double C = 15 * ell->e2 * ell->e2 / 64 + 105 * ell->e2 * ell->e2 * ell->e2 / 256 + 2205 *pow(ell->e2,4) / 4096;
	double D = 35 * pow(ell->e2,3) / 512 + 315 * pow(ell->e2,4) / 2048;
	double x= ell->a * (1 - ell->e2) * (A * (b2 - b1) - B * (sin(2 * b2)- sin(2 * b1)) / 2 + C * (sin(4 * b2) - sin(4 * b1)) / 4 - D * (sin(6 * b2) - sin(6 * b1)) / 6);
	return x;
}

//高斯投影正算（从经纬度到平面坐标）
struct Cartisian *guassPositive(double b,double l)
{
	double m=l * cos(b);
	double t = tan(b);
	double sita=cos(b) * sqrt(ell->e2 / (1 - ell->e2));
	double n=normalLine(b);
	double x=meridianArcLen(0,b) + n * t * (m * m / 2 + (5 - t * t + 9 * sita * sita + 4 * pow(sita,4)) * pow(m,4) / 24 + pow(m,6) * (61 - 58 * t * t + pow(t,4)/720));
	double y=n * (m + (1 - t * t + sita * sita) * pow(m,3) / 6 + pow(m,5) * (5 - 18 * t * t + pow(t,4) + 14 * sita * sita - 58 * t * t * sita * sita) / 180);
	struct Cartisian *cart = (struct Cartisian*)malloc(sizeof(struct Cartisian));
	cart->x=x;
	cart->y=y;
	cart->z=0;

	return cart;
}

//求垂足纬度
double pedalLat(double x)
{
	double A=1 + 3 * ell->e2 / 4 + 45 * ell->e2 * ell->e2 / 64 + 175 * ell->e2 * ell->e2 * ell->e2 / 256 + 11025 * ell->e2 * ell->e2 * ell->e2 * ell->e2 / 16384;
	double b0=0;
	double b1=x / ell->a / (1 - ell->e2) / A;
	while(fabs(b1 - b0) * 180 / PI * 3600 > 0.0001)
	{
		b0 = b1;
		double arcLen=meridianArcLen(0,b1);
		double delta=x - arcLen;
		b1 = b1 + delta / ell->a / (1 - ell->e2) / A;
	}
	return b1;
}

//高丝投影反算
struct Geodetic *guassNegative(double x,double y)
{
	double b0=pedalLat(x);
	double n=normalLine(b0);
	double t=tan(b0);
	double v2=1 + ell->e2 * cos(b0) * cos(b0) / (1 - ell->e2);
	double sita=cos(b0) * sqrt(ell->e2 / (1 - ell->e2));
	double b=b0 -v2 * t * (y * y / n / n  - pow(y/n,4) * (5 + 3 * t *t + sita * sita - 9 *sita *sita * t *t) / 12 + pow(y/n,6) * (61+90*t*t+45*pow(t,4)) / 360) / 2;
	double l=(y / n - pow(y / n,3) * (1 + 2 * t * t + sita * sita) / 6 + pow(y / n,5)*(5 + 28 * t * t + 24 * pow(t,4) + 6 * sita * sita + 8 * t * t * sita * sita) / 120) / cos(b0);
	struct Geodetic *geo=(struct Geodetic*)malloc(sizeof(struct Geodetic));
	geo->b=b;
	geo->l=l;
	geo->h=0;
	return geo;

}

//大地问题正算
struct Bessel *besselPositive(double b1,double l1,double A1,double S)
{
	double u1=atan(tan(b1) * sqrt(1 - ell->e2));
	double m=asin(cos(u1) * sin(A1));
	if(m<0) m=m+2*PI;
	double M=atan(tan(u1) / cos(A1));
	if(M<0) M=M+PI;
	double k2=cos(m) * cos(m) * ell->e2 / (1 - ell->e2);
	double alpha=(1 - k2/4 + 7*k2*k2/64 - 15*pow(k2,3)/256)*sqrt(1+ell->e2/(1-ell->e2))/ell->a;
	double beta=k2/4 - k2*k2/8 +37*pow(k2,3)/512;
	double gammar=k2*k2/128 - pow(k2,3)/128;
	double sigma0=alpha * S;
	double sigma1=0;
	while(fabs(sigma0 - sigma1) * 3600*180/PI > 0.001)
	{
		sigma1=sigma0;
		sigma0=alpha*S + beta*sin(sigma1)*cos(2*M+sigma1) + gammar*sin(2*sigma1)*cos(4*M+2*sigma1);

	}
	double A2=atan(tan(m) / cos(M+sigma0));
	if(A2<0) A2=A2+PI;
	if(A1<PI) A2=A2+PI;
	double u2=atan(-cos(A2) * tan(M+sigma0));
	double b2=atan(tan(u2) * sqrt(1 + ell->e2/(1-ell->e2)));
	double lambd1=atan(sin(u1) * tan(A1));
	if(lambd1<0) lambd1=lambd1+PI;
	if(m>PI) lambd1=lambd1+PI;
	double lambd2=atan(sin(u2) * tan(A2));
	if(lambd2<0) lambd2=lambd2+PI;
	if((m<PI && M+sigma0>PI) || (m>=PI && M+sigma0<=PI))
		lambd2=lambd2+PI;
	double k2dot=ell->e2 * cos(m) * cos(m);
	double alphadot=ell->e2/2 + ell->e2*ell->e2/8 + pow(ell->e2,3)/16 - k2dot*ell->e2*(1+ell->e2)/16 + 3*ell->e2*k2dot*k2dot/128;
	double betadot=ell->e2*(1+ell->e2)*k2dot/16 - ell->e2*k2dot*k2dot/32;
	double gammardot=ell->e2*k2dot*k2dot/256;
	double l0=lambd2-lambd1 - sin(m)*(alphadot*sigma0 + betadot*sin(sigma0)*cos(2*M+sigma0) + gammardot*sin(2*sigma0)*cos(4*M+2*sigma0));
	double l2=l0+l1;
	struct Bessel *bessel=(struct Bessel*)malloc(sizeof(struct Bessel));
	bessel->b=b2;
	bessel->l=l2;
	bessel->A2=A2;

	return bessel;
}

//大地问题反算
//大地线和书上的例子还有米级的误差
struct Bessel *besselNegative(double b1,double l1,double b2,double l2)
{
	double u1=atan(tan(b1) * sqrt(1-ell->e2));
	double u2=atan(tan(b2) * sqrt(1-ell->e2));
	double alp=ell->e2*(0.5 + ell->e2 / 8 + ell->e2 * ell->e2 /16);

	double deltalam=0;
	double deltasig=0;
	double sigma0=acos(sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(l2-l1+deltalam));
	double m0=asin(cos(u1) * cos(u2) * sin(l2-l1+deltalam) / sin(sigma0));

	deltalam=alp * sigma0 * sin(m0);
	deltasig=sin(m0) * deltalam;
	sigma0=deltasig + sigma0;
	m0=asin(cos(u1) * cos(u2) * sin(l2-l1+deltalam) / sin(sigma0));
	double A10=atan(sin(l2-l1+deltalam) /(cos(u1) * tan(u2) - sin(u1)*cos(l2-l1+deltalam)));
	if(A10<0) A10=A10+PI;
	if(m0<0) A10=A10+PI;
	double M=atan(sin(u1)*tan(A10)/sin(m0));
	if(M<0) M=M+PI;
	double k2=ell->e2 * cos(m0)*cos(m0);
	double alpha=alp + pow(ell->e2,3)/16 - ell->e2*k2*(1+ell->e2)/16 + 3*ell->e2*k2*k2/128;
	double beta=ell->e2*k2*((1+ell->e2)/16 - k2/32);
	double lambd=l2-l1 + sin(m0)*(alpha*sigma0 + beta*sin(sigma0)*cos(2*M+sigma0));
	sigma0=acos(sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lambd));
	printf("sigma=%.10lf\n",sigma0);

	double k21=cos(m0) * cos(m0) * ell->e2 / (1 - ell->e2);
	double alpha1=(1 - k21/4 + 7*k21*k21/64 - 15*pow(k21,3)/256)*sqrt(1+ell->e2/(1-ell->e2))/ell->a;
	double beta1=k21/4 - k21*k21/8 +37*pow(k21,3)/512;
	double gammar1=k21*k21/128 - pow(k21,3)/128;

	double S=(sigma0 - beta1*sin(sigma0)*cos(2*M+sigma0) - gammar1*sin(2*sigma0)*cos(4*M+2*sigma0))/alpha1;
	double A1=atan(sin(lambd) / (cos(u1)*tan(u2) -sin(u1)*cos(lambd)));
	if(A1<0) A1=A1+PI;
	if(m0<=0) A1=A1+PI;
	double A2=atan(sin(lambd)/(sin(u2)*cos(lambd)- tan(u1)*cos(u2)));
	if(A2<0) A2=A2+PI;
	if(m0>=0) A2=A2+PI;

	struct Bessel *bessel=(struct Bessel*)malloc(sizeof(struct Bessel));
	bessel->s=S;
	bessel->A1=A1;
	bessel->A2=A2;

	return bessel;
}

//平面四参数计算
/*
oldAxis 转换前的坐标
newAxis 转换后的坐标
n 公共点个数
返回值共七个数: x方向平移，y方向平移，xyz旋转角度（秒），缩放因子
*/

/*
测试说明：测试一组数据，y方向误差在0.5左右，z方向误差在0.1，x方向在小数点后两位
*/
struct Matrix *planeTransfer(struct Cartisian *oldAxis,struct Cartisian *newAxis,int n)
{
	int i;
	if(n<2) return NULL;
	struct Matrix *l=mallocMatrix(2*n,1);
	struct Matrix *B=mallocMatrix(2*n,4);
	for(i=0;i<n;i++)
	{
		l->data[2*i][0]=newAxis[i].x;
		l->data[2*i+1][0] = newAxis[i].y;

		B->data[2*i][0]=1;
		B->data[2*i][1]=0;
		B->data[2*i][2]=oldAxis[i].x;
		B->data[2*i][3]=-oldAxis[i].y;
		B->data[2*i+1][0]=0;
		B->data[2*i+1][1]=1;
		B->data[2*i+1][2]=oldAxis[i].y;
		B->data[2*i+1][3]=oldAxis[i].x;
	}
	struct Matrix *Btrans=transposition(B);
	struct Matrix *normal=matrixMultMatrix(2,Btrans,B);
	struct Matrix *N=inverseMatrix(normal);

	struct Matrix *x=matrixMultMatrix(3,N,Btrans,l);
	float a=x->data[2][0];
	float b=x->data[3][0];
	float m=sqrt(a*a + b*b);
	float alpha=atan(b/a);
	x->data[2][0]=alpha;
	x->data[3][0]=m;

	freeMatrix(l);
	freeMatrix(B);
	freeMatrix(Btrans);
	freeMatrix(normal);
	freeMatrix(N);

	return x;
}

//平面七参数计算
/*
返回值：x方向平移，y方向平移，z方向平移，x旋转因子，y旋转因子，z旋转因子，缩放因子
*/
struct Matrix *spaceTransfer(struct Cartisian *oldAxis,struct Cartisian *newAxis,int n)
{
	if(n<3) return NULL;

	struct Matrix *l=mallocMatrix(3*n,1);
	struct Matrix *B=mallocMatrix(3*n,7);
	for(int i=0;i<n;i++)
	{
		l->data[3*i][0]=newAxis[i].x;
		l->data[3*i+1][0] = newAxis[i].y;
		l->data[3*i+2][0]=newAxis[i].z;

		B->data[3*i][0]=1;
		B->data[3*i][1]=0;
		B->data[3*i][2]=0;
		B->data[3*i][3]=oldAxis[i].x;
		B->data[3*i][4]=0;
		B->data[3*i][5]=-oldAxis[i].z;
		B->data[3*i][6]=oldAxis[i].y;

		B->data[3*i+1][0]=0;
		B->data[3*i+1][1]=1;
		B->data[3*i+1][2]=0;
		B->data[3*i+1][3]=oldAxis[i].y;
		B->data[3*i+1][4]=oldAxis[i].z;
		B->data[3*i+1][5]=0;
		B->data[3*i+1][6]=-oldAxis[i].x;

		B->data[3*i+2][0]=0;
		B->data[3*i+2][1]=0;
		B->data[3*i+2][2]=1;
		B->data[3*i+2][3]=oldAxis[i].z;
		B->data[3*i+2][4]=-oldAxis[i].y;
		B->data[3*i+2][5]=oldAxis[i].x;
		B->data[3*i+2][6]=0;
	}
	struct Matrix *Btrans=transposition(B);
	struct Matrix *normal=matrixMultMatrix(2,Btrans,B);

	struct Matrix *N=inverseMatrix(normal);

	struct Matrix *x=matrixMultMatrix(3,N,Btrans,l);
	float a=x->data[3][0];
	float b=x->data[4][0];
	float c=x->data[5][0];
	float d=x->data[6][0];

	x->data[6][0]=a ;
	x->data[3][0]=b/a *180 / PI *3600;
	x->data[4][0]=c/a *180 / PI *3600;
	x->data[5][0]=d/a *180 / PI *3600;

	freeMatrix(B);
	freeMatrix(l);
	freeMatrix(Btrans);
	freeMatrix(normal);
	freeMatrix(N);

	return x;
}
