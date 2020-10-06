#ifndef GEODESY_H_
#define GEODESY_H_

#include "commons.h"
/*
参考书籍：大地测量学基础（测绘出版社）2010年3月第一版
*/

//定义椭球
struct Ellipsoid
{
	float a;
	float b;
	double e2;
};

//笛卡尔坐标系
struct Cartisian
{
	double x;
	double y;
	double z;
};

//大地坐标系
struct Geodetic
{
	double b;
	double l;
	float  h;
};

//白塞尔大地问题解算返回类型
struct Bessel
{
	double b;
	double l;
	double s;
	double A1;
	double A2;
};


//全局变量 特定椭球
extern struct Ellipsoid *ell;

//global 
extern struct Ellipsoid ellips[3];

//法线长
double normalLine(double b);

//大地坐标转空间直角坐标
struct Cartisian *geoToCart(double b,double l,float h);

//空间直角坐标转大地坐标
struct Geodetic *cartToGeo(float x,float y,float z);

//子午线弧长
double meridianArcLen(double b1,double b2);

//高斯投影正算
struct Cartisian *guassPositive(double b,double l);

//高丝投影反算
struct Geodetic *guassNegative(double x,double y);

//大地问题正算
struct Bessel *besselPositive(double b,double l,double A1,double S);

//大地问题反算
struct Bessel *besselNegative(double b1,double l1,double b2,double l2);

//平面四参数计算
/*
返回值：x方向平移，y方向平移，旋转因子，缩放因子
*/
struct Matrix *planeTransfer(struct Cartisian *oldAxis,struct Cartisian *newAxis,int n);

//平面七参数计算
struct Matrix *spaceTransfer(struct Cartisian *oldAxis,struct Cartisian *newAxis,int n);
#endif

