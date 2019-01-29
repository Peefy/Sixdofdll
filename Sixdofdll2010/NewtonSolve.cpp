
#include "stdafx.h"
#include "NewtonSolve.h"

#include <iostream>
#include <cmath>

using namespace std;

#define pi 3.141592653589793

#define N 6
#define EPSILON  0.0001
#define ITER_MAX 1

const int N2 = 2 * N;

double UpZ = 100;
double DownZ = 700;
double UpR = 680;
double DownR = 840;
double Dis = 190;

void ff(double xx[N], double yy[N],double dLen1, double dLen2, double dLen3, double dLen4, double dLen5, double dLen6);
void ffjacobian(double xx[N], double yy[N][N]);
void inv_jacobian(double yy[N][N], double inv[N][N]);
void newdundiedai(double x0[N], double inv[N][N], double y0[N], double x1[N]);

void ff(double xx[N], double yy[N], double dLen1, double dLen2, double dLen3, double dLen4, double dLen5, double dLen6)
{
	double acosUpR = acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR));
	double acosDownR = acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR));
	double acosDisUpR = acosUpR / 2.0 - pi / 6.0;
	double acosDisDownR = acosDownR / 2.0 - pi / 2.0;
	double cos_acosDisUpR = cos(acosDisUpR);
	double sin_acosDisUpR = sin(acosDisUpR);
	double cos_acosDisDownR = cos(acosDisDownR);
	double sin_acosDisDownR = sin(acosDisDownR);

	double x = xx[0];
	double y = xx[1];
	double z = xx[2];
	double a = xx[3];
	double b = xx[4];
	double c = xx[5];

	yy[0] = sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) /
		(2.0 * DownR * DownR)) / 2.0 - pi / 2.0) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
		2.0 - pi / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) /
		(2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*sin(c)), 2) + pow((UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x +
		DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 2.0) + UpR*
		sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) -
		UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*cos(c)), 2) + pow((DownZ + z -
		UpZ*cos(a)*cos(b) - UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*sin(b) +
		UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*sin(a)), 2))) - dLen1 -
		sqrt((pow((UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2 * UpR * UpR)) / 2.0 - pi / 6.0) -
		DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2 * DownR * DownR)) / 2.0 - pi / 2.0)), 2) +
		pow((DownZ - UpZ), 2) + pow((UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0) +
		DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 2.0)), 2)));
	yy[1] = sqrt((pow((UpZ*cos(a)*cos(b) - z - DownZ + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - 
		(5.0 * pi) / 6.0)*sin(b) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		cos(b)*sin(a)), 2) + 
		pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / 
		(2.0 * DownR * DownR)) / 2.0 - pi / 3.0) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / 
		(2.0 * DownR * DownR)) / 2.0 - pi / 3.0) - UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		cos(b)*sin(c)), 2))) - dLen2 - sqrt((pow((DownZ - UpZ), 2) + pow((UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - 
		(5.0 * pi) / 6.0) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 3.0)), 2) + 
		pow((UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0) + DownR*cos(acos((Dis * Dis - 2.0 * 
		DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 3.0)), 2)));
	yy[2] = sqrt((pow((UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) /
		(2.0 * DownR * DownR)) / 2.0 - (2.0 * pi) / 3.0) - UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
		2.0 - pi / 2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
		2.0 - pi / 2.0)*cos(b)*cos(c)), 2) + pow((DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 -
		(2.0 * pi) / 3.0) - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - y + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) /
		(2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) /
		(2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*cos(b)*sin(c)), 2) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*cos(acos((Dis * Dis -
		2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*sin(b) - UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
		2.0 - pi / 2.0)*cos(b)*sin(a)), 2))) - dLen3 - sqrt((pow((DownZ - UpZ), 2) + pow((UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) /
		(2.0 * UpR * UpR)) / 2.0 - pi / 2.0) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 -
		(2.0 * pi) / 3.0)), 2) + pow((UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0) + DownR*
		cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - (2.0 * pi) / 3.0)), 2)));
	yy[3] = sqrt((pow((UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) /
		2.0 - (2.0 * pi) / 3.0) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*(cos(a)*sin(c) - cos(c)*
		sin(a)*sin(b)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*cos(b)*cos(c)), 2) 
		+ pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*sin(b) + 
		UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*cos(b)*sin(a)), 2) + pow((y + UpZ*(cos(c)*
		sin(a) - cos(a)*sin(b)*sin(c)) + DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - (2.0 * pi) / 3.0) + 
		UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - 
		UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0)*cos(b)*sin(c)), 2))) - dLen4 - sqrt((pow((DownZ - UpZ), 2) +
		pow((UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 2.0) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) /
		(2.0 * DownR * DownR)) / 2.0 - (2.0 * pi) / 3.0)), 2) + pow((UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
		2.0 - pi / 2.0) + DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - (2 * pi) / 3)), 2)));
	yy[4] = sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - 
		(5.0 * pi) / 6.0)*sin(b) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*sin(acos((Dis * Dis - 2 * DownR * DownR) / 
		(2.0 * DownR * DownR)) / 2.0 - pi / 3.0) - UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*cos(b)*
		cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) /
		(2.0 * DownR * DownR)) / 2.0 - pi / 3.0) + UpR*sin(acos((Dis * Dis - 2 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0)*
		cos(b)*sin(c)), 2))) - dLen5 - sqrt((pow((DownZ - UpZ), 2) + pow((UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 -
		(5.0 * pi) / 6.0) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 3.0)), 2) +
		pow((UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - (5.0 * pi) / 6.0) + DownR*
		cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 3.0)), 2)));
	yy[5] = sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) /
		2.0 - pi / 2) - UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) +
		UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*sin(c)), 2) + pow((UpZ*cos(a)*cos(b) - z -
		DownZ + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*sin(b) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) /
		(2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 2.0) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) /
		(2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
		2.0 - pi / 6.0)*cos(b)*cos(c)), 2))) - dLen6 - sqrt((pow((UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2 - pi / 6) -
		DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 2.0)), 2) + pow((DownZ - UpZ), 2) +
		pow((UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) /
		(2.0 * DownR * DownR)) / 2.0 - pi / 2.0)), 2)));
}

void ffjacobian(double xx[N], double yy[N][N])
{


	double acosUpR = acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR));
	double acosDownR = acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR));
	double acosDisUpR = acosUpR / 2.0 - pi / 6.0;
	double acosDisDownR = acosDownR / 2.0 - pi / 2.0;
	double cos_acosDisUpR = cos(acosDisUpR);
	double sin_acosDisUpR = sin(acosDisUpR);
	double cos_acosDisDownR = cos(acosDisDownR);
	double sin_acosDisDownR = sin(acosDisDownR);

	double x = xx[0];
	double y = xx[1];
	double z = xx[2];
	double a = xx[3];
	double b = xx[4];
	double c = xx[5];

	yy[0][0] = -(2 * UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2.0 * x + 2.0 * 
		DownR*cos_acosDisDownR + 
		2.0 * UpR*sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - 2.0 *
		UpR*cos_acosDisUpR*cos(b)*cos(c)) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*
		sin(c)) + DownR*sin_acosDisDownR + UpR*sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 
		UpR*cos_acosDisUpR*cos(b)*sin(c)), 2) +
		pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*cos_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*cos_acosDisUpR*cos(b)*cos(c)), 2) + 
		pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*cos_acosDisUpR*sin(b) + 
		UpR*sin_acosDisUpR*cos(b)*
		sin(a)), 2))));
	yy[0][1] =  (2.0 * y + 2.0 * UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 2.0 * 
		DownR*sin_acosDisDownR + 2.0 * UpR*sin_acosDisUpR*(cos(a)*cos(c) + 
		sin(a)*sin(b)*sin(c)) + 2.0 * 
		UpR*cos_acosDisUpR*cos(b)*sin(c)) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*sin(a) - 
		cos(a)*sin(b)*
		sin(c)) + DownR*sin_acosDisDownR + UpR*sin_acosDisUpR*(cos(a)*cos(c) + 
		sin(a)*sin(b)*sin(c)) + UpR*cos_acosDisUpR*cos(b)*sin(c)), 2) + 
		pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*cos_acosDisDownR + 
		UpR*sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*
		sin(a)*sin(b)) - UpR*cos_acosDisUpR*cos(b)*cos(c)), 2) + pow((DownZ + z - UpZ*
		cos(a)*cos(b) - UpR*cos_acosDisUpR*sin(b) + UpR*sin_acosDisUpR*cos(b)*sin(a)), 2))));
	yy[0][2] = (2.0 * DownZ + 2.0 * z - 2.0 * UpZ*cos(a)*cos(b) - 2.0 * UpR*cos_acosDisUpR*
		sin(b) + 2.0 * UpR*sin_acosDisUpR*cos(b)*sin(a)) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*
		sin(a) - cos(a)*sin(b)*sin(c)) + DownR*sin_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c)), 2) + pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*
		cos(c)*sin(b)) + DownR*cos_acosDisDownR + UpR*sin_acosDisUpR*(cos(a)*sin(c) - 
		cos(c)*sin(a)*sin(b)) - UpR*cos_acosDisUpR*cos(b)*cos(c)), 2) + 
		pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*cos_acosDisUpR*sin(b) + 
		UpR*sin_acosDisUpR*cos(b)*sin(a)), 2))));
	yy[0][3] = (2.0 * (UpZ*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*sin_acosDisUpR*
		(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 
		DownR*sin_acosDisDownR + UpR*sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*
		sin(b)*sin(c)) + UpR*cos_acosDisUpR*cos(b)*sin(c)) + 
		2.0 * (UpZ*cos(b)*sin(a) + 
		UpR*sin_acosDisUpR*cos(a)*cos(b))*(DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos_acosDisUpR*sin(b) + UpR*sin_acosDisUpR*cos(b)*sin(a)) + 2.0 * 
		(UpZ*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*sin_acosDisUpR*(sin(a)*sin(c) + 
		cos(a)*cos(c)*sin(b)))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*cos_acosDisDownR + 
		UpR*sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos_acosDisUpR*cos(b)*cos(c))) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 
		DownR*sin_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c)), 2) + pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*
		cos(c)*sin(b)) + DownR*cos_acosDisDownR + UpR*sin_acosDisUpR*(cos(a)*sin(c) - 
		cos(c)*sin(a)*sin(b)) - UpR*cos_acosDisUpR*cos(b)*cos(c)), 2) + 
		pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*cos_acosDisUpR*sin(b) + 
		UpR*sin_acosDisUpR*cos(b)*sin(a)), 2))));
	yy[0][4] = -(2 * (UpR*cos_acosDisUpR*cos(b) - UpZ*cos(a)*sin(b) + UpR*
		sin_acosDisUpR*sin(a)*sin(b))*(DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos_acosDisUpR*sin(b) + UpR*sin_acosDisUpR*cos(b)*sin(a)) + 2.0 * (UpZ*cos(a)*cos(b)*sin(c) + UpR*
		cos_acosDisUpR*sin(b)*sin(c) - UpR*
		sin_acosDisUpR*cos(b)*sin(a)*sin(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		sin_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c)) - 2.0 * (UpZ*cos(a)*cos(b)*cos(c) + UpR*
		cos_acosDisUpR*cos(c)*sin(b) - UpR*
		sin_acosDisUpR*cos(b)*cos(c)*sin(a))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		cos_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos_acosDisUpR*cos(b)*cos(c))) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		sin_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c)), 2) + pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		cos_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos_acosDisUpR*cos(b)*cos(c)), 2) + pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos_acosDisUpR*sin(b) + UpR*
		sin_acosDisUpR*cos(b)*sin(a)), 2))));
	yy[0][5] = -(2 * (UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + UpR*
		sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos_acosDisUpR*cos(b)*cos(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		sin_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c)) - 2.0 * (UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		cos_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos_acosDisUpR*cos(b)*cos(c))) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		sin_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos_acosDisUpR*cos(b)*sin(c)), 2) + pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		cos_acosDisDownR + UpR*
		sin_acosDisUpR*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos_acosDisUpR*cos(b)*cos(c)), 2) + pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos_acosDisUpR*sin(b) + UpR*
		sin_acosDisUpR*cos(b)*sin(a)), 2))));

	yy[1][0] = (2.0 * x - 2.0 * UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + 2.0 * DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + 2.0 * UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2.0 * UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)) / (2.0 * sqrt((pow((-DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)), 2))));
	yy[1][1] = (2.0 * y + 2 * UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - 2 * DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - 2.0 * UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 2 * UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)) / (2 * sqrt((pow((-DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR/ 2.0- (5.0 * pi) / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)), 2))));
	yy[1][2] = -(2.0 * UpZ*cos(a)*cos(b) - 2.0 * z - 2.0 * DownZ + 2.0 * UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + 2.0 * UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)) / (2.0 * sqrt((pow((-DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)), 2))));
	yy[1][3] = -(2.0 * (UpZ*cos(b)*sin(a) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(a)*cos(b))*(UpZ*cos(a)*cos(b) - z - DownZ + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)) + 2.0 * (UpZ*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)) - 2.0 * (UpZ*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c))) / (2.0 * sqrt((pow((-DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)), 2))));
	yy[1][4] = -(2.0 * (UpZ*cos(a)*sin(b) - UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(a)*sin(b))*(UpZ*cos(a)*cos(b) - z - DownZ + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)) + 2.0 * (UpZ*cos(a)*cos(b)*cos(c) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(c)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)*sin(a))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)) + 2.0 * (UpZ*cos(a)*cos(b)*sin(c) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b)*sin(c) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)*sin(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c))) / (2.0 * sqrt((pow((UpZ*cos(a)*cos(b) - z - DownZ + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)), 2))));
	yy[1][5] = (2.0 * (UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)) - 2.0 * (UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c))) / (2.0 * sqrt((pow((-DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*sin(b) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(a)), 2) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR / 2.0 - pi / 3.0) + UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*cos(c)), 2) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		cos(acosDownR / 2.0 - pi / 3.0) - UpR*
		sin(acosUpR / 2.0 - (5.0 * pi) / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR / 2.0 - (5.0 * pi) / 6.0)*cos(b)*sin(c)), 2))));

	yy[2][0] = -(2.0*UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2.0*x + 2.0*DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - 2.0*UpR*sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((- y - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0))));
	yy[2][1] = -(2.0*DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) - 2.0*UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - 2.0*y + 2.0*UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((- y - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0))));
	yy[2][2] = (2.0*DownZ + 2.0*z - 2.0*UpZ*cos(a)*cos(b) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - 2.0*UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((- y - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0))));
	yy[2][3] = (2.0*(UpZ*cos(b)*sin(a) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(a)*cos(b))*(DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)) + 2.0*(UpZ*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)) - 2.0*(UpZ*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - y + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((- y - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0))));
	yy[2][4] = (2.0*(UpZ*cos(a)*cos(b)*sin(c) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b)*sin(c) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)*sin(c))*(DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - y + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)) + 2.0*(UpZ*cos(a)*sin(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*sin(a)*sin(b))*(DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)) + 2.0*(UpZ*cos(a)*cos(b)*cos(c) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(c)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)*sin(a))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)))/(2.0*sqrt((pow((UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - y + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0))));
	yy[2][5] = (2.0*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c))*(DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - y + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)) - 2.0*(UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((- y - UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0))));

	yy[3][0] = -(2.0*UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2.0*x + 2.0*DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + 2.0*UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0))));
	yy[3][1] =(2.0*y + 2.0*UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 2.0*DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + 2.0*UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - 2.0*UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0)+ pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0))));
	yy[3][2] = (2.0*DownZ + 2.0*z - 2.0*UpZ*cos(a)*cos(b) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + 2.0*UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0))));
	yy[3][3] = (2.0*(UpZ*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)) + 2.0*(UpZ*cos(b)*sin(a) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(a)*cos(b))*(DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)) + 2.0*(UpZ*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0))));
	yy[3][4] = (2.0*(UpZ*cos(a)*sin(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b) - UpR*
		sin(acosUpR/2.0 - pi/2.0)*sin(a)*sin(b))*(DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)) + 2.0*(UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b)*sin(c) - UpZ*cos(a)*cos(b)*sin(c) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)*sin(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)) - 2.0*(UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(c)*sin(b) - UpZ*cos(a)*cos(b)*cos(c) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)*sin(a))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0))));
	yy[3][5] = -(2.0*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)) - 2.0*(UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c))*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - x + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)))/(2.0*sqrt((pow((- x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*cos(c)),2.0) + pow((DownZ + z - UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/2.0)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*cos(b)*sin(a)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - (2.0*pi)/3.0) + UpR*
		sin(acosUpR/2.0 - pi/2.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		cos(acosUpR/2.0 - pi/2.0)*cos(b)*sin(c)),2.0))));

	yy[4][0] =  (2.0*x - 2.0*UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + 2.0*DownR*
		sin(acosDownR/2.0 - pi/3.0) - 2.0*UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2.0*UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c))/(2.0*sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)),2.0))));
	yy[4][1] = (2.0*y + 2.0*UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 2.0*DownR*
		cos(acosDownR/2.0 - pi/3.0) + 2.0*UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 2.0*UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c))/(2.0*sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)),2.0))));
	yy[4][2] = (2.0*DownZ + 2.0*z - 2.0*UpZ*cos(a)*cos(b) - 2.0*UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + 2.0*UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a))/(2.0*sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)),2.0))));
	yy[4][3] = (2.0*(UpZ*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)) - 2.0*(UpZ*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)) + 2.0*(UpZ*cos(b)*sin(a) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(a)*cos(b))*(DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)))/(2.0*sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)),2.0))));
	yy[4][4] = -(2.0*(UpZ*cos(a)*cos(b)*cos(c) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(c)*sin(b) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)*sin(a))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)) + 2.0*(UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b) - UpZ*cos(a)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*sin(a)*sin(b))*(DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)) + 2.0*(UpZ*cos(a)*cos(b)*sin(c) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b)*sin(c) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)*sin(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)))/(2.0*sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)),2.0))));
	yy[4][5] = -(2.0*(UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)) + 2.0*(UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)))/(2.0*sqrt((pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*
		sin(acosDownR/2.0 - pi/3.0) - UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*cos(c)),2.0) + pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + DownR*
		cos(acosDownR/2.0 - pi/3.0) + UpR*
		sin(acosUpR/2.0 - (5.0*pi)/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - (5.0*pi)/6)*cos(b)*sin(c)),2.0))));

	yy[5][0] = (2.0*x - 2.0*UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2.0*DownR*
		cos(acosDownR/2.0 - pi/2.0) + 2.0*UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c))/(2.0*sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)),2.0) + pow((- DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)),2.0))));
	yy[5][1] = (2.0*y + 2.0*UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - 2.0*DownR*
		sin(acosDownR/2.0 - pi/2.0) - 2.0*UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 2.0*UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c))/(2.0*sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)),2.0) + pow((- DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)),2.0))));
	yy[5][2] = -(2.0*UpZ*cos(a)*cos(b) - 2.0*z - 2.0*DownZ + 2.0*UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + 2.0*UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a))/(2.0*sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)),2.0) + pow((- DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)),2.0))));
	yy[5][3] = -(2.0*(UpZ*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		sin(acosUpR/2.0 - pi/6)*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)) + 2.0*(UpZ*cos(b)*sin(a) - UpR*
		sin(acosUpR/2.0 - pi/6)*cos(a)*cos(b))*(UpZ*cos(a)*cos(b) - z - DownZ + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)) - 2.0*(UpZ*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)))/(2.0*sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)),2.0) + pow((- DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)),2.0))));
	yy[5][4] = -(2.0*(UpZ*cos(a)*sin(b) - UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*sin(a)*sin(b))*(UpZ*cos(a)*cos(b) - z - DownZ + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)) + 2.0*(UpZ*cos(a)*cos(b)*cos(c) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(c)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*cos(c)*sin(a))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)) + 2.0*(UpZ*cos(a)*cos(b)*sin(c) + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b)*sin(c) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)*sin(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)))/(2.0*sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)),2.0) + pow((UpZ*cos(a)*cos(b) - z - DownZ + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)),2.0))));
	yy[5][5] = -(2.0*(UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c))*(x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)) - 2.0*(UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c))*(y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)))/(2.0*sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - DownR*
		sin(acosDownR/2.0 - pi/2.0) - UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*sin(c)),2.0) + pow((- DownZ - z + UpZ*cos(a)*cos(b) + UpR*
		cos(acosUpR/2.0 - pi/6)*sin(b) + UpR*
		sin(acosUpR/2.0 - pi/6)*cos(b)*sin(a)),2.0) + pow((x - UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - DownR*
		cos(acosDownR/2.0 - pi/2.0) + UpR*
		sin(acosUpR/2.0 - pi/6)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + UpR*
		cos(acosUpR/2.0 - pi/6)*cos(b)*cos(c)),2.0))));

	double mu = 0.0002;
	for (int i = 0;i < N;++i)
	{
		yy[i][i] += mu;
	}

}

void inv_jacobian(double yy[N][N], double inv[N][N])
{
	double aug[N][N2];
	double L = 0;
	int i, j, k;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			aug[i][j] = yy[i][j];
		}
		for (j = N; j < N2; j++)
		{
			aug[i][j] = (j == i + N) ? 1 : 0;
		}
	}
	for (i = 0; i < N; i++)
	{
		for (k = i + 1; k < N; k++)
		{
			L = -aug[k][i] / aug[i][i];
			for (j = i; j < N2; j++)
			{
				aug[k][j] = aug[k][j] + L * aug[i][j];
			}
		}
	}
	for (i = N - 1; i > 0; i--)
	{
		for (k = i - 1; k >= 0; k--)
		{
			L = -aug[k][i] / aug[i][i];
			for (j = N2 - 1; j >= 0; j--)
			{
				aug[k][j] = aug[k][j] + L * aug[i][j];
			}
		}
	}
	for (i = N - 1; i > 0; i--)
	{
		for (j = N2 - 1; j >= 0; j--)
		{
			aug[i][j] = aug[i][j] / aug[i][i];
		}
	}
	for (i = 0; i < N; i++)
	{
		for (j = N; j < N2; j++)
		{
			inv[i][j - N] = aug[i][j];
		}
	}

}

void newdundiedai(double x0[N], double inv[N][N], double y0[N], double x1[N])
{
	int i, j;
	double sum = 0;
	for (i = 0; i < N; i++)
	{
		sum = 0;
		for (j = 0; j < N; j++)
		{
			sum += inv[i][j] * y0[j];
		}
		x1[i] = x0[i] - sum;
	}
}

double result[N] = { 0, 0, 0, 0, 0, 0 };

double* ForwardKinematics(double dLen1, double dLen2, double dLen3, double dLen4, double dLen5, double dLen6, 
						  double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
						  double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom) 
{
	UpZ = planeAboveHingeLength;
	DownZ = planeAboveBottomLength;
	UpR = circleTopRadius;
	DownR = circleBottomRadius;
	Dis = (distanceBetweenHingeTop + distanceBetweenHingeBottom) / 2.0;
	double x0[N] = { 0, 0, 0, 0, 0, 0 };
	double x1[N] = { 0, 0, 0, 0, 0, 0 };
	double y0[N] = { 0, 0, 0, 0, 0, 0 };
	double jacobian[N][N] = { { 0 },{ 0 },{ 0 },{ 0 },{ 0 },{ 0 } };
	double invjacobian[N][N] = { { 0 },{ 0 },{ 0 },{ 0 },{ 0 },{ 0 } };
	double errornorm = 0;
	int i, j, iter = 0;
	do
	{
		iter++;
		ff(x0, y0, dLen1, dLen2, dLen3, dLen4, dLen5, dLen6);
		ffjacobian(x0, jacobian);
		inv_jacobian(jacobian, invjacobian);
		newdundiedai(x0, invjacobian, y0, x1);
		errornorm = 0;
		for (i = 0; i < N; i++)
		{
			errornorm += fabs(x1[i] - x0[i]);
		}
		if (errornorm < EPSILON)
		{
			break;
		}
		memcpy(x0, x1, sizeof(double) * N);	
	} while (iter < ITER_MAX);
	memcpy(result, x1, sizeof(double) * N);
	return result;
}


