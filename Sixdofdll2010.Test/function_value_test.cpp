
#include "stdafx.h"
#include "function_value_test.h"

const double pi = 3.141592653589793;

double BuildFunctionValue(double x, double y, double z, double a, double b, double c)
{
    double UpZ = 100;
    double DownZ = 700;
    double UpR = 680;
    double DownR = 840;
    double Dis = 190;
    double f = -(2 * UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2.0 * x + 2.0 * DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) /
        2.0 - pi / 2.0) + 2.0 * UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - 2.0 *
        UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*cos(c)) / (2.0 * sqrt((pow((y + UpZ*(cos(c)*sin(a) - cos(a)*sin(b)*
        sin(c)) + DownR*sin(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 2.0) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) /
        2 - pi / 6.0)*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + UpR*cos(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*sin(c)), 2) +
        pow((-x + UpZ*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + DownR*cos(acos((Dis * Dis - 2.0 * DownR * DownR) / (2.0 * DownR * DownR)) / 2.0 - pi / 2.0) + UpR*
        sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) - UpR*cos(acos((Dis * Dis -
        2.0 * UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*cos(c)), 2) + pow((DownZ + z - UpZ*cos(a)*cos(b) - UpR*cos(acos((Dis * Dis - 2.0 *
        UpR * UpR) / (2.0 * UpR * UpR)) / 2.0 - pi / 6.0)*sin(b) + UpR*sin(acos((Dis * Dis - 2.0 * UpR * UpR) / (2 * UpR * UpR)) / 2.0 - pi / 6.0)*cos(b)*
        sin(a)), 2))));
    return f;
}
