#pragma once

#include <math.h>
#include "config.h"

#define LENGTH_COUNT 6
#define POSE_COUNT 6

#define INIT_LENGTH 0
#define INIT_POSE 0

class Point3D
{
public:
    double X;
    double Y;
    double Z;

    Point3D() {}

    Point3D(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    Point3D(double* val)
    {
        if (val == nullptr)
            return;
        X = val[0];
        Y = val[1];
        Z = val[2];
    }

    double* ToArray()
    {       
        val[0] = X;
        val[1] = Y;
        val[2] = Z;
        return val;
    }

    double Abs()
    {
        return sqrt(X * X + Y * Y + Z * Z);
    }

    double Distance(Point3D & right)
    {
        auto x = (X - right.X) * (X - right.X);
        auto y = (Y - right.Y) * (Y - right.Y);
        auto z = (Z - right.Z) * (Z - right.Z);
        return sqrt(x + y + z);
    }

    Point3D operator+(Point3D & right)
    {
        return Point3D(X + right.X, Y + right.Y, Z + right.Z);
    }

    Point3D operator-(Point3D & right)
    {
        return Point3D(X - right.X, Y - right.Y, Z - right.Z);
    }

private:
    double val[3];
};

class Platform
{
public:
    Platform();
    ~Platform();
    Point3D ZeroPoint;
    Point3D HingePositions[AXIS_COUNT];
    Point3D BottomPositions[AXIS_COUNT];
    double ConversionMatrix[3][3];
	double TsMatrix[3][3];
    double AxisInitLength[AXIS_COUNT];
    double AxisDeltaLength[AXIS_COUNT];
	void SetPlatformPara(double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
		double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom);
    double* Control(double x, double y, double z, double roll, double yaw, double pitch);
    double* FromLengthToPose(double * lengths);
    double* WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch);
	double* WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch,
		double xacc, double yacc, double zacc, double rollSpeed, double yawSpeed, double pitchSpeed);
private:
	double lengths[LENGTH_COUNT];
	double poses[POSE_COUNT];
	double CircleTopRadius;
	double CircleBottomRadius;
	double DistanceBetweenHingeTop;
	double DistanceBetweenHingeBottom;
	double PlaneAboveHingeLength;
	double PlaneAboveBottomLength;
	void PlatformParaInit();
    void AxisInit();
    void PositonsInit();
    void BuildConversionMatrix(double yaw, double roll, double pitch);
	void BuildTsMatrix(double yaw, double roll, double pitch);
    double CosineTheorem(double a, double b, double c);
};

