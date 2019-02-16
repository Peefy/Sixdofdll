#pragma once

#include <math.h>
#include "config.h"

// �˶�ѧ����������

// �����ɶ�ƽ̨�׵�����
#define LENGTH_COUNT 6
// �����ɶ�ƽ̨�໥���������ɶ�����
#define POSE_COUNT 6
// ��ʼ����Ϊ0
#define INIT_LENGTH 0
// ��ʼ��̬Ϊ0
#define INIT_POSE 0
// �������ά��
#define MATRIX_DIMENSION 3

// ϴ���㷨

// ���ٶȵ�����
#define ACC_NUM          3       

// ���ٶȵ�����
#define ANGLE_SPEED_NUM  3       

// �������ٶ�ȡ9.8 m/s^2
#define EARTH_G          9.8     

// �����������Ƿ����������ٶ�
#define IS_ADD_EARTH_G   0

// ���ٶ��޷�
#define ACC_UP_RANGE        (0.02 * EARTH_G) //��λ m/s^2 

// ���ٶ����Ʒ���
#define ANGLE_VEL_UP_RANGE  3 // m/s^2

// ��ά����
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
	// ��ƽ��������Ϊ�ѿ���ֱ������ϵԭ������
    Point3D ZeroPoint;
	// �����Ͻ�֧������
    Point3D HingePositions[AXIS_COUNT];
	// �����½�֧������
    Point3D BottomPositions[AXIS_COUNT];
	// ��������ϵת��������ϵ�Ĵ��ݾ���
    double ConversionMatrix[MATRIX_DIMENSION][MATRIX_DIMENSION];
	// ϴ���㷨���ٶȱ任����
	double TsMatrix[MATRIX_DIMENSION][MATRIX_DIMENSION];
	// ϴ���㷨���ٶȱ任����
	double LsMatrix[MATRIX_DIMENSION][MATRIX_DIMENSION];
	// ƽ̨�����׵ĳ�ʼ����
    double AxisInitLength[AXIS_COUNT];
	// ƽ̨�����׵���������
    double AxisDeltaLength[AXIS_COUNT];
	void SetPlatformPara(double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
		double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom);
	void SetWashOutFilterPara(double hpfAccWn, double lpfAccWn, double hpfAngleSpdWn);
    double* Control(double x, double y, double z, double roll, double yaw, double pitch);
    double* FromLengthToPose(double * lengths);
	double* WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch,
		double xacc, double yacc, double zacc, double rollSpeed, double yawSpeed, double pitchSpeed);
private:
	double lengths[LENGTH_COUNT];
	// [0]:x [1]:y [2]:z [3]:roll [4]:pitch [5]:yaw
	double poses[POSE_COUNT];
	// ��Բ�뾶
	double CircleTopRadius;
	// ��Բ�뾶
	double CircleBottomRadius;
	// ÿ���Ͻ�֧�������ļ��
	double DistanceBetweenHingeTop;
	// ÿ���½�֧�������ļ��
	double DistanceBetweenHingeBottom;
	// ��ƽ�浽�Ͻ�֧�����ĵĴ�ֱ����
	double PlaneAboveHingeLength;
	// ��ƽ�浽�½�֧�����ĵĴ�ֱ����
	double PlaneAboveBottomLength;
	double HpfAccWn;
	double LpfAccWn;
	double HpfAngleSpdWn;
	void PlatformParaInit();
	void WashOutFilterParaInit();
    void AxisInit();
    void PositonsInit();
    void BuildConversionMatrix(double yaw, double roll, double pitch);
	void BuildTsMatrix(double yaw, double roll, double pitch);
	void BuildLsMatrix(double yaw, double roll, double pitch);
    double CosineTheorem(double a, double b, double c);
};

