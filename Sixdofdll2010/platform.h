#pragma once

#include <math.h>
#include "config.h"
#include "filters.h"

// 运动学正、反解算

// 六自由度平台缸的数量
#define LENGTH_COUNT 6
// 六自由度平台相互独立的自由度数量
#define POSE_COUNT 6
// 初始长度为0
#define INIT_LENGTH 0
// 初始姿态为0
#define INIT_POSE 0
// 方阵矩阵维数
#define MATRIX_DIMENSION 3

// 洗出算法

// 加速度的数量
#define ACC_NUM          3       

// 角速度的数量
#define ANGLE_SPEED_NUM  3       

// 重力加速度取9.8 m/s^2
#define EARTH_G          9.8     

// 比力传感器是否考虑重力加速度
#define IS_ADD_EARTH_G   0

// 洗出算法是否使用传递矩阵
#define IS_USE_TRANS_MATRIX 1

// 是否增加协调转弯
#define IS_ADD_COOR_TURN_GAIN 1

// 加速度限幅
#define ACC_UP_RANGE        (0.02 * EARTH_G) //单位 m/s^2 

// 角速度限制幅度
#define ANGLE_VEL_UP_RANGE 6.5 // m/s^2

// 洗出算法滤波器阶数+1
#define WASHOUT_FILTER_ORDER_PLUS_ONE 3

// 三维坐标
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

    explicit Point3D(double* val)
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
	// 上平面中心作为笛卡尔直角坐标系原点坐标
    Point3D ZeroPoint;
	// 六个上铰支座坐标
    Point3D HingePositions[AXIS_COUNT];
	// 六个下铰支座坐标
    Point3D BottomPositions[AXIS_COUNT];
	// 地球坐标系转惯性坐标系的传递矩阵
    double ConversionMatrix[MATRIX_DIMENSION][MATRIX_DIMENSION];
	// 洗出算法加速度变换矩阵
	double TsMatrix[MATRIX_DIMENSION][MATRIX_DIMENSION];
	// 洗出算法角速度变换矩阵
	double LsMatrix[MATRIX_DIMENSION][MATRIX_DIMENSION];
	// 平台六个缸的初始长度
    double AxisInitLength[AXIS_COUNT];
	// 平台六个缸的增量长度
    double AxisDeltaLength[AXIS_COUNT];
	// 设置平台机械参数
	void SetPlatformPara(double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
		double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom);
	// 设置洗出滤波算法参数
	void SetWashOutFilterPara(double hpfAccWn, double lpfAccWn, double hpfAngleSpdWn, double sampleTime);
	// 运动学反解
    double* Control(double x, double y, double z, double roll, double yaw, double pitch);
	// 运动学正解
    double* FromLengthToPose(double * lengths);
	// 洗出滤波算法
	double* WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch,
		double xacc, double yacc, double zacc, double rollSpeed, double yawSpeed, double pitchSpeed);
private:
	double lengths[LENGTH_COUNT];
	// [0]:x [1]:y [2]:z [3]:roll [4]:pitch [5]:yaw
	double poses[POSE_COUNT];
	// 上圆半径
	double CircleTopRadius;
	// 下圆半径
	double CircleBottomRadius;
	// 每对上铰支座的中心间距
	double DistanceBetweenHingeTop;
	// 每对下铰支座的中心间距
	double DistanceBetweenHingeBottom;
	// 上平面到上铰支座中心的垂直距离
	double PlaneAboveHingeLength;
	// 上平面到下铰支座中心的垂直距离
	double PlaneAboveBottomLength;
	// 加速度高通滤波器截止频率
	double HpfAccWn;
	// 加速度低通滤波器截止频率
	double LpfAccWn;
	// 角速度高通滤波器截止频率
	double HpfAngleSpdWn;
	// 采样周期
	double SampleTime;
	void PlatformParaInit();
	void WashOutFilterParaInit();
    void AxisInit();
    void PositonsInit();
    void BuildConversionMatrix(double yaw, double roll, double pitch);
	void BuildTsMatrix(double yaw, double roll, double pitch);
	void BuildLsMatrix(double yaw, double roll, double pitch);
    double CosineTheorem(double a, double b, double c);
protected:
	// 洗出算法滤波器
	AccHighPassFilter accHighPassFilters[ACC_NUM];
	AccIntZtrans accIntZtrans[ACC_NUM];
	AccLowPassFilter accLowPassFilter[ACC_NUM];
	AngleSpeedHighPassFilterAndInt angleHpfAndInt[ANGLE_SPEED_NUM];
	double nums[WASHOUT_FILTER_ORDER_PLUS_ONE];
	double dens[WASHOUT_FILTER_ORDER_PLUS_ONE];
};

