
#include "stdafx.h"

#include <memory.h>

#include "platform.h"
#include "filters.h"
#include "MatrixOpreate.h"
#include "NewtonSolve.h"

Platform::Platform()
{
	PlatformParaInit();
    PositonsInit();
    AxisInit();
}

Platform::~Platform()
{

}

void Platform::SetPlatformPara(double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
							   double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom)
{
	// 单位：mm
	PlaneAboveHingeLength = planeAboveHingeLength;
	PlaneAboveBottomLength = planeAboveBottomLength;
	CircleTopRadius = circleTopRadius;
	CircleBottomRadius = circleBottomRadius;
	DistanceBetweenHingeTop = distanceBetweenHingeTop;
	distanceBetweenHingeBottom = distanceBetweenHingeBottom;
	PositonsInit();
	AxisInit();
}

double* Platform::Control(double x, double y, double z, double roll, double yaw, double pitch)
{
    BuildConversionMatrix(yaw, roll, pitch);
    for (auto i = 0; i < AXIS_COUNT; ++i)
    {
        double positionUforOArray[3];
        auto positionU = HingePositions[i];
        auto deltaPoint = Point3D(x, y, z);
        MatrixMultiplyVector(ConversionMatrix, positionU.ToArray(), positionUforOArray);
        auto positionUforO = Point3D(positionUforOArray);
        auto nowPosition = ZeroPoint + deltaPoint + positionUforO - BottomPositions[i];
        AxisDeltaLength[i] = nowPosition.Abs() - AxisInitLength[i];
    }
    return AxisDeltaLength;
}

void Platform::PlatformParaInit()
{
	// 单位：mm
	PlaneAboveHingeLength = 100;
	PlaneAboveBottomLength = 700;
	CircleTopRadius = 680;
	CircleBottomRadius = 840;
	DistanceBetweenHingeTop = 190;
	DistanceBetweenHingeBottom = 190;
}

void Platform::AxisInit()
{
    for (auto i = 0; i < AXIS_COUNT; ++i)
    {
        AxisInitLength[i] = HingePositions[i].Distance(BottomPositions[i]);
        AxisDeltaLength[i] = 0;
    }
}

void Platform::PositonsInit()
{
	memset(this->lengths, INIT_LENGTH, sizeof(double) * LENGTH_COUNT);
	memset(this->poses, INIT_POSE, sizeof(double) * POSE_COUNT);

	double CIRCLE_TOP_RADIUS = CircleTopRadius;
	double CIRCLE_BOTTOM_RADIUS = CircleBottomRadius;
	double DIS_BETWEEN_HINGE_TOP = DistanceBetweenHingeTop;
	double DIS_BETWEEN_HINGE_BOTTOM = DistanceBetweenHingeBottom;
	double PLANE_ABOVE_HINGE_LENGTH = PlaneAboveHingeLength;
	double PLANE_ABOVE_BOTTOM_LENGTH = PlaneAboveBottomLength;

    for (auto i = 0;i < AXIS_COUNT;++i)
    {
		HingePositions[i].Z = -PLANE_ABOVE_HINGE_LENGTH;
    }
	for (auto i = 0;i < AXIS_COUNT;++i)
	{
		BottomPositions[i].Z = -PLANE_ABOVE_BOTTOM_LENGTH;
	}

    auto top_between_hinge_rad = CosineTheorem(DIS_BETWEEN_HINGE_TOP, CIRCLE_TOP_RADIUS, CIRCLE_TOP_RADIUS) / 2.0;
    auto bottom_between_hinge_rad = CosineTheorem(DIS_BETWEEN_HINGE_BOTTOM, CIRCLE_BOTTOM_RADIUS, CIRCLE_BOTTOM_RADIUS) / 2.0;
    auto rad30 = RAD_30;
    auto rad60 = RAD_60;

    BottomPositions[0].X = CIRCLE_BOTTOM_RADIUS * cos(bottom_between_hinge_rad);
    BottomPositions[0].Y = CIRCLE_BOTTOM_RADIUS * sin(bottom_between_hinge_rad);

    BottomPositions[1].X = -CIRCLE_BOTTOM_RADIUS * sin(rad30 - bottom_between_hinge_rad);
    BottomPositions[1].Y = CIRCLE_BOTTOM_RADIUS * cos(rad30 - bottom_between_hinge_rad);

    BottomPositions[2].X = -CIRCLE_BOTTOM_RADIUS * sin(rad30 + bottom_between_hinge_rad);
    BottomPositions[2].Y = CIRCLE_BOTTOM_RADIUS * cos(rad30 + bottom_between_hinge_rad);

    BottomPositions[3].X = BottomPositions[2].X;
    BottomPositions[3].Y = -BottomPositions[2].Y;

    BottomPositions[4].X = BottomPositions[1].X;
    BottomPositions[4].Y = -BottomPositions[1].Y;

    BottomPositions[5].X = BottomPositions[0].X;
    BottomPositions[5].Y = -BottomPositions[0].Y;

    HingePositions[0].X = CIRCLE_TOP_RADIUS * cos(rad60 - top_between_hinge_rad);
    HingePositions[0].Y = CIRCLE_TOP_RADIUS * sin(rad60 - top_between_hinge_rad);

    HingePositions[1].X = CIRCLE_TOP_RADIUS * cos(rad60 + top_between_hinge_rad);
    HingePositions[1].Y = CIRCLE_TOP_RADIUS * sin(rad60 + top_between_hinge_rad);

    HingePositions[2].X = -CIRCLE_TOP_RADIUS * cos(top_between_hinge_rad);
    HingePositions[2].Y = CIRCLE_TOP_RADIUS * sin(top_between_hinge_rad);

    HingePositions[3].X = HingePositions[2].X;
    HingePositions[3].Y = -HingePositions[2].Y;

    HingePositions[4].X = HingePositions[1].X;
    HingePositions[4].Y = -HingePositions[1].Y;

    HingePositions[5].X = HingePositions[0].X;
    HingePositions[5].Y = -HingePositions[0].Y;
}

void Platform::BuildTsMatrix(double yaw, double roll, double pitch)
{
	double y = yaw;    //psi
	double a = roll;   //theta
	double b = pitch;  //phi

	TsMatrix[0][0] = 1;
	TsMatrix[0][1] = sin(b) * tan(a);
	TsMatrix[0][2] = cos(b) * tan(a);

	TsMatrix[1][0] = 0;
	TsMatrix[1][1] = cos(b);
	TsMatrix[1][2] = -sin(b);

	TsMatrix[2][0] = 0;
	TsMatrix[2][1] = sin(b) * 1.0 / cos(a);
	TsMatrix[2][2] = cos(b) * 1.0 / cos(a);

}

void Platform::BuildLsMatrix(double yaw, double roll, double pitch)
{
	double y = yaw;    //psi
	double a = roll;   //theta
	double b = pitch;  //phi

	LsMatrix[0][0] = cos(a) * cos(y);
	LsMatrix[0][1] = sin(b) * cos(a) * cos(y) - cos(b) * sin(y);
	LsMatrix[0][2] = cos(b) * sin(a) * cos(y) + sin(b) * sin(y);

	LsMatrix[1][0] = cos(a) * sin(y);	
	LsMatrix[1][1] = sin(b) * sin(a) * sin(y) + cos(b) * cos(y);
	LsMatrix[1][2] = cos(b) * sin(a) * sin(y) - sin(b) * cos(y);

	LsMatrix[2][0] = -sin(a);
	LsMatrix[2][1] = sin(b) * cos(a);
	LsMatrix[2][2] = cos(b) * cos(a);
}

void Platform::BuildConversionMatrix(double yaw, double roll, double pitch)
{
	double y = yaw;    //psi
	double a = roll;   //theta
	double b = pitch;  //phi

    ConversionMatrix[0][0] = cos(y) * cos(b);
    ConversionMatrix[0][1] = -sin(y) * cos(a) + cos(y) * sin(b) * sin(a);
    ConversionMatrix[0][2] = sin(y) * sin(a) + cos(y) * sin(b) * cos(a);

    ConversionMatrix[1][0] = sin(y) * cos(b);
    ConversionMatrix[1][1] = cos(y) * cos(a) + sin(y) * sin(b) * sin(a);
    ConversionMatrix[1][2] = -cos(y) * sin(a) + sin(y) * sin(b) * cos(a);

    ConversionMatrix[2][0] = -sin(b);
    ConversionMatrix[2][1] = cos(b) * sin(a);
    ConversionMatrix[2][2] = cos(b) * cos(a);
}

/*
用余弦定理求三角形边a对应的夹角
*/
double Platform::CosineTheorem(double a, double b, double c)
{
    return acos((b * b + c * c - a * a) / (2 * b * c));
}
/*
运动学正解，由位姿信息求解杆长
*/
double * Platform::FromLengthToPose(double * lengths)
{
	/*Matlab Code as follows*/
	/*
	解六元非线性方程组
	*/
	memcpy(this->lengths, lengths, sizeof(double) * LENGTH_COUNT);
	auto result = ForwardKinematics(this->lengths[0], this->lengths[1], this->lengths[2], this->lengths[3], this->lengths[4], this->lengths[5], 
		PlaneAboveHingeLength, PlaneAboveBottomLength, CircleTopRadius, CircleBottomRadius, DistanceBetweenHingeTop, DistanceBetweenHingeBottom);
	memcpy(this->poses, result, sizeof(double) * LENGTH_COUNT);
	return this->poses;
}

/*
* 洗出算法又称体感模拟算法,算法引入了经典滤波算法、惯性坐标转换、限制环节等，
*/
double * Platform::WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch,
						 double xacc, double yacc, double zacc, double rollSpeed, double yawSpeed, double pitchSpeed)
{
	static double acc_scale = 1.0;
	static double angleSpd_scale = 1.0;
	static double coor_turn_gain = 0.1;
	static AccHighPassFilter accHighPassFilters[ACC_NUM];
	static AccIntZtrans accIntZtrans[ACC_NUM];
	static AccLowPassFilter accLowPassFilter[ACC_NUM];
	static AngleSpeedHighPassFilterAndInt angleHpfAndInt[ANGLE_SPEED_NUM];
	BuildLsMatrix(yaw, roll, pitch);
	BuildTsMatrix(yaw, roll, pitch);
	//double y = yaw;    //psi
	double a = roll;   //theta
	double b = pitch;  //phi
	double fAA[ACC_NUM] = {xacc * acc_scale, yacc * acc_scale, zacc * acc_scale};
	double wAA[ANGLE_SPEED_NUM] = {a * angleSpd_scale, b * angleSpd_scale, yaw * angleSpd_scale};
	double f2[ACC_NUM];
	double flow[ACC_NUM];
	double beta2[ANGLE_SPEED_NUM];
	double betahigh[ANGLE_SPEED_NUM];
	double betalow[ANGLE_SPEED_NUM];
	double ahigh[ACC_NUM];
	double betaS[ANGLE_SPEED_NUM];
	MatrixMultiplyVector(LsMatrix, fAA, f2);  
	MatrixMultiplyVector(TsMatrix, wAA, beta2);
	double a2[3] = {0};
	memcpy(a2, f2, sizeof(double) * ACC_NUM);
#if IS_ADD_EARTH_G
	a2[ACC_NUM - 1] = f2[ACC_NUM - 1] + EARTH_G;
#endif
	for (int i = 0; i < ACC_NUM; ++i)
	{
		ahigh[i] = accHighPassFilters[i].Update(a2[i]);
		poses[i] = accIntZtrans[i].Update(ahigh[i]) / 1000.0; //变成mm
		flow[i] = accLowPassFilter[i].Update(fAA[i]);
		betalow[i] = LIMITER(flow[i] * coor_turn_gain, -ANGLE_VEL_UP_RANGE, +ANGLE_VEL_UP_RANGE);
		betaS[i] = betalow[i] + angleHpfAndInt->Update(beta2[i]);
	}
	memcpy(poses + ACC_NUM, betaS, sizeof(double) * ANGLE_SPEED_NUM);
	return poses;
}
