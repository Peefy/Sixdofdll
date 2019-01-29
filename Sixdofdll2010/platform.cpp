
#include "stdafx.h"

#include <memory.h>

#include "platform.h"
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
	// ��λ��mm
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
	// ��λ��mm
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
	double y = yaw;    //pothi
	double a = roll;   //theta
	double b = pitch;  //phi

	TsMatrix[0][0] = 1;
	TsMatrix[0][1] = -sin(y) * cos(a) + cos(y) * sin(b) * sin(a);
	TsMatrix[0][2] = sin(y) * sin(a) + cos(y) * sin(b) * cos(a);

	TsMatrix[1][0] = 0;
	TsMatrix[1][1] = cos(y) * cos(a) + sin(y) * sin(b) * sin(a);
	TsMatrix[1][2] = -cos(y) * sin(a) + sin(y) * sin(b) * cos(a);

	TsMatrix[2][0] = 0;
	TsMatrix[2][1] = cos(b) * sin(a);
	TsMatrix[2][2] = cos(b) * cos(a);

}

void Platform::BuildConversionMatrix(double yaw, double roll, double pitch)
{
	double y = yaw;    //pothi
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
�����Ҷ����������α�a��Ӧ�ļн�
*/
double Platform::CosineTheorem(double a, double b, double c)
{
    return acos((b * b + c * c - a * a) / (2 * b * c));
}
/*
�˶�ѧ���⣬��λ����Ϣ���˳�
*/
double * Platform::FromLengthToPose(double * lengths)
{
	/*Matlab Code as follows*/
	/*
	����Ԫ�����Է�����
	*/
	memcpy(this->lengths, lengths, sizeof(double) * LENGTH_COUNT);
	auto result = ForwardKinematics(this->lengths[0], this->lengths[1], this->lengths[2], this->lengths[3], this->lengths[4], this->lengths[5], 
		PlaneAboveHingeLength, PlaneAboveBottomLength, CircleTopRadius, CircleBottomRadius, DistanceBetweenHingeTop, DistanceBetweenHingeBottom);
	memcpy(this->poses, result, sizeof(double) * LENGTH_COUNT);
	return this->poses;
}

double * Platform::WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch)
{
    return nullptr;
}

double * Platform::WashOutFiltering(double x, double y, double z, double roll, double yaw, double pitch,
						 double xacc, double yacc, double zacc, double rollSpeed, double yawSpeed, double pitchSpeed)
{
	return nullptr;
}