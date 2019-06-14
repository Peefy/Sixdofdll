// Sixdofdll2010.Test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "../Sixdofdll2010/Sixdofdll2010.h"
#include "function_value_test.h"

#include <fstream>

#define PlaneAboveHingeLength       225.0
#define PlaneAboveBottomLength      2105.0
#define CircleTopRadius             880.7
#define CircleBottomRadius          1519.0
#define DistanceBetweenHingeTop     200.0
#define DistanceBetweenHingeBottom  300.0

#define FILENAME "illusiondata.txt"

using namespace std;

void PrintControlData(double * length)
{
	if (length == NULL)
		return;
	printf("%2.3f %2.3f %2.3f %2.3f %2.3f %2.3f\r\n", length[0], length[1], length[2],
		length[3], length[4], length[5]);
}

void PrintPosition(double* pos)
{
	if (pos == NULL)
		return;
	printf("X:%f;y:%f;z:%f\r\n", pos[0], pos[1], pos[2]);
}

void WashoutTest()
{
	int datacount = 200;
	double x, y, z, xacc, yacc, zacc, rollSpeed, pitchSpeed, yawSpeed, roll, pitch, yaw;
	ifstream fin(FILENAME);
	double poses[6] = {0};
	SetWashOutFilterPara(1.0, 0.5, 0.5, 0.047);
	for (int i = 0;i < datacount;++i)
	{
		fin >> xacc;
		fin >> yacc;
		fin >> zacc;
		fin >> rollSpeed;
		fin >> pitchSpeed;
		fin >> yawSpeed;
		fin >> roll;
		fin >> yaw;
		fin >> pitch;
		auto* tmp = WashOutFiltering(poses[0], poses[1], poses[2], poses[3], poses[5], poses[4], 
			xacc, yacc, zacc, rollSpeed, yawSpeed, pitchSpeed);
		memcpy(poses, tmp, sizeof(double) * 6);
		printf("the %d poses is %f \r\n", i, poses[0]);
		Sleep(47);
	}
	fin.close();
}

int main()
{
	SetPlatformPara(PlaneAboveHingeLength, PlaneAboveBottomLength, 
		CircleTopRadius, CircleBottomRadius, DistanceBetweenHingeTop,
		DistanceBetweenHingeBottom);
	PrintControlData(Control(0, 0, 0, 0, 0, 0));
	PrintControlData(Control(0, 0, 200, 0, 0, 0));
	PrintControlData(Control(0, 0, -200, 0, 0, 0));
	PrintControlData(Control(200, 0, 0, 0, 0, 0));
	PrintControlData(Control(-200, 0, 0, 0, 0, 0));
	PrintControlData(Control(0, 200, 0, 0, 0, 0));
	PrintControlData(Control(0, -200, 0, 0, 0, 0));
	PrintPosition(GetTopPosition(0));
	PrintControlData(Control(0, 0, 0, 10, 0, 0));
	double lengths[6] = {30, 40, 40, 40, 40, 40};
		
	auto poses = FromLengthToPose(lengths);
	for (int i = 0;i < 6;++i)
	{
		printf("the %d poses is %f \r\n", i, poses[i]);
	}
	double i = 0; 
	while (i < 3)
	{
		i += 0.003;
		auto angle = 250 * sin(i);
		PrintControlData(Control(0, 0, angle, 0, 0, 0));
		Sleep(5);
	}
	WashoutTest();
	getchar();
	return 0;
}


