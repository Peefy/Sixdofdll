// Sixdofdll2010.Test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "../Sixdofdll2010/Sixdofdll2010.h"
#include "function_value_test.h"

#define DIS_PER_R 5.36
#define PULSE_PER_R 1280000
#define LENGTH_TO_PULSE_SCALE PULSE_PER_R / DIS_PER_R

using namespace std;

void PrintControlData(double * length)
{
	if (length == NULL)
		return;
	printf("%2.3f %2.3f %2.3f %2.3f %2.3f %2.3f\r\n", length[0] * LENGTH_TO_PULSE_SCALE, length[1], length[2],
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
	int datacount = 8000;
	WashOutFiltering(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Sleep(47);
}

int main()
{
	PrintControlData(Control(0, 0, 0, 0, 0, 0));
	PrintControlData(Control(0, 0, 200, 0, 0, 0));
	PrintControlData(Control(0, 0, -200, 0, 0, 0));
	PrintControlData(Control(200, 0, 0, 0, 0, 0));
	PrintControlData(Control(-200, 0, 0, 0, 0, 0));
	PrintControlData(Control(0, 200, 0, 0, 0, 0));
	PrintControlData(Control(0, -200, 0, 0, 0, 0));
	PrintPosition(GetTopPosition(0));
	PrintControlData(Control(0, 0, 0, 10, 0, 0));
	double lengths[6] = {10,20,20,20,20,20};
	void WashoutTest();
	
	auto poses = FromLengthToPose(lengths);
	for (int i = 0;i < 6;++i)
	{
		printf("the %d poses is %f \r\n", i, poses[i]);
	}
	double i = 0;
	while (i < 3)
	{
		i += 0.005;
		auto angle = 30 * sin(i);
		PrintControlData(Control(0, 0, angle, 0, 0, 0));
		Sleep(5);
	}
	getchar();
	return 0;
}


