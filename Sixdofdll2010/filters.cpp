
#include "stdafx.h"
#include "filters.h"

Ztrans::Ztrans(double dT, int order)
{
	this->dT = dT;
	Order = order + 1;
	for (int i = 0; i < Order; i++)
	{
		input.push_front(0);
		output.push_front(0);
	}
}

Ztrans::Ztrans(double dT, int order, double* nums, double* dens)
{
	this->dT = dT;
	Order = order + 1;
	for (int i = 0; i < Order; i++)
	{
		input.push_front(0);
		output.push_front(0);
	}
	SetNumsAndDens(nums, dens);
}

void Ztrans::SetNumsAndDens(double* nums, double* dens)
{
	if (nums == nullptr || dens == nullptr)
	{
		return;
	}
	memcpy(inner_nums, nums, sizeof(double) * Order);
	memcpy(inner_dens, dens, sizeof(double) * Order);
}

double Ztrans::Update(double now)
{
	double out = 0;
	input.push_front(now);
	input.pop_back();
	for (int i = 0; i < Order - 1; ++i)
	{
		out -= inner_dens[i + 1] * output[i];
	}
	for (int i = 0; i < Order; ++i)
	{
		out += inner_nums[i] * input[i];
	}
	output.push_front(out);
	output.pop_back();
	return out;
}

AccHighPassFilter::AccHighPassFilter(double dT) : Ztrans(dT, 2)
{
	//double nums[3] = {1, -2, 0.9996};
	//double dens[3] = {1, -1.93, 0.9305};
	double nums[3] = {1, -1.886, 0.8856};
	double dens[3] = {1, -1.052, 0.2369};
	SetNumsAndDens(nums, dens);
}

AccIntZtrans::AccIntZtrans(double dT) : Ztrans(dT, 2)
{
	//double nums[3] = {0, 4.983e-5, 4.967e-5};
	//double dens[3] = {1, -1.99, 0.99};
	double nums[3] = {0, 0.01873, 0.01752};
	double dens[3] = {1, -1.819, 0.8187};
	SetNumsAndDens(nums, dens);
}

AccLowPassFilter::AccLowPassFilter(double dT) : Ztrans(dT, 2)
{
	//double nums[3] = {0, 0.0003073, 0.0003023};
	//double dens[3] = {1, -1.951, 0.9512};
	double nums[3] = {0, 0.0902, 0.06461};
	double dens[3] = {1, -1.213, 0.3679};
	SetNumsAndDens(nums, dens);
}

AngleSpeedHighPassFilterAndInt::AngleSpeedHighPassFilterAndInt(double dT) : Ztrans(dT, 2)
{
	//double nums[3] = {0, 0.009753, -0.009753};
	//double dens[3] = {1, -1.951, 0.9512};
	double nums[3] = {0, 0.1213, -0.1213};
	double dens[3] = {1, -1.213, 0.3679};
	SetNumsAndDens(nums, dens);
}

