
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
	for (int i = 0;i < Order - 1;++i)
	{
		out -= inner_dens[i + 1] * output[i];
	}
	for (int i = 0;i < Order;++i)
	{
		out += inner_nums[i] * input[i];
	}
	output.push_front(out);
	output.pop_back();
	return out;
}

AccHighPassFilter::AccHighPassFilter(double dT) : Ztrans(dT, 1)
{
	double nums[2] = {0, 0.01};
	double dens[2] = {1, -1};
	SetNumsAndDens(nums, dens);
}

AccLowPassFilter::AccLowPassFilter(double dT) : Ztrans(dT, 1)
{
	double nums[2] = {0, 0.01};
	double dens[2] = {1, -1};
	SetNumsAndDens(nums, dens);
}

AccIntZtrans::AccIntZtrans(double dT) : Ztrans(dT, 1)
{
	double nums[2] = {0, 0.01};
	double dens[2] = {1, -1};
	SetNumsAndDens(nums, dens);
}

AngleSpeedHighPassFilter::AngleSpeedHighPassFilter(double dT) : Ztrans(dT, 1)
{
	double nums[2] = {0, 0.01};
	double dens[2] = {1, -1};
	SetNumsAndDens(nums, dens);
}
