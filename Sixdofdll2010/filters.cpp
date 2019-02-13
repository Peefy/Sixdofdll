
#include "filters.h"

Filter::Filter(int order, double* nums, double* dens)
{
	Order = order + 1;
	for (int i = 0; i < Order; i++)
	{
		input.push_front(0);
		output.push_front(0);
	}
	memcpy(inner_nums, nums, sizeof(double) * Order);
	memcpy(inner_dens, dens, sizeof(double) * Order);
}

Filter::~Filter()
{
}

double Filter::Update(double now)
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

