
#ifndef _FILTERS_H_
#define _FILTERS_H_

#define MAX_ORDER 50

#include <queue>
#include <vector>
#include <deque>

using namespace std;

class Filter
{
public:
	Filter(int order, double* nums, double* dens);
	~Filter();
	int Order;
	double Update(double now);
private:
	deque<double> input;
	deque<double> output;
	double inner_nums[MAX_ORDER];
	double inner_dens[MAX_ORDER];
};

#endif // !_FILTERS_H_

