
#ifndef _FILTERS_H_
#define _FILTERS_H_

#include <queue>
#include <vector>
#include <deque>

using namespace std;

#define LIMITER(x, min, max)     (((x)<(min) ? (min) : ( (x)>(max) ? (max):(x) )))

#define MAX_ORDER 10
//#define INIT_DT   0.01
#define INIT_DT   0.047

double RateLimiter(double dataIn, double risingSlewRate, double fallingSlewRate, double dT);

class Ztrans
{
public:
	Ztrans(double dT, int order, double* nums, double* dens);
	Ztrans(double dT, int order);
	~Ztrans() {};
	double dT;
	int Order;
	double Update(double now);
	void SetNumsAndDens(double* nums, double* dens);
private:
	deque<double> input;
	deque<double> output;
	double inner_nums[MAX_ORDER];
	double inner_dens[MAX_ORDER];
};

class AccHighPassFilter : public Ztrans
{
public:
	AccHighPassFilter(double dT = INIT_DT); 
	~AccHighPassFilter() {};

private:

};

class AccIntZtrans : public Ztrans
{
public:
	AccIntZtrans(double dT = INIT_DT); 
	~AccIntZtrans() {};

private:

};

class AccLowPassFilter : public Ztrans
{
public:
	AccLowPassFilter(double dT = INIT_DT); 
	~AccLowPassFilter() {};

private:

};

class AngleSpeedHighPassFilterAndInt : public Ztrans
{
public:
	AngleSpeedHighPassFilterAndInt(double dT = INIT_DT); 
	~AngleSpeedHighPassFilterAndInt() {};

private:

};



#endif // !_FILTERS_H_
