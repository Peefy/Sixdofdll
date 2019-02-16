
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

#define PARA_IN
#define PARA_OUT

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
	void SetNumsAndDensZtrans(double* nums, double* dens);
	void SetNumsAndDensLaplace(double* nums, double* dens, double fs);
	void Bilinear(PARA_IN double* b, PARA_IN double* a, double fs, int dimensions, 
		PARA_OUT double* bprimes, PARA_OUT double* aprimes);
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

