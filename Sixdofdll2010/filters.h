
#ifndef _FILTERS_H_
#define _FILTERS_H_

#include <queue>
#include <vector>
#include <deque>

using namespace std;

// 限幅
#define LIMITER(x, min, max)     (((x)<(min) ? (min) : ( (x)>(max) ? (max):(x) )))

// 滤波器最大阶数
#define MAX_ORDER 10
// 数字洗出滤波算法采样周期
#define INIT_DT   0.047

// 标识数据输入
#define PARA_IN
// 标识数据输出
#define PARA_OUT

// 斜率限制器
double RateLimiter(double dataIn, double risingSlewRate, double fallingSlewRate, double dT);

// z变换传递函数
class Ztrans
{
public:
	Ztrans(double dT, int order, double* nums, double* dens);
	Ztrans(double dT, int order);
	~Ztrans() {};
	int Order;
	// 设置采样周期
	void SetSampleTime(double sampleTime);
	// 输入经传递函数得到输出
	double Update(double now);
	// 设置z传递函数的分子分母系数
	void SetNumsAndDensZtrans(double* nums, double* dens);
	// 设置s传递函数的分子分母系数
	void SetNumsAndDensLaplace(double* nums, double* dens, double fs);
	// s域到z域的双线性变换
	void Bilinear(PARA_IN double* b, PARA_IN double* a, double fs, int dimensions, 
		PARA_OUT double* bprimes, PARA_OUT double* aprimes);
private:
	double dT;
	double fs;
	deque<double> input;
	deque<double> output;
	double inner_nums[MAX_ORDER];
	double inner_dens[MAX_ORDER];
};

// 加速度高通滤波器
class AccHighPassFilter : public Ztrans
{
public:
	AccHighPassFilter(double dT = INIT_DT); 
	~AccHighPassFilter() {};

private:

};

// 加速度积分成位移
class AccIntZtrans : public Ztrans
{
public:
	AccIntZtrans(double dT = INIT_DT); 
	~AccIntZtrans() {};

private:

};

// 加速度滤波器
class AccLowPassFilter : public Ztrans
{
public:
	AccLowPassFilter(double dT = INIT_DT); 
	~AccLowPassFilter() {};

private:

};

// 角速度高通滤波器
class AngleSpeedHighPassFilterAndInt : public Ztrans
{
public:
	AngleSpeedHighPassFilterAndInt(double dT = INIT_DT); 
	~AngleSpeedHighPassFilterAndInt() {};

private:

};



#endif // !_FILTERS_H_

