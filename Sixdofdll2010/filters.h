
#ifndef _FILTERS_H_
#define _FILTERS_H_

#include <queue>
#include <vector>
#include <deque>

using namespace std;

// �޷�
#define LIMITER(x, min, max)     (((x)<(min) ? (min) : ( (x)>(max) ? (max):(x) )))

// �˲���������
#define MAX_ORDER 10
// ����ϴ���˲��㷨��������
#define INIT_DT   0.047

// ��ʶ��������
#define PARA_IN
// ��ʶ�������
#define PARA_OUT

// б��������
double RateLimiter(double dataIn, double risingSlewRate, double fallingSlewRate, double dT);

// z�任���ݺ���
class Ztrans
{
public:
	Ztrans(double dT, int order, double* nums, double* dens);
	Ztrans(double dT, int order);
	~Ztrans() {};
	int Order;
	// ���ò�������
	void SetSampleTime(double sampleTime);
	// ���뾭���ݺ����õ����
	double Update(double now);
	// ����z���ݺ����ķ��ӷ�ĸϵ��
	void SetNumsAndDensZtrans(double* nums, double* dens);
	// ����s���ݺ����ķ��ӷ�ĸϵ��
	void SetNumsAndDensLaplace(double* nums, double* dens, double fs);
	// s��z���˫���Ա任
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

// ���ٶȸ�ͨ�˲���
class AccHighPassFilter : public Ztrans
{
public:
	AccHighPassFilter(double dT = INIT_DT); 
	~AccHighPassFilter() {};

private:

};

// ���ٶȻ��ֳ�λ��
class AccIntZtrans : public Ztrans
{
public:
	AccIntZtrans(double dT = INIT_DT); 
	~AccIntZtrans() {};

private:

};

// ���ٶ��˲���
class AccLowPassFilter : public Ztrans
{
public:
	AccLowPassFilter(double dT = INIT_DT); 
	~AccLowPassFilter() {};

private:

};

// ���ٶȸ�ͨ�˲���
class AngleSpeedHighPassFilterAndInt : public Ztrans
{
public:
	AngleSpeedHighPassFilterAndInt(double dT = INIT_DT); 
	~AngleSpeedHighPassFilterAndInt() {};

private:

};



#endif // !_FILTERS_H_

