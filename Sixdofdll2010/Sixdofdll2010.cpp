// Sixdofdll2010.cpp : ���� DLL Ӧ�ó���ĵ���������
//

#include "stdafx.h"
#include "Sixdofdll2010.h"
#include "platform.h"
#include "MatrixOpreate.h"

Platform platform;

/*
����̬�ǻ�ȡ�׵��쳤�г̣�λ�Ƶ�λΪ����(mm)���Ƕȵ�λΪ�Ƕ�(deg)
@para
    x : X��λ�� mm
    y : y��λ�� mm
    z : z��λ�� mm
    roll : ����� deg
    yaw : ƫ���� deg
    pitch : ������ deg
@return
    double[0] ��1���쳤�� 
    double[1] ��2���쳤��
    double[2] ��3���쳤��
    double[3] ��4���쳤��
    double[4] ��5���쳤��
    double[5] ��6���쳤��
*/
SIXDOFDLL2010_API double* Control(double x, double y, double z, double roll, double yaw, double pitch)
{
    // deg 2 rad
    roll = roll * pi / 180.0;
    yaw = yaw * pi / 180.0;
    pitch = pitch * pi / 180.0;
    return platform.Control(x, y, z, roll, yaw, pitch);
}

/*
�ɸ׵��쳤�г̻�ȡ��̬�ǣ�λ�Ƶ�λΪ����(mm)���Ƕȵ�λΪ�Ƕ�(deg)
@para
lengths[0] ��1���쳤�� 
lengths[1] ��2���쳤��
lengths[2] ��3���쳤��
lengths[3] ��4���쳤��
lengths[4] ��5���쳤��
lengths[5] ��6���쳤��
@return
x : X��λ�� mm
y : y��λ�� mm
z : z��λ�� mm
roll : ����� deg
yaw : ƫ���� deg
pitch : ������ deg
*/
SIXDOFDLL2010_API double* FromLengthToPose(double * lengths)
{
	auto poses = platform.FromLengthToPose(lengths);
	// rad 2 deg
	poses[3] = poses[3] * 180.0 / pi;
	poses[4] = poses[4] * 180.0 / pi;
	poses[5] = poses[5] * 180.0 / pi;
	return poses;
}

/*
��ȡ ���������Ŀռ����꣬����ԭ��Ϊƽ̨����󶥲�ƽ̨���ģ���λΪ����(mm)
@return
double[0] xλ�� mm
double[1] yλ�� mm
double[2] zλ�� mm
*/
SIXDOFDLL2010_API double* GetTopPosition(int index)
{
    if (index < 0 || index >= AXIS_COUNT)
    {
        return nullptr;
    }
    return platform.HingePositions[index].ToArray();
}

/*
��ȡ �ײ������Ŀռ����꣬����ԭ��Ϊƽ̨����󶥲�ƽ̨���ģ���λΪ����(mm)
@para
index : �˵�������
@return
double[0] xλ�� mm
double[1] yλ�� mm
double[2] zλ�� mm
*/
SIXDOFDLL2010_API double* GetBottomPosition(int index)
{
    if (index < 0 || index >= AXIS_COUNT)
    {
        return nullptr;
    }
    return platform.BottomPositions[index].ToArray();
}

/*
���� ���������Ŀռ����꣬����ԭ��Ϊƽ̨����󶥲�ƽ̨���ģ���λΪ����(mm)
@para
index : �˵�������
x : x����
y : y����
z : z����
@return
bool �Ƿ����óɹ�
*/
SIXDOFDLL2010_API bool SetTopPosition(int index, double x, double y, double z)
{
    if (index < 0 || index >= AXIS_COUNT)
    {
        return false;
    }
    platform.HingePositions[index] = Point3D(x, y, z);
    return true;
}

/*
���� �ײ������Ŀռ����꣬����ԭ��Ϊƽ̨����󶥲�ƽ̨���ģ���λΪ����(mm)
@para
index : �˵�������
x : x����
y : y����
z : z����
@return
bool �Ƿ����óɹ�
*/
SIXDOFDLL2010_API bool SetBottomPosition(int index, double x, double y, double z)
{
    if (index < 0 || index >= AXIS_COUNT)
    {
        return false;
    }
    platform.BottomPositions[index] = Point3D(x, y, z);
    return true;
}

/*
���� �����ɶ�ƽ̨�Ľṹ��������λΪ����(mm)
@para
planeAboveHingeLength : ��ƽ̨�ϱ��������ƽ̨�����Ĵ�ֱ����
planeAboveBottomLength : ��ƽ̨�ϱ���������Ĵ�ֱ����
circleTopRadius : ��ƽ̨ԲȦ�뾶
circleBottomRadius : ��ƽ̨ԲȦ�뾶
distanceBetweenHinge : ��ƽ̨ͬһ���������������ľ���
@return
bool �Ƿ����óɹ�
*/
SIXDOFDLL2010_API void SetPlatformPara(double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
					 double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom)
{
	platform.SetPlatformPara(planeAboveHingeLength, planeAboveBottomLength, circleTopRadius, 
					 circleBottomRadius, distanceBetweenHingeTop, distanceBetweenHingeBottom);
}
