// Sixdofdll2010.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "Sixdofdll2010.h"
#include "platform.h"
#include "MatrixOpreate.h"

Platform platform;

/*
由姿态角获取缸的伸长行程，位移单位为毫米(mm)，角度单位为角度(deg)
@para
    x : X轴位移 mm
    y : y轴位移 mm
    z : z轴位移 mm
    roll : 横滚角 deg
    yaw : 偏航角 deg
    pitch : 俯仰角 deg
@return
    double[0] 杆1的伸长量 
    double[1] 杆2的伸长量
    double[2] 杆3的伸长量
    double[3] 杆4的伸长量
    double[4] 杆5的伸长量
    double[5] 杆6的伸长量
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
由缸的伸长行程获取姿态角，位移单位为毫米(mm)，角度单位为角度(deg)
@para
lengths[0] 杆1的伸长量 
lengths[1] 杆2的伸长量
lengths[2] 杆3的伸长量
lengths[3] 杆4的伸长量
lengths[4] 杆5的伸长量
lengths[5] 杆6的伸长量
@return
x : X轴位移 mm
y : y轴位移 mm
z : z轴位移 mm
roll : 横滚角 deg
yaw : 偏航角 deg
pitch : 俯仰角 deg
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
获取 顶部铰链的空间坐标，坐标原点为平台升起后顶部平台中心，单位为毫米(mm)
@return
double[0] x位移 mm
double[1] y位移 mm
double[2] z位移 mm
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
获取 底部铰链的空间坐标，坐标原点为平台升起后顶部平台中心，单位为毫米(mm)
@para
index : 杆的索引号
@return
double[0] x位移 mm
double[1] y位移 mm
double[2] z位移 mm
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
设置 顶部铰链的空间坐标，坐标原点为平台升起后顶部平台中心，单位为毫米(mm)
@para
index : 杆的索引号
x : x坐标
y : y坐标
z : z坐标
@return
bool 是否设置成功
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
设置 底部铰链的空间坐标，坐标原点为平台升起后顶部平台中心，单位为毫米(mm)
@para
index : 杆的索引号
x : x坐标
y : y坐标
z : z坐标
@return
bool 是否设置成功
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
设置 六自由度平台的结构参数，单位为毫米(mm)
@para
planeAboveHingeLength : 上平台上表面距离上平台铰链的垂直距离
planeAboveBottomLength : 上平台上表面距离地面的垂直距离
circleTopRadius : 上平台圆圈半径
circleBottomRadius : 下平台圆圈半径
distanceBetweenHinge : 上平台同一组两个铰链的中心距离
@return
bool 是否设置成功
*/
SIXDOFDLL2010_API void SetPlatformPara(double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
					 double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom)
{
	platform.SetPlatformPara(planeAboveHingeLength, planeAboveBottomLength, circleTopRadius, 
					 circleBottomRadius, distanceBetweenHingeTop, distanceBetweenHingeBottom);
}
