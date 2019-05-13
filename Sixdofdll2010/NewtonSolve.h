
#ifndef _NEWTON_SOLVE_H_
#define _NEWTON_SOLVE_H_

// 六自由度平台正解算法：迭代法解非线性方程组
double* ForwardKinematics(double dlen1, double dlen2, double dlen3, double dlen4, double dlen5, double dlen6, 
						  double planeAboveHingeLength, double planeAboveBottomLength, double circleTopRadius, 
						  double circleBottomRadius, double distanceBetweenHingeTop, double distanceBetweenHingeBottom);

#endif
