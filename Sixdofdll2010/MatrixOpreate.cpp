#include "stdafx.h"
#include "MatrixOpreate.h"

bool MatrixMultiplyVector(double matrix1[3][3], double vector[3], double(&result)[3]) {
    for (int k = 0; k < 3; k++) {
        result[k] = 0.0;
        for (int j = 0; j < 3; j++) {
            result[k] = result[k] + matrix1[k][j] * vector[j];
        }
    }
    return true;
}
