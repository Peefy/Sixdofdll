#pragma once

#define AXIS_COUNT 6
// 所有单位 毫米
// Star Comming Para

#define DIS_PER_R 5.36
#define PULSE_PER_R 1280000
#define LENGTH_TO_PULSE_SCALE PULSE_PER_R / DIS_PER_R

#define SIN_30_DEG 0.5
#define COS_30_DEG 0.86602540378

#define SIN_60_DEG 0.86602540378
#define COS_60_DEG 0.5

#define RAD_30 (30.0 / 180.0 * pi)
#define RAD_60 (60.0 / 180.0 * pi)
