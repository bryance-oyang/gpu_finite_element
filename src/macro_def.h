#ifndef MACRO_DEF_H
#define MACRO_DEF_H

#include <float.h>

#ifdef DOUBLE_PRECISION
typedef double number;
#define NUMBER_MAX DBL_MAX
#define EPSILON 1e-12
#else
typedef float number;
#define NUMBER_MAX FLT_MAX
#define EPSILON 1e-12
#endif

#define DIM 2

#define SQR(x) ((x)*(x))

#define IMAGE_WIDTH 1024
#define IMAGE_HEIGHT 1024

#endif /* MACRO_DEF_H */
