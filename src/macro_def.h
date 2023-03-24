#ifndef MACRO_DEF_H
#define MACRO_DEF_H

#include <float.h>

// watch conjugate gradient solve the system (CPU only, no gpu)
#ifndef GPU_COMPUTE
#define ANIMATE 1
#endif

typedef float number;
#define NUMBER_MAX FLT_MAX
#define EPSILON 1e-12

#define DIM 2

#define SQR(x) ((x)*(x))

#define IMAGE_WIDTH 1024
#define IMAGE_HEIGHT 1024

#endif /* MACRO_DEF_H */
