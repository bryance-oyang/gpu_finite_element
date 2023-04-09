/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 */

#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include "mesh.h"

struct vis;

#ifdef GPU_COMPUTE
#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>
#endif /* GPU_COMPUTE */

/* Ac = b */
struct finite_element_problem {
	struct mesh *mesh;
	/** stiffness */
	struct sparse A;
	/** forces */
	struct vec b;
	/** solution */
	struct vec c;

#ifdef GPU_COMPUTE

	cublasHandle_t blas_handle;
	cusparseHandle_t sparse_handle;

	int32_t *gpu_rows;
	int32_t *gpu_cols;
	double *gpu_A;

	double *gpu_b;
	double *gpu_c;

	double *gpu_r;
	double *gpu_d;

	double *gpu_A_d;

	cusparseSpMatDescr_t descr_A;
	cusparseDnVecDescr_t descr_d;
	cusparseDnVecDescr_t descr_A_d;

	void *gpu_scratch;

#endif /* GPU_COMPUTE */
};

int fep_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh);
void fep_destroy(struct finite_element_problem *restrict p);
int fep_solve(struct finite_element_problem *restrict p, double tolerance, struct vis *restrict vis);

#endif /* FINITE_ELEMENT_H */
