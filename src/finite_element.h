/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include "mesh.h"
#include "visualize.h"

#ifdef GPU_COMPUTE
#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>
#endif /* GPU_COMPUTE */

/* Ac = b */
struct finite_element_problem {
	struct mesh *mesh;
	struct sparse A;
	struct vec b;
	struct vec c;

#ifdef GPU_COMPUTE

	cublasHandle_t blas_handle;
	cusparseHandle_t sparse_handle;

	float *gpu_rows;
	float *gpu_cols;
	float *gpu_A;

	float *gpu_b;
	float *gpu_c;

	float *gpu_r;
	float *gpu_d;
	float *gpu_tmp;

	float *gpu_A_d;
	float *gpu_A_alpha_d;

	cusparseSpMatDescr_t descr_A;
	cusparseDnVecDescr_t descr_A_d;
	cusparseDnVecDescr_t descr_A_alpha_d;

#endif /* GPU_COMPUTE */
};

int fep_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh);
void fep_destroy(struct finite_element_problem *restrict p);
void fep_solve(struct finite_element_problem *restrict p, float tolerance, struct vis *vis);

void fep_scalar_stress(struct finite_element_problem *restrict p);

#endif /* FINITE_ELEMENT_H */
