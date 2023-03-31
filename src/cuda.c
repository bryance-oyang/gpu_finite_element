/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "finite_element.h"

#ifdef GPU_COMPUTE

#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>

int cuda_init(struct finite_element_problem *restrict p)
{
	cublasStatus_t blas_status;
	cusparseStatus_t sparse_status = CUSPARSE_STATUS_SUCCESS;
	cudaError_t cerr = cudaSuccess;

	int dim = p->b.dim;
	int nnz = p->A.len;

	sparse_status = cusparseCreate(&p->sparse_handle);
	blas_status = cublasCreate_v2(&p->blas_handle);

	cerr = cudaMalloc((void**)&p->gpu_rows, nnz*sizeof(*p->gpu_rows));
	cerr = cudaMalloc((void**)&p->gpu_cols, nnz*sizeof(*p->gpu_cols));
	cerr = cudaMalloc((void**)&p->gpu_A, nnz*sizeof(*p->gpu_A));

	cerr = cudaMalloc((void**)&p->gpu_b, dim*sizeof(*p->gpu_b));
	cerr = cudaMalloc((void**)&p->gpu_c, dim*sizeof(*p->gpu_c));

	cerr = cudaMalloc((void**)&p->gpu_r, dim*sizeof(*p->gpu_r));
	cerr = cudaMalloc((void**)&p->gpu_d, dim*sizeof(*p->gpu_d));
	cerr = cudaMalloc((void**)&p->gpu_tmp, dim*sizeof(*p->gpu_tmp));

	cerr = cudaMalloc((void**)&p->gpu_A_d, dim*sizeof(*p->gpu_A_d));
	cerr = cudaMalloc((void**)&p->gpu_A_alpha_d, dim*sizeof(*p->gpu_A_alpha_d));

	/* A */
	cerr = cudaMemcpy(p->gpu_rows, p->A.row, nnz*sizeof(p->gpu_rows), cudaMemcpyHostToDevice);
	cerr = cudaMemcpy(p->gpu_cols, p->A.col, nnz*sizeof(p->gpu_cols), cudaMemcpyHostToDevice);
	cerr = cudaMemcpy(p->gpu_A, p->A.A, nnz*sizeof(p->gpu_A), cudaMemcpyHostToDevice);

	/* b */
	cerr = cudaMemcpy(p->gpu_b, p->b.x, dim*sizeof(p->gpu_b), cudaMemcpyHostToDevice);

	/* sparse descr */
	sparse_status = cusparseCreateCoo(
		&p->descr_A, dim, dim, nnz, p->gpu_rows, p->gpu_cols, p->gpu_A,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
	sparse_status = cusparseCreateDnVec(
		&p->descr_A_d, dim, p->gpu_A_d, CUDA_R_32F);
	sparse_status = cusparseCreateDnVec(
		&p->descr_A_alpha_d, dim, p->gpu_A_alpha_d, CUDA_R_32F);
}

void cuda_destroy(struct finite_element_problem *restrict p)
{
}

#endif /* GPU_COMPUTE */
