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
	cublasStatus_t blas_status = CUBLAS_STATUS_SUCCESS;
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
		&p->descr_d, dim, p->gpu_d, CUDA_R_32F);
	sparse_status = cusparseCreateDnVec(
		&p->descr_A_d, dim, p->gpu_A_d, CUDA_R_32F);
	sparse_status = cusparseCreateDnVec(
		&p->descr_A_alpha_d, dim, p->gpu_A_alpha_d, CUDA_R_32F);

	/* scratch buffer */
	float one = 1, zero = 0;
	size_t scratch_size;
	sparse_status = cusparseSpMV_bufferSize(
		p->sparse_handle,
		CUSPARSE_OPERATION_NON_TRANSPOSE,
		&one,
		p->descr_A,
		p->descr_d,
		&zero,
		p->descr_A_d,
		CUDA_R_32F,
		CUSPARSE_SPMV_ALG_DEFAULT,
		&scratch_size
	);
	cerr = cudaMalloc(&p->gpu_scratch, scratch_size);
}

void cuda_destroy(struct finite_element_problem *restrict p)
{
}

void gpu_conj_gradient(struct finite_element_problem *restrict p, float tolerance)
{
	const float zero = 0;
	const float one = 1;
	const float neg_one = -1;

	float bsquared;
	float alpha, beta, old_r2, dSd;
	int dim = p->b.dim;

	cublasSdot_v2(p->blas_handle, dim, p->gpu_b, 1, p->gpu_b, 1, &bsquared);
	cublasSscal_v2(p->blas_handle, dim, &zero, p->gpu_c, 1);
	cublasScopy_v2(p->blas_handle, dim, p->gpu_b, 1, p->gpu_r, 1);
	cublasScopy_v2(p->blas_handle, dim, p->gpu_b, 1, p->gpu_d, 1);

	for (int k = 1; k <= dim; k++) {
		cudaDeviceSynchronize();

		cublasSdot_v2(p->blas_handle, dim, p->gpu_r, 1, p->gpu_r, 1, &old_r2);
		if (bsquared > 0 && old_r2 / bsquared <= tolerance) {
			break;
		} else if (bsquared == 0 && old_r2 <= tolerance) {
			break;
		}

		/* Ad */
		cusparseSpMV(
			p->sparse_handle,
			CUSPARSE_OPERATION_NON_TRANSPOSE,
			&one,
			p->descr_A,
			p->descr_d,
			&zero,
			p->descr_A_d,
			CUDA_R_32F,
			CUSPARSE_SPMV_ALG_DEFAULT,
			p->gpu_scratch
		);
		cudaDeviceSynchronize();

		/* dAd */
		cublasSdot_v2(p->blas_handle, dim, p->gpu_d, 1, p->gpu_A_d, 1, &dSd);
		alpha = old_r2 / dSd;

		cublasSscal_v2(p->blas_handle, dim, &alpha, p->gpu_d, 1);
		cudaDeviceSynchronize();
		cublasSaxpy_v2(p->blas_handle, dim, &one, p->gpu_d, 1, p->gpu_c, 1);
		cudaDeviceSynchronize();

		/* A alpha d */
		cusparseSpMV(
			p->sparse_handle,
			CUSPARSE_OPERATION_NON_TRANSPOSE,
			&one,
			p->descr_A,
			p->descr_d,
			&zero,
			p->descr_A_alpha_d,
			CUDA_R_32F,
			CUSPARSE_SPMV_ALG_DEFAULT,
			p->gpu_scratch
		);
		cudaDeviceSynchronize();
		cublasSaxpy_v2(p->blas_handle, dim, &neg_one, p->gpu_A_alpha_d, 1, p->gpu_r, 1);
		cudaDeviceSynchronize();

		cublasSdot_v2(p->blas_handle, dim, p->gpu_r, 1, p->gpu_r, 1, &beta);
		beta /= old_r2 * alpha;
		cublasSscal_v2(p->blas_handle, dim, &beta, p->gpu_d, 1);
		cudaDeviceSynchronize();
		cublasSaxpy_v2(p->blas_handle, dim, &one, p->gpu_r, 1, p->gpu_d, 1);
	}

	cudaDeviceSynchronize();
	cudaMemcpy(p->c.x, p->gpu_c, dim*sizeof(*p->gpu_c), cudaMemcpyDeviceToHost);
}

#endif /* GPU_COMPUTE */
