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

	blas_status = cublasCreate_v2(&p->blas_handle);
	if (blas_status != CUBLAS_STATUS_SUCCESS) {
		goto err_blas;
	}
	sparse_status = cusparseCreate(&p->sparse_handle);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_sparse;
	}

	cerr = cudaMalloc((void**)&p->gpu_rows, nnz*sizeof(*p->gpu_rows));
	if (cerr != cudaSuccess) {
		goto err_rows;
	}
	cerr = cudaMalloc((void**)&p->gpu_cols, nnz*sizeof(*p->gpu_cols));
	if (cerr != cudaSuccess) {
		goto err_cols;
	}
	cerr = cudaMalloc((void**)&p->gpu_A, nnz*sizeof(*p->gpu_A));
	if (cerr != cudaSuccess) {
		goto err_A;
	}

	cerr = cudaMalloc((void**)&p->gpu_b, dim*sizeof(*p->gpu_b));
	if (cerr != cudaSuccess) {
		goto err_b;
	}
	cerr = cudaMalloc((void**)&p->gpu_c, dim*sizeof(*p->gpu_c));
	if (cerr != cudaSuccess) {
		goto err_c;
	}

	cerr = cudaMalloc((void**)&p->gpu_r, dim*sizeof(*p->gpu_r));
	if (cerr != cudaSuccess) {
		goto err_r;
	}
	cerr = cudaMalloc((void**)&p->gpu_d, dim*sizeof(*p->gpu_d));
	if (cerr != cudaSuccess) {
		goto err_d;
	}

	cerr = cudaMalloc((void**)&p->gpu_A_d, dim*sizeof(*p->gpu_A_d));
	if (cerr != cudaSuccess) {
		goto err_A_d;
	}

	/* A */
	cerr = cudaMemcpy(p->gpu_rows, p->A.row, nnz*sizeof(*p->gpu_rows), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}
	cerr = cudaMemcpy(p->gpu_cols, p->A.col, nnz*sizeof(*p->gpu_cols), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}
	cerr = cudaMemcpy(p->gpu_A, p->A.A, nnz*sizeof(*p->gpu_A), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}

	/* b */
	cerr = cudaMemcpy(p->gpu_b, p->b.x, dim*sizeof(*p->gpu_b), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}

	/* sparse descr */
	sparse_status = cusparseCreateCoo(
		&p->descr_A, dim, dim, nnz, p->gpu_rows, p->gpu_cols, p->gpu_A,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_descr_A;
	}
	sparse_status = cusparseCreateDnVec(
		&p->descr_d, dim, p->gpu_d, CUDA_R_32F);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_descr_d;
	}
	sparse_status = cusparseCreateDnVec(
		&p->descr_A_d, dim, p->gpu_A_d, CUDA_R_32F);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_descr_A_d;
	}

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
		CUSPARSE_SPMV_COO_ALG2,
		&scratch_size
	);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_scratch_size;
	}
	cerr = cudaMalloc(&p->gpu_scratch, scratch_size);
	if (cerr != cudaSuccess) {
		goto err_scratch;
	}

	return 0;

err_scratch:
err_scratch_size:
	cusparseDestroyDnVec(p->descr_A_d);
err_descr_A_d:
	cusparseDestroyDnVec(p->descr_d);
err_descr_d:
	cusparseDestroySpMat(p->descr_A);
err_descr_A:
err_memcpy:
	cudaFree(p->gpu_A_d);
err_A_d:
	cudaFree(p->gpu_d);
err_d:
	cudaFree(p->gpu_r);
err_r:
	cudaFree(p->gpu_c);
err_c:
	cudaFree(p->gpu_b);
err_b:
	cudaFree(p->gpu_A);
err_A:
	cudaFree(p->gpu_cols);
err_cols:
	cudaFree(p->gpu_rows);
err_rows:
	cusparseDestroy(p->sparse_handle);
err_sparse:
	cublasDestroy_v2(p->blas_handle);
err_blas:
	return -1;
}

void cuda_destroy(struct finite_element_problem *restrict p)
{
	cudaFree(p->gpu_scratch);
	cusparseDestroyDnVec(p->descr_A_d);
	cusparseDestroyDnVec(p->descr_d);
	cusparseDestroySpMat(p->descr_A);
	cudaFree(p->gpu_A_d);
	cudaFree(p->gpu_d);
	cudaFree(p->gpu_r);
	cudaFree(p->gpu_c);
	cudaFree(p->gpu_b);
	cudaFree(p->gpu_A);
	cudaFree(p->gpu_cols);
	cudaFree(p->gpu_rows);
	cusparseDestroy(p->sparse_handle);
	cublasDestroy_v2(p->blas_handle);
}

int gpu_conj_gradient(struct finite_element_problem *restrict p, float tolerance)
{
	cublasStatus_t blas_status = CUBLAS_STATUS_SUCCESS;
	cusparseStatus_t sparse_status = CUSPARSE_STATUS_SUCCESS;
	cudaError_t cerr = cudaSuccess;

	const float zero = 0;
	const float one = 1;
	const float neg_one = -1;

	float bsquared;
	float alpha, beta, old_r2, dAd;
	int dim = p->b.dim;

	blas_status |= cublasSdot_v2(p->blas_handle, dim, p->gpu_b, 1, p->gpu_b, 1, &bsquared);
	blas_status |= cublasSscal_v2(p->blas_handle, dim, &zero, p->gpu_c, 1);
	blas_status |= cublasScopy_v2(p->blas_handle, dim, p->gpu_b, 1, p->gpu_r, 1);
	blas_status |= cublasScopy_v2(p->blas_handle, dim, p->gpu_b, 1, p->gpu_d, 1);

	for (int k = 1; k <= dim; k++) {
		blas_status |= cublasSdot_v2(p->blas_handle, dim, p->gpu_r, 1, p->gpu_r, 1, &old_r2);
		if (bsquared > 0 && old_r2 / bsquared <= tolerance) {
			break;
		} else if (bsquared == 0 && old_r2 <= tolerance) {
			break;
		}

		/* Ad */
		sparse_status |= cusparseSpMV(
			p->sparse_handle,
			CUSPARSE_OPERATION_NON_TRANSPOSE,
			&one,
			p->descr_A,
			p->descr_d,
			&zero,
			p->descr_A_d,
			CUDA_R_32F,
			CUSPARSE_SPMV_COO_ALG2,
			p->gpu_scratch
		);

		/* dAd */
		blas_status |= cublasSdot_v2(p->blas_handle, dim, p->gpu_d, 1, p->gpu_A_d, 1, &dAd);
		alpha = old_r2 / dAd;

		blas_status |= cublasSscal_v2(p->blas_handle, dim, &alpha, p->gpu_A_d, 1);
		blas_status |= cublasSscal_v2(p->blas_handle, dim, &alpha, p->gpu_d, 1);
		blas_status |= cublasSaxpy_v2(p->blas_handle, dim, &one, p->gpu_d, 1, p->gpu_c, 1);

		blas_status |= cublasSaxpy_v2(p->blas_handle, dim, &neg_one, p->gpu_A_d, 1, p->gpu_r, 1);

		blas_status |= cublasSdot_v2(p->blas_handle, dim, p->gpu_r, 1, p->gpu_r, 1, &beta);
		beta /= old_r2 * alpha;
		blas_status |= cublasSscal_v2(p->blas_handle, dim, &beta, p->gpu_d, 1);
		blas_status |= cublasSaxpy_v2(p->blas_handle, dim, &one, p->gpu_r, 1, p->gpu_d, 1);
	}

	cerr = cudaMemcpy(p->c.x, p->gpu_c, dim*sizeof(*p->gpu_c), cudaMemcpyDeviceToHost);

	if (blas_status || sparse_status || cerr) {
		return -1;
	}

	return 0;
}

#endif /* GPU_COMPUTE */
