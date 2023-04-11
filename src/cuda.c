/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @brief GPU accelerated solve
 */

#include "solver.h"

#ifdef GPU_COMPUTE

#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>

int cuda_init(struct solver *restrict solver)
{
	cublasStatus_t blas_status = CUBLAS_STATUS_SUCCESS;
	cusparseStatus_t sparse_status = CUSPARSE_STATUS_SUCCESS;
	cudaError_t cerr = cudaSuccess;

	int dim = solver->b.dim;
	int nnz = solver->A.len;

	blas_status = cublasCreate_v2(&solver->blas_handle);
	if (blas_status != CUBLAS_STATUS_SUCCESS) {
		goto err_blas;
	}
	sparse_status = cusparseCreate(&solver->sparse_handle);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_sparse;
	}

	cerr = cudaMalloc((void**)&solver->gpu_rows, nnz*sizeof(*solver->gpu_rows));
	if (cerr != cudaSuccess) {
		goto err_rows;
	}
	cerr = cudaMalloc((void**)&solver->gpu_cols, nnz*sizeof(*solver->gpu_cols));
	if (cerr != cudaSuccess) {
		goto err_cols;
	}
	cerr = cudaMalloc((void**)&solver->gpu_A, nnz*sizeof(*solver->gpu_A));
	if (cerr != cudaSuccess) {
		goto err_A;
	}

	cerr = cudaMalloc((void**)&solver->gpu_b, dim*sizeof(*solver->gpu_b));
	if (cerr != cudaSuccess) {
		goto err_b;
	}
	cerr = cudaMalloc((void**)&solver->gpu_c, dim*sizeof(*solver->gpu_c));
	if (cerr != cudaSuccess) {
		goto err_c;
	}

	cerr = cudaMalloc((void**)&solver->gpu_r, dim*sizeof(*solver->gpu_r));
	if (cerr != cudaSuccess) {
		goto err_r;
	}
	cerr = cudaMalloc((void**)&solver->gpu_d, dim*sizeof(*solver->gpu_d));
	if (cerr != cudaSuccess) {
		goto err_d;
	}

	cerr = cudaMalloc((void**)&solver->gpu_A_d, dim*sizeof(*solver->gpu_A_d));
	if (cerr != cudaSuccess) {
		goto err_A_d;
	}

	/* A */
	cerr = cudaMemcpy(solver->gpu_rows, solver->A.row, nnz*sizeof(*solver->gpu_rows), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}
	cerr = cudaMemcpy(solver->gpu_cols, solver->A.col, nnz*sizeof(*solver->gpu_cols), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}
	cerr = cudaMemcpy(solver->gpu_A, solver->A.A, nnz*sizeof(*solver->gpu_A), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}

	/* b */
	cerr = cudaMemcpy(solver->gpu_b, solver->b.x, dim*sizeof(*solver->gpu_b), cudaMemcpyHostToDevice);
	if (cerr != cudaSuccess) {
		goto err_memcpy;
	}

	/* sparse descr */
	sparse_status = cusparseCreateCoo(
		&solver->descr_A, dim, dim, nnz, solver->gpu_rows, solver->gpu_cols, solver->gpu_A,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_descr_A;
	}
	sparse_status = cusparseCreateDnVec(
		&solver->descr_d, dim, solver->gpu_d, CUDA_R_64F);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_descr_d;
	}
	sparse_status = cusparseCreateDnVec(
		&solver->descr_A_d, dim, solver->gpu_A_d, CUDA_R_64F);
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_descr_A_d;
	}

	/* scratch buffer */
	double one = 1, zero = 0;
	size_t scratch_size;
	sparse_status = cusparseSpMV_bufferSize(
		solver->sparse_handle,
		CUSPARSE_OPERATION_NON_TRANSPOSE,
		&one,
		solver->descr_A,
		solver->descr_d,
		&zero,
		solver->descr_A_d,
		CUDA_R_64F,
		CUSPARSE_SPMV_COO_ALG2,
		&scratch_size
	);
	scratch_size++; /* CUDA doesn't want NULL: prevent NULL ptr */
	if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
		goto err_scratch_size;
	}
	cerr = cudaMalloc(&solver->gpu_scratch, scratch_size);
	if (cerr != cudaSuccess) {
		goto err_scratch;
	}

	return 0;

err_scratch:
err_scratch_size:
	cusparseDestroyDnVec(solver->descr_A_d);
err_descr_A_d:
	cusparseDestroyDnVec(solver->descr_d);
err_descr_d:
	cusparseDestroySpMat(solver->descr_A);
err_descr_A:
err_memcpy:
	cudaFree(solver->gpu_A_d);
err_A_d:
	cudaFree(solver->gpu_d);
err_d:
	cudaFree(solver->gpu_r);
err_r:
	cudaFree(solver->gpu_c);
err_c:
	cudaFree(solver->gpu_b);
err_b:
	cudaFree(solver->gpu_A);
err_A:
	cudaFree(solver->gpu_cols);
err_cols:
	cudaFree(solver->gpu_rows);
err_rows:
	cusparseDestroy(solver->sparse_handle);
err_sparse:
	cublasDestroy_v2(solver->blas_handle);
err_blas:
	return -1;
}

void cuda_destroy(struct solver *restrict solver)
{
	cudaFree(solver->gpu_scratch);
	cusparseDestroyDnVec(solver->descr_A_d);
	cusparseDestroyDnVec(solver->descr_d);
	cusparseDestroySpMat(solver->descr_A);
	cudaFree(solver->gpu_A_d);
	cudaFree(solver->gpu_d);
	cudaFree(solver->gpu_r);
	cudaFree(solver->gpu_c);
	cudaFree(solver->gpu_b);
	cudaFree(solver->gpu_A);
	cudaFree(solver->gpu_cols);
	cudaFree(solver->gpu_rows);
	cusparseDestroy(solver->sparse_handle);
	cublasDestroy_v2(solver->blas_handle);
}

int gpu_conj_gradient(struct solver *restrict solver, double tolerance)
{
	cublasStatus_t blas_status = CUBLAS_STATUS_SUCCESS;
	cusparseStatus_t sparse_status = CUSPARSE_STATUS_SUCCESS;
	cudaError_t cerr = cudaSuccess;

	const double zero = 0;
	const double one = 1;
	const double neg_one = -1;

	double bsquared;
	double alpha, beta, old_r2, dAd;
	int dim = solver->b.dim;

	blas_status |= cublasDdot_v2(solver->blas_handle, dim, solver->gpu_b, 1, solver->gpu_b, 1, &bsquared);
	blas_status |= cublasDscal_v2(solver->blas_handle, dim, &zero, solver->gpu_c, 1);
	blas_status |= cublasDcopy_v2(solver->blas_handle, dim, solver->gpu_b, 1, solver->gpu_r, 1);
	blas_status |= cublasDcopy_v2(solver->blas_handle, dim, solver->gpu_b, 1, solver->gpu_d, 1);

	for (int k = 1; k <= dim; k++) {
		blas_status |= cublasDdot_v2(solver->blas_handle, dim, solver->gpu_r, 1, solver->gpu_r, 1, &old_r2);
		if (bsquared > 0 && old_r2 / bsquared <= tolerance) {
			break;
		} else if (bsquared == 0 && old_r2 <= tolerance) {
			break;
		}

		/* Ad */
		sparse_status |= cusparseSpMV(
			solver->sparse_handle,
			CUSPARSE_OPERATION_NON_TRANSPOSE,
			&one,
			solver->descr_A,
			solver->descr_d,
			&zero,
			solver->descr_A_d,
			CUDA_R_64F,
			CUSPARSE_SPMV_COO_ALG2,
			solver->gpu_scratch
		);

		/* dAd */
		blas_status |= cublasDdot_v2(solver->blas_handle, dim, solver->gpu_d, 1, solver->gpu_A_d, 1, &dAd);

		/* Ad = alpha Ad; d = alpha d; */
		alpha = old_r2 / dAd;
		blas_status |= cublasDscal_v2(solver->blas_handle, dim, &alpha, solver->gpu_A_d, 1);
		blas_status |= cublasDscal_v2(solver->blas_handle, dim, &alpha, solver->gpu_d, 1);

		/* c += alpha d */
		blas_status |= cublasDaxpy_v2(solver->blas_handle, dim, &one, solver->gpu_d, 1, solver->gpu_c, 1);

		/* r -= alpha Ad */
		blas_status |= cublasDaxpy_v2(solver->blas_handle, dim, &neg_one, solver->gpu_A_d, 1, solver->gpu_r, 1);

		blas_status |= cublasDdot_v2(solver->blas_handle, dim, solver->gpu_r, 1, solver->gpu_r, 1, &beta);
		beta /= old_r2 * alpha;
		blas_status |= cublasDscal_v2(solver->blas_handle, dim, &beta, solver->gpu_d, 1);
		blas_status |= cublasDaxpy_v2(solver->blas_handle, dim, &one, solver->gpu_r, 1, solver->gpu_d, 1);
	}

	cerr = cudaMemcpy(solver->c.x, solver->gpu_c, dim*sizeof(*solver->gpu_c), cudaMemcpyDeviceToHost);

	if (blas_status || sparse_status || cerr) {
		return -1;
	}

	return 0;
}

#endif /* GPU_COMPUTE */
