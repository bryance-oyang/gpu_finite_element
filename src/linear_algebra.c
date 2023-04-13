/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @brief non-GPU accelerated and unoptimized versions
 */

#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linear_algebra.h"
#include "mesh.h"
#include "visualize.h"

/*
 * comparison function to keep sparse matrix sorted:
 * (row, col) < (row, col)?
 */
static int sparse_row_col_cmp(int r0, int c0, int r1, int c1)
{
	if (r0 != r1) {
		return (r0 > r1) - (r0 < r1);
	} else {
		return (c0 > c1) - (c0 < c1);
	}
}

int sparse_init(struct sparse *restrict S)
{
	int size = 8;
	S->len = 0;
	S->size = size;

	S->row = malloc(size * sizeof(*S->row));
	if (S->row == NULL) {
		goto err_norow;
	}

	S->col = malloc(size * sizeof(*S->col));
	if (S->col == NULL) {
		goto err_nocol;
	}

	S->A = malloc(size * sizeof(*S->A));
	if (S->A == NULL) {
		goto err_noentry;
	}
	return 0;

err_noentry:
	free(S->col);
err_nocol:
	free(S->row);
err_norow:
	return -1;
}

void sparse_destroy(struct sparse *restrict S)
{
	S->len = 0;
	S->size = 0;
	free(S->row);
	free(S->col);
	free(S->A);
}

/* insert nonexistent */
int sparse_add(struct sparse *restrict S, int row, int col, double entry)
{
	if (S->len == S->size) {
		int new_size = 2 * S->size;
		void *tmp_row, *tmp_col, *tmp_A;

		if ((tmp_row = realloc(S->row, new_size * sizeof(*S->row))) == NULL) {
			return -1;
		}
		S->row = tmp_row;

		if ((tmp_col = realloc(S->col, new_size * sizeof(*S->col))) == NULL) {
			return -1;
		}
		S->col = tmp_col;

		if ((tmp_A = realloc(S->A, new_size * sizeof(*S->A))) == NULL) {
			return -1;
		}
		S->A = tmp_A;

		S->size = new_size;
	}

	S->row[S->len] = row;
	S->col[S->len] = col;
	S->A[S->len] = entry;
	S->len++;

	return 0;
}

static inline void swap_int(int *a, int *b)
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

static inline void swap_double(double *a, double *b)
{
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

/* partition function for quicksort */
static int sparse_partition(struct sparse *restrict S, int lo, int hi)
{
	int i, j;
	int prow = S->row[hi - 1];
	int pcol = S->col[hi - 1];

	for (i = lo, j = lo; j < hi - 1; j++) {
		if (sparse_row_col_cmp(S->row[j], S->col[j], prow, pcol) < 0) {
			swap_int(&S->row[i], &S->row[j]);
			swap_int(&S->col[i], &S->col[j]);
			swap_double(&S->A[i], &S->A[j]);
			i++;
		}
	}

	swap_int(&S->row[i], &S->row[hi - 1]);
	swap_int(&S->col[i], &S->col[hi - 1]);
	swap_double(&S->A[i], &S->A[hi - 1]);
	return i;
}

/* sort into row-major format */
static void sparse_qsort(struct sparse *restrict S, int lo, int hi)
{
	if (hi - lo < 2) {
		return;
	}

	int mid = sparse_partition(S, lo, hi);
	sparse_qsort(S, lo, mid);
	sparse_qsort(S, mid + 1, hi);
}

/* sort into row-major format */
void sparse_sort(struct sparse *restrict S)
{
	sparse_qsort(S, 0, S->len);
}

/* S should be sorted: duplicate entries for row, col are summed into one */
void sparse_consolidate(struct sparse *restrict S)
{
	int i, j;
	int row, col;

	if (S->len < 2) {
		return;
	}

	row = S->row[0];
	col = S->col[0];
	for (i = 0, j = 1; j < S->len; j++) {
		if (row == S->row[j] && col == S->col[j]) {
			S->A[i] += S->A[j];
		} else {
			i++;
			row = S->row[j];
			col = S->col[j];
			S->row[i] = row;
			S->col[i] = col;
			S->A[i] = S->A[j];
		}
	}

	S->len = i + 1;
}

static void lu_swap_row(double *restrict lu, int dim, int *restrict row_idx,
	int row0, int row1)
{
	/* swap rows */
	for (int j = 0; j < dim; j++) {
		double entry = lu[row0*dim + j];
		lu[row0*dim + j] = lu[row1*dim + j];
		lu[row1*dim + j] = entry;
	}

	/* swap row indices */
	int tmp;
	tmp = row_idx[row0];
	row_idx[row0] = row_idx[row1];
	row_idx[row1] = tmp;
}

/**
 * LU Decomposition
 *
 * matrix will be replaced by its decomposition
 *
 * row permute of matrix = A B where A is lower triangular with 1's on diag and
 * B is upper and combine both A and B like [A\B] in original matrix
 *
 * row_idx should be an input of row indices [0, 1, 2, ...] and will be permuted
 * to indicate the row permutation applied to matrix
 *
 * row permutations are necessary for numerical stability (div by small number)
 * and prevent div by 0
 *
 * returns (-1)^(# of swaps)
 */
int lu_decomp(double *restrict matrix, int dim, int *restrict row_idx)
{
	int parity = 1;

	for (int j = 0; j < dim; j++) {
		/* determine b_ij for i <= j */
		for (int i = 0; i <= j; i++) {
			double ab = 0;
			for (int k = 0; k < i; k++) {
				ab += matrix[i*dim + k] * matrix[k*dim + j];
			}
			matrix[i*dim + j] = matrix[i*dim + j] - ab;
		}

		/* determine a_ij without dividing by b_jj for i > j */
		double max_abs_b = fabs(matrix[j*dim + j]);
		int max_row = j;
		for (int i = j + 1; i < dim; i++) {
			double ab = 0;
			for (int k = 0; k < j; k++) {
				ab += matrix[i*dim + k] * matrix[k*dim + j];
			}
			double b = matrix[i*dim + j] - ab;
			double abs_b = fabs(b);
			matrix[i*dim + j] = b;
			if (abs_b > max_abs_b) {
				max_abs_b = abs_b;
				max_row = i;
			}
		}

		if (max_row != j) {
			parity *= -1;
			lu_swap_row(matrix, dim, row_idx, max_row, j);
		}

		/* divide a_ij by b_jj */
		for (int i = j + 1; i < dim; i++) {
			matrix[i*dim + j] /= matrix[j*dim + j];
		}
	}

	return parity;
}

/* perform LU Decomposition and do 2 linear solves */
void inverse_matrix(double *restrict matrix, int dim, double *restrict inverse)
{
	double *restrict tmp = malloc(dim * dim * sizeof(*tmp));
	int *restrict row_idx = malloc(dim * sizeof(*row_idx));
	if (tmp == NULL || row_idx == NULL) {
		raise(SIGSEGV);
	}

	for (int i = 0; i < dim; i++) {
		row_idx[i] = i;
	}

	lu_decomp(matrix, dim, row_idx);

	/* solve Ly = I */
	for (int j = 0; j < dim; j++) {
		for (int i = 0; i < j; i++) {
			tmp[i*dim + j] = 0;
		}

		tmp[j*dim + j] = 1.0;

		for (int i = j + 1; i < dim; i++) {
			double ay = 0;
			for (int k = j; k < i; k++) {
				ay += matrix[i*dim + k] * tmp[k*dim + j];
			}
			tmp[i*dim + j] = -ay;
		}
	}

	/* solve Ux = y */
	for (int j = 0; j < dim; j++) {
		tmp[(dim - 1)*dim + j] = tmp[(dim - 1)*dim + j] / matrix[(dim - 1)*dim + (dim - 1)];
		for (int i = dim - 2; i >= 0; i--) {
			double bx = 0;
			for (int k = i + 1; k < dim; k++) {
				bx += matrix[i*dim + k] * tmp[k*dim + j];
			}
			tmp[i*dim + j] = (tmp[i*dim + j] - bx) / matrix[i*dim + i];
		}
	}

	/* permute columns back according to the way rows were permuted in LU */
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			inverse[i*dim + row_idx[j]] = tmp[i*dim + j];
		}
	}

	free(tmp);
	free(row_idx);
}

/* 2x2 determinant */
double matrix_det2(double *restrict matrix)
{
	return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

/* invert 2x2 */
void inverse_matrix2(double *restrict matrix, double *restrict inverse)
{
	double inv_det = 1.0 / matrix_det2(matrix);
	inverse[0] = matrix[3] * inv_det;
	inverse[3] = matrix[0] * inv_det;
	inverse[1] = -matrix[1] * inv_det;
	inverse[2] = -matrix[2] * inv_det;
}

int vec_init(struct vec *restrict v, int dim)
{
	v->dim = dim;
	v->x = malloc(dim * sizeof(*v->x));
	if (v->x == NULL) {
		return -1;
	}
	for (int i = 0; i < dim; i++) {
		v->x[i] = 0;
	}
	return 0;
}

void vec_destroy(struct vec *restrict v)
{
	v->dim = 0;
	free(v->x);
}

/* out = Sv */
void sparse_mult_vec(const struct sparse *restrict S, const struct vec *restrict v, struct vec *restrict out)
{
	for (int i = 0; i < out->dim; i++) {
		out->x[i] = 0;
	}

	for (int n = 0; n < S->len; n++) {
		int i = S->row[n];
		int k = S->col[n];
		out->x[i] += S->A[n] * v->x[k];
	}
}

double vec_dot(const struct vec *a, const struct vec *b)
{
	double result = 0;
	int dim = a->dim;

	for (int i = 0; i < dim; i++) {
		result += a->x[i] * b->x[i];
	}
	return result;
}

/* they must have the same dimensions */
void vec_copy(const struct vec *restrict in, struct vec *restrict out)
{
	for (int i = 0; i < in->dim; i++) {
		out->x[i] = in->x[i];
	}
}

/* they must have the same dimensions */
void vec_add(struct vec *a, struct vec *b, struct vec *out)
{
	int dim = out->dim;
	for (int i = 0; i < dim; i++) {
		out->x[i] = a->x[i] + b->x[i];
	}
}

/* they must have the same dimensions */
void vec_sub(struct vec *a, struct vec *b, struct vec *out)
{
	int dim = out->dim;
	for (int i = 0; i < dim; i++) {
		out->x[i] = a->x[i] - b->x[i];
	}
}

/* scalar multiply */
void vec_scale(double scalar, struct vec *restrict v)
{
	for (int i = 0; i < v->dim; i++) {
		v->x[i] *= scalar;
	}
}

void vec2_copy(struct vec2 *restrict in, struct vec2 *restrict out)
{
	for (int i = 0; i < 2; i++) {
		out->x[i] = in->x[i];
	}
}

void vec2_add(struct vec2 *a, struct vec2 *b, struct vec2 *out)
{
	for (int i = 0; i < 2; i++) {
		out->x[i] = a->x[i] + b->x[i];
	}
}

void vec2_sub(struct vec2 *a, struct vec2 *b, struct vec2 *out)
{
	for (int i = 0; i < 2; i++) {
		out->x[i] = a->x[i] - b->x[i];
	}
}

void vec2_scale(double scalar, struct vec2 *restrict v)
{
	for (int i = 0; i < 2; i++) {
		v->x[i] *= scalar;
	}
}

double vec2_dot(struct vec2 *a, struct vec2 *b)
{
	double result = 0;
	for (int i = 0; i < 2; i++) {
		result += a->x[i] * b->x[i];
	}
	return result;
}

void vec2_normalize(struct vec2 *restrict v)
{
	double norm = sqrt(vec2_dot(v, v));
	vec2_scale(1.0/norm, v);
}

void vec2_midpoint(struct vec2 *a, struct vec2 *b, struct vec2 *out)
{
	vec2_add(a, b, out);
	vec2_scale(0.5, out);
}

int sparse_conj_grad(struct sparse *restrict A, struct vec *restrict b,
	struct vec *restrict c, double tolerance, struct vis *restrict vis,
	struct mesh *restrict mesh)
{
	double bsquared;
	double alpha, beta, old_r2, dAd;
	struct vec r, d, A_d;

	bsquared = vec_dot(b, b);

	vec_init(&r, c->dim);
	vec_init(&d, c->dim);
	vec_init(&A_d, c->dim);

	vec_scale(0, c);
	vec_copy(b, &r);
	vec_copy(b, &d);

	for (int k = 1; k <= c->dim; k++) {
#ifdef ANIMATE
		vis_fill(vis, mesh, c);
		vis_send(vis);
		usleep(15000);
#else /* ANIMATE */
		(void)vis;
		(void)mesh;
#endif /* ANIMATE */

		old_r2 = vec_dot(&r, &r);
		if (bsquared > 0 && old_r2 / bsquared <= tolerance) {
			break;
		} else if (bsquared == 0 && old_r2 <= tolerance) {
			break;
		}

		/* dAd */
		sparse_mult_vec(A, &d, &A_d);
		dAd = vec_dot(&d, &A_d);

		/* Ad = alpha Ad; d = alpha d; */
		alpha = old_r2 / dAd;
		vec_scale(alpha, &A_d);
		vec_scale(alpha, &d);

		/* c += alpha d */
		vec_add(&d, c, c);

		/* r -= alpha Ad */
		vec_sub(&r, &A_d, &r);

		beta = vec_dot(&r, &r) / old_r2;
		vec_scale(beta / alpha, &d);
		vec_add(&r, &d, &d);
	}

	vec_destroy(&r);
	vec_destroy(&d);
	vec_destroy(&A_d);

	return 0;
}
