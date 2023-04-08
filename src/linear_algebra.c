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
int sparse_add(struct sparse *restrict S, int row, int col, float entry)
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

static inline void swap_float(float *a, float *b)
{
	float tmp = *a;
	*a = *b;
	*b = tmp;
}

static int sparse_partition(struct sparse *restrict S, int lo, int hi)
{
	int i, j;
	int prow = S->row[hi - 1];
	int pcol = S->col[hi - 1];

	for (i = lo, j = lo; j < hi - 1; j++) {
		if (sparse_row_col_cmp(S->row[j], S->col[j], prow, pcol) < 0) {
			swap_int(&S->row[i], &S->row[j]);
			swap_int(&S->col[i], &S->col[j]);
			swap_float(&S->A[i], &S->A[j]);
			i++;
		}
	}

	swap_int(&S->row[i], &S->row[hi - 1]);
	swap_int(&S->col[i], &S->col[hi - 1]);
	swap_float(&S->A[i], &S->A[hi - 1]);
	return i;
}

static void sparse_qsort(struct sparse *restrict S, int lo, int hi)
{
	if (hi - lo < 2) {
		return;
	}

	int mid = sparse_partition(S, lo, hi);
	sparse_qsort(S, lo, mid);
	sparse_qsort(S, mid + 1, hi);
}

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

float vec_dot(const struct vec *a, const struct vec *b)
{
	float result = 0;
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
void vec_scale(float scalar, struct vec *restrict v)
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

void vec2_scale(float scalar, struct vec2 *restrict v)
{
	for (int i = 0; i < 2; i++) {
		v->x[i] *= scalar;
	}
}

float vec2_dot(struct vec2 *a, struct vec2 *b)
{
	float result = 0;
	for (int i = 0; i < 2; i++) {
		result += a->x[i] * b->x[i];
	}
	return result;
}

void vec2_normalize(struct vec2 *restrict v)
{
	float norm = sqrtf(vec2_dot(v, v));
	vec2_scale(1.0f/norm, v);
}

void vec2_midpoint(struct vec2 *a, struct vec2 *b, struct vec2 *out)
{
	vec2_add(a, b, out);
	vec2_scale(0.5, out);
}

static float minor(float *restrict matrix, int dim, int i, int j)
{
	float *restrict sub = malloc((dim - 1) * (dim - 1) * sizeof(*sub));
	if (sub == NULL) {
		raise(SIGSEGV);
	}

	for (int ii = 0, iii = 0; ii < dim; ii++) {
		if (ii == i) {
			continue;
		}
		for (int jj = 0, jjj = 0; jj < dim; jj++) {
			if (jj == j) {
				continue;
			}
			sub[iii*(dim - 1) + jjj] = matrix[ii*dim + jj];
			jjj++;
		}
		iii++;
	}

	float sub_det = det(sub, dim - 1);
	free(sub);
	return sub_det;
}

float det(float *restrict matrix, int dim)
{
	if (dim == 1) {
		return *matrix;
	}

	int sgn = 1;
	float result = 0;
	for (int j = 0; j < dim; j++) {
		result += sgn * matrix[j] * minor(matrix, dim, 0, j);
		sgn *= -1;
	}
	return result;
}

void get_inverse(float *restrict matrix, int dim, float *restrict inverse)
{
	float inv_det = 1.0f / det(matrix, dim);

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			int sgn = 1;
			if ((i + j) % 2 != 0) {
				sgn = -1;
			}
			inverse[i*dim + j] = sgn * minor(matrix, dim, j, i) * inv_det;
		}
	}
}

int sparse_conj_grad(struct sparse *restrict A, struct vec *restrict b,
	struct vec *restrict c, float tolerance, struct vis *restrict vis,
	struct mesh *restrict mesh)
{
	float bsquared;
	float alpha, beta, old_r2, dAd;
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
		mesh_scalar_stress(mesh, c);
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
