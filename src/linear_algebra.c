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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linear_algebra.h"

/*
 * comparison function to keep sparse matrix sorted:
 * (row, col) < (row, col)? in column major format (same as CUDA)
 */
static int sparse_row_col_cmp(int r0, int c0, int r1, int c1)
{
	if (c0 != c1) {
		return (c0 > c1) - (c0 < c1);
	} else {
		return (r0 > r1) - (r0 < r1);
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

/* returns index of element at (row,col) or the index it should be inserted at */
int sparse_get_idx(struct sparse *restrict S, int row, int col, int *found)
{
	int lo, mid, hi, cmp;

	if (S->len == 0) {
		if (found != NULL) {
			*found = 0;
		}
		return 0;
	}

	lo = 0;
	hi = S->len;
	while (hi - lo > 1) {
		mid = (lo + hi) / 2;

		cmp = sparse_row_col_cmp(row, col, S->row[mid], S->col[mid]);
		if (cmp < 0) {
			hi = mid;
		} else if (cmp > 0) {
			lo = mid;
		} else {
			if (found != NULL) {
				*found = 1;
			}
			return mid;
		}
	}
	if (sparse_row_col_cmp(row, col, S->row[lo], S->col[lo]) == 0) {
		if (found != NULL) {
			*found = 1;
		}
		return lo;
	} else {
		if (found != NULL) {
			*found = 0;
		}
		return lo + 1;
	}
}

/* insert nonexistent */
static int sparse_insert(struct sparse *restrict S, int idx, int row,
	int col, number entry)
{
	if (S->len == S->size) {
		S->size *= 2;
		S->row = reallocarray(S->row, S->size, sizeof(*S->row));
		S->col = reallocarray(S->col, S->size, sizeof(*S->col));
		S->A = reallocarray(S->A, S->size, sizeof(*S->A));

		if (S->row == NULL || S->col == NULL || S->A == NULL) {
			return -1;
		}
	}

	memmove(&S->row[idx + 1], &S->row[idx], (S->len - idx) * sizeof(*S->row));
	memmove(&S->col[idx + 1], &S->col[idx], (S->len - idx) * sizeof(*S->col));
	memmove(&S->A[idx + 1], &S->A[idx], (S->len - idx) * sizeof(*S->A));
	S->row[idx] = row;
	S->col[idx] = col;
	S->A[idx] = entry;
	S->len++;
	return 0;
}

int sparse_add(struct sparse *restrict S, int row, int col, number entry)
{
	int found, idx;

	idx = sparse_get_idx(S, row, col, &found);
	if (found) {
		S->A[idx] += entry;
		return 0;
	} else {
		return sparse_insert(S, idx, row, col, entry);
	}
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

number vec_dot(const struct vec *a, const struct vec *b)
{
	number result = 0;
	int dim = a->dim;

	for (int i = 0; i < dim; i++) {
		result += a->x[i] * b->x[i];
	}
	return result;
}

number vec_S_dot(const struct vec *a, const struct sparse *restrict S, const struct vec *b)
{
	number result = 0;

	for (int n = 0; n < S->len; n++) {
		int i = S->row[n];
		int j = S->col[n];
		result += a->x[i] * S->A[n] * b->x[j];
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
void vec_scale(number scalar, struct vec *restrict v)
{
	for (int i = 0; i < v->dim; i++) {
		v->x[i] *= scalar;
	}
}

void sparse_conj_grad(const struct sparse *restrict S,
	const struct vec *restrict b, struct vec *restrict out, number tolerance)
{
	number bsquared;
	number alpha, beta, old_r2;
	struct vec r, d, tmp, A_alpha_d;

	bsquared = vec_dot(b, b);

	vec_init(&r, out->dim);
	vec_init(&d, out->dim);
	vec_init(&tmp, out->dim);
	vec_init(&A_alpha_d, out->dim);

	vec_scale(0, out);
	vec_copy(b, &r);
	vec_copy(b, &d);

	for (int k = 1; k <= out->dim; k++) {
		old_r2 = vec_dot(&r, &r);
		if (bsquared > 0 && old_r2 / bsquared <= tolerance) {
			break;
		} else if (bsquared == 0 && old_r2 <= tolerance) {
			break;
		}

		number dSd = vec_S_dot(&d, S, &d);
		alpha = old_r2 / dSd;

		vec_copy(&d, &tmp);
		vec_scale(alpha, &tmp);
		vec_add(out, &tmp, out);

		sparse_mult_vec(S, &tmp, &A_alpha_d);
		vec_sub(&r, &A_alpha_d, &r);

		beta = vec_dot(&r, &r) / old_r2;
		vec_scale(beta / alpha, &tmp);
		vec_add(&r, &tmp, &d);
	}

	vec_destroy(&r);
	vec_destroy(&d);
	vec_destroy(&tmp);
	vec_destroy(&A_alpha_d);
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

void vec2_scale(number scalar, struct vec2 *restrict v)
{
	for (int i = 0; i < 2; i++) {
		v->x[i] *= scalar;
	}
}

number vec2_dot(struct vec2 *a, struct vec2 *b)
{
	number result = 0;
	for (int i = 0; i < 2; i++) {
		result += a->x[i] * b->x[i];
	}
	return result;
}

void vec2_normalize(struct vec2 *restrict v)
{
	number norm = sqrt(vec2_dot(v, v));
	vec2_scale(1.0/norm, v);
}
