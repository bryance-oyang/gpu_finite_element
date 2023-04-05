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

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <stdint.h>

#include "macro_def.h"

struct mesh;
struct vis;

/** A_{row[k], col[k]} = A[k] */
struct sparse {
	/** number of nonzero matrix elements */
	int len;
	/** allocated memory in bytes */
	int size;
	int32_t *row;
	int32_t *col;
	float *A;
};

int sparse_init(struct sparse *restrict S);
void sparse_destroy(struct sparse *restrict S);
int sparse_add(struct sparse *restrict S, int row, int col, float entry);
void sparse_sort(struct sparse *restrict S);
void sparse_consolidate(struct sparse *restrict S);

struct vec {
	int dim;
	float *x;
};

int vec_init(struct vec *restrict v, int dim);
void vec_destroy(struct vec *restrict v);

void sparse_mult_vec(const struct sparse *restrict S, const struct vec *restrict v, struct vec *restrict out);
float vec_dot(const struct vec *a, const struct vec *b);
void vec_copy(const struct vec *restrict in, struct vec *restrict out);
void vec_add(struct vec *a, struct vec *b, struct vec *out);
void vec_sub(struct vec *a, struct vec *b, struct vec *out);
void vec_scale(float scalar, struct vec *restrict v);

struct vec2 {
	float x[2];
};

void vec2_copy(struct vec2 *restrict in, struct vec2 *restrict out);
void vec2_add(struct vec2 *a, struct vec2 *b, struct vec2 *out);
void vec2_sub(struct vec2 *a, struct vec2 *b, struct vec2 *out);
void vec2_scale(float scalar, struct vec2 *restrict v);
float vec2_dot(struct vec2 *a, struct vec2 *b);
void vec2_normalize(struct vec2 *restrict v);

int sparse_conj_grad(struct sparse *restrict A, struct vec *restrict b,
	struct vec *restrict c, float tolerance, struct vis *restrict vis,
	struct mesh *restrict mesh);

#endif /* LINEAR_ALGEBRA_H */
