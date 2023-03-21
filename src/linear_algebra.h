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

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#if SINGLE_PRECISION == 1
typedef float number;
#else
typedef double number;
#endif

struct sparse {
	int len;
	int size;
	int *row;
	int *col;
	number *A;
};

struct vec {
	int dim;
	number *x;
};

int sparse_init(struct sparse *restrict S);
void sparse_destroy(struct sparse *restrict S);
int sparse_get_idx(struct sparse *restrict S, int row, int col, int *found);
int sparse_add(struct sparse *restrict S, int row, int col, number entry);

int vec_init(struct vec *restrict v, int dim);
void vec_destroy(struct vec *restrict v);

void sparse_mult_vec(const struct sparse *restrict S, const struct vec *restrict v, struct vec *restrict out);
number vec_dot(const struct vec *a, const struct vec *b);
number vec_S_dot(const struct vec *a, const struct sparse *restrict S, const struct vec *b);
void vec_copy(const struct vec *restrict in, struct vec *restrict out);
void vec_add(struct vec *a, struct vec *b, struct vec *out);
void vec_sub(struct vec *a, struct vec *b, struct vec *out);
void vec_scale(number scalar, struct vec *restrict v);

void sparse_conj_grad(number tolerance, const struct sparse *restrict S, const struct vec *restrict b, struct vec *restrict out);

#endif /* LINEAR_ALGEBRA_H */
