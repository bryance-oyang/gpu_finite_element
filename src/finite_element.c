/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include "finite_element.h"
#include "cuda.h"

int fep_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh)
{
	if (sparse_init(&p->A) != 0) {
		goto err_noA;
	}
	if (vec_init(&p->b, DIM * mesh->nenabled) != 0) {
		goto err_nob;
	}
	if (vec_init(&p->c, DIM * mesh->nenabled) != 0) {
		goto err_noc;
	}

	p->mesh = mesh;
	mesh_assign_vertex_ids(mesh);
	mesh_construct_problem(mesh, &p->A, &p->b);

#ifdef GPU_COMPUTE
	if (cuda_init(p) != 0) {
		goto err_nocuda;
	}
#endif /* GPU_COMPUTE */

	return 0;

#ifdef GPU_COMPUTE
err_nocuda:
	vec_destroy(&p->c);
#endif /* GPU_COMPUTE */

err_noc:
	vec_destroy(&p->b);
err_nob:
	sparse_destroy(&p->A);
err_noA:
	return -1;
}

void fep_destroy(struct finite_element_problem *restrict p)
{
#ifdef GPU_COMPUTE
	cuda_destroy(p);
#endif
	vec_destroy(&p->c);
	vec_destroy(&p->b);
	sparse_destroy(&p->A);
}

int sparse_conj_grad(struct finite_element_problem *restrict p,
	float tolerance, struct vis *vis)
{
	struct sparse *restrict A = &p->A;
	struct vec *restrict b = &p->b;
	struct vec *restrict c = &p->c;

	float bsquared;
	float alpha, beta, old_r2;
	struct vec r, d, A_alpha_d;

	bsquared = vec_dot(b, b);

	vec_init(&r, c->dim);
	vec_init(&d, c->dim);
	vec_init(&A_alpha_d, c->dim);

	vec_scale(0, c);
	vec_copy(b, &r);
	vec_copy(b, &d);

	for (int k = 1; k <= c->dim; k++) {
#ifdef ANIMATE
		fep_scalar_stress(p);
		vis_fill(vis, p->mesh);
		vis_send(vis);
		usleep(15000);
#else /* ANIMATE */
		(void)vis;
#endif /* ANIMATE */

		old_r2 = vec_dot(&r, &r);
		if (bsquared > 0 && old_r2 / bsquared <= tolerance) {
			break;
		} else if (bsquared == 0 && old_r2 <= tolerance) {
			break;
		}

		float dSd = vec_S_dot(&d, A, &d);
		alpha = old_r2 / dSd;

		vec_scale(alpha, &d);
		vec_add(c, &d, c);

		sparse_mult_vec(A, &d, &A_alpha_d);
		vec_sub(&r, &A_alpha_d, &r);

		beta = vec_dot(&r, &r) / old_r2;
		vec_scale(beta / alpha, &d);
		vec_add(&r, &d, &d);
	}

	vec_destroy(&r);
	vec_destroy(&d);
	vec_destroy(&A_alpha_d);

	return 0;
}

int fep_solve(struct finite_element_problem *restrict p, float tolerance, struct vis *vis)
{
	int retval;
	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);

#ifdef GPU_COMPUTE
	retval = gpu_conj_gradient(p, tolerance);
#else /* GPU_COMPUTE */
	retval = sparse_conj_grad(p, tolerance, vis);
#endif /* GPU_COMPUTE */

	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	float msec = (end.tv_sec - start.tv_sec)*1e3 + (end.tv_nsec - start.tv_nsec)*1e-6;
	printf("benchmarked solve time: %g ms\n", msec);

	return retval;
}

void fep_scalar_stress(struct finite_element_problem *restrict p)
{
	struct mesh *restrict mesh = p->mesh;
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
	for (int i = 0; i < mesh->nelements; i++) {
		mesh->elements[i]->vtable->scalar_stress(&p->c, mesh->elements[i]);
	}
}
