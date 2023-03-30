/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <unistd.h>
#include <math.h>
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

	if (cuda_init(p) != 0) {
		goto err_nocuda;
	}

	return 0;

err_nocuda:
	vec_destroy(&p->c);
err_noc:
	vec_destroy(&p->b);
err_nob:
	sparse_destroy(&p->A);
err_noA:
	return -1;
}

void fep_destroy(struct finite_element_problem *restrict p)
{
	cuda_destroy(p);
	vec_destroy(&p->c);
	vec_destroy(&p->b);
	sparse_destroy(&p->A);
}

void sparse_conj_grad(struct finite_element_problem *restrict p,
	number tolerance, struct vis *vis)
{
	struct sparse *restrict A = &p->A;
	struct vec *restrict b = &p->b;
	struct vec *restrict c = &p->c;

	number bsquared;
	number alpha, beta, old_r2;
	struct vec r, d, tmp, A_alpha_d;

	bsquared = vec_dot(b, b);

	vec_init(&r, c->dim);
	vec_init(&d, c->dim);
	vec_init(&tmp, c->dim);
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

		number dSd = vec_S_dot(&d, A, &d);
		alpha = old_r2 / dSd;

		vec_copy(&d, &tmp);
		vec_scale(alpha, &tmp);
		vec_add(c, &tmp, c);

		sparse_mult_vec(A, &tmp, &A_alpha_d);
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

void fep_solve(struct finite_element_problem *restrict p, number tolerance, struct vis *vis)
{
#ifdef GPU_COMPUTE
#else /* GPU_COMPUTE */
	sparse_conj_grad(p, tolerance, vis);
#endif /* GPU_COMPUTE */
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
