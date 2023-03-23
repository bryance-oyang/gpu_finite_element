/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "finite_element.h"

int finite_element_problem_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh)
{
	mesh_assign_vertex_ids(mesh);

	if (sparse_init(&p->A) != 0) {
		goto err_noA;
	}
	if (vec_init(&p->b, DIM * mesh->nenabled) != 0) {
		goto err_nob;
	}
	if (vec_init(&p->c, DIM * mesh->nenabled) != 0) {
		goto err_noc;
	}

	mesh_construct_problem(mesh, &p->A, &p->b);
	return 0;

err_noc:
	vec_destroy(&p->b);
err_nob:
	sparse_destroy(&p->A);
err_noA:
	return -1;
}

void finite_element_problem_destroy(struct finite_element_problem *restrict p)
{
	vec_destroy(&p->c);
	vec_destroy(&p->b);
	sparse_destroy(&p->A);
}

void finite_element_problem_solve(struct finite_element_problem *restrict p, number tolerance)
{
	sparse_conj_grad(&p->A, &p->b, &p->c, tolerance);
}
