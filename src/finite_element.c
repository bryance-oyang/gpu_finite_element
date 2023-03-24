/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <math.h>
#include "finite_element.h"

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
	return 0;

err_noc:
	vec_destroy(&p->b);
err_nob:
	sparse_destroy(&p->A);
err_noA:
	return -1;
}

void fep_destroy(struct finite_element_problem *restrict p)
{
	vec_destroy(&p->c);
	vec_destroy(&p->b);
	sparse_destroy(&p->A);
}

void fep_solve(struct finite_element_problem *restrict p, number tolerance)
{
	sparse_conj_grad(&p->A, &p->b, &p->c, tolerance);
}

static void triangle_scalar_stress(struct vec *restrict c, struct triangle *restrict triangle, struct vertex *restrict vertices)
{
	number sxx = 0;
	number sxy = 0;
	number syy = 0;

	for (int m = 0; m < 3; m++) {
		if (!vertices[triangle->vertices[m]].enabled) {
			continue;
		}

		int i = vertices[triangle->vertices[m]].id;

		sxx += c->x[2*i] * triangle->dof_grad[m].x[0];
		syy += c->x[2*i + 1] * triangle->dof_grad[m].x[1];
		sxy += 0.5 * (c->x[2*i] * triangle->dof_grad[m].x[1] + c->x[2*i + 1] * triangle->dof_grad[m].x[0]);
	}

	sxx *= triangle->elasticity;
	sxy *= triangle->elasticity;
	syy *= triangle->elasticity;

	triangle->scalar_stress = sqrt(SQR(sxx) + 2*SQR(sxy) + SQR(syy));
}

void fep_scalar_stress(struct finite_element_problem *restrict p)
{
	struct mesh *mesh = p->mesh;
	for (int i = 0; i < mesh->ntriangles; i++) {
		triangle_scalar_stress(&p->c, &mesh->triangles[i], mesh->vertices);
	}
}
