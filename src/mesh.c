/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @brief mesh
 */

#include <stdlib.h>
#include <math.h>

#include "mesh.h"
#include "container_of.h"

static struct element_vtable triangle_vtable = {
	.stiffness_add = stiffness_add_triangle,
	.forces_add = forces_add_triangle,
	.scalar_stress = triangle_scalar_stress,
};

int mesh_init(struct mesh *restrict mesh)
{
	mesh->nvertices = 0;
	mesh->vertices = NULL;
	mesh->nelements = 0;
	mesh->elements = NULL;
	mesh->nenabled = 0;
	return 0;
}

void mesh_destroy(struct mesh *restrict mesh)
{
	if (mesh->vertices != NULL) {
		free(mesh->vertices);
		mesh->vertices = NULL;
	}
	if (mesh->elements != NULL) {
		for (int i = 0; i < mesh->nelements; i++) {
			free(mesh->elements[i]);
		}
		free(mesh->elements);
		mesh->elements = NULL;
	}
}

struct vertex *mesh_add_vertex(struct mesh *restrict mesh, float x, float y, bool enabled)
{
	mesh->nvertices++;
	mesh->vertices = realloc(mesh->vertices, mesh->nvertices * sizeof(*mesh->vertices));
	if (mesh->vertices == NULL) {
		return NULL;
	}
	int idx = mesh->nvertices - 1;

	mesh->vertices[idx].pos.x[0] = x;
	mesh->vertices[idx].pos.x[1] = y;
	mesh->vertices[idx].enabled = enabled;
	if (enabled) {
		mesh->nenabled++;
	}

	return &mesh->vertices[idx];
}

struct triangle *mesh_add_triangle(struct mesh *restrict mesh, struct vertex *v0,
	struct vertex *v1, struct vertex *v2, float density, float elasticity)
{
	mesh->nelements++;
	mesh->elements = realloc(mesh->elements, mesh->nelements * sizeof(*mesh->elements));
	if (mesh->elements == NULL) {
		return NULL;
	}

	struct triangle *triangle = malloc(sizeof(*triangle));
	if (triangle == NULL) {
		return NULL;
	}
	mesh->elements[mesh->nelements - 1] = &triangle->element;

	triangle->element.vtable = &triangle_vtable;
	triangle->element.mesh = mesh;

	triangle->element.nvertices = 3;
	triangle->element.vertices[0] = v0;
	triangle->element.vertices[1] = v1;
	triangle->element.vertices[2] = v2;

	triangle->element.density = density;
	triangle->element.elasticity = elasticity;

	triangle_compute_area(triangle);
	triangle_compute_dof(triangle);

	return triangle;
}

void mesh_assign_vertex_ids(struct mesh *restrict mesh)
{
	int next_id = 0;

	for (int i = 0; i < mesh->nvertices; i++) {
		if (mesh->vertices[i].enabled) {
			mesh->vertices[i].id = next_id;
			next_id++;
		} else {
			mesh->vertices[i].id = -1;
		}
	}
}

void triangle_compute_area(struct triangle *restrict triangle)
{
	struct vec2 v10, v20;

	struct vec2 *v0 = &triangle->element.vertices[0]->pos;
	struct vec2 *v1 = &triangle->element.vertices[1]->pos;
	struct vec2 *v2 = &triangle->element.vertices[2]->pos;

	vec2_sub(v1, v0, &v10);
	vec2_sub(v2, v0, &v20);
	triangle->area = fabsf(0.5f * (v10.x[0]*v20.x[1] - v10.x[1]*v20.x[0]));
}

void triangle_compute_dof(struct triangle *restrict triangle)
{
	struct vec2 v01, v21;

	for (int i = 0; i < 3; i++) {
		struct vec2 *v0 = &triangle->element.vertices[(i+0)%3]->pos;
		struct vec2 *v1 = &triangle->element.vertices[(i+1)%3]->pos;
		struct vec2 *v2 = &triangle->element.vertices[(i+2)%3]->pos;

		// gram schmidt
		vec2_sub(v0, v1, &v01);
		vec2_sub(v2, v1, &v21);
		vec2_scale(vec2_dot(&v01, &v21) / vec2_dot(&v21, &v21), &v21);
		vec2_sub(&v01, &v21, &triangle->dof_grad[(i+0)%3]);

		float s = 1.0f / vec2_dot(&v01, &triangle->dof_grad[(i+0)%3]);
		vec2_scale(s, &triangle->dof_grad[(i+0)%3]);
	}
}

float triangle_dof(struct triangle *restrict triangle, int vertex, float x, float y)
{
	struct vec2 vxy;
	vxy.x[0] = x;
	vxy.x[1] = y;
	return 1 + vec2_dot(&triangle->dof_grad[vertex], &vxy) - vec2_dot(&triangle->dof_grad[vertex], &triangle->element.vertices[vertex]->pos);
}

/*
 * integral of (E/2) * (D_i u_j D_i u_j + D_i u_i + D_j u_j)
 */
void stiffness_add_triangle(struct sparse *restrict A, struct element *restrict element)
{
	struct triangle *triangle = container_of(element, struct triangle, element);

	/* triangle vertex loops */
	for (int m = 0; m < 3; m++) {
		struct vertex *vm = element->vertices[m];
		if (!vm->enabled) {
			continue;
		}

		for (int n = 0; n < 3; n++) {
			struct vertex *vn = element->vertices[n];
			if (!vn->enabled) {
				continue;
			}

			int i = vm->id;
			int j = vn->id;

			/* xy loops */
			for (int r = 0; r < DIM; r++) {
				for (int s = 0; s < DIM; s++) {
					float entry = 0;
					if (r == s) {
						entry += vec2_dot(&triangle->dof_grad[m], &triangle->dof_grad[n]);
					}
					entry += triangle->dof_grad[m].x[r] * triangle->dof_grad[n].x[s];

					entry *= triangle->area * element->elasticity / 2;
					sparse_add(A, DIM*i + r, DIM*j + s, entry);
				}
			}
		}
	}
}

void forces_add_triangle(struct vec *restrict b, struct element *restrict element)
{
	struct triangle *triangle = container_of(element, struct triangle, element);

	float third_weight = triangle->area * element->density / 3;
	for (int n = 0; n < 3; n++) {
		struct vertex *v = element->vertices[n];
		if (v->enabled) {
			int i = v->id;
			b->x[DIM*i + 1] += -third_weight;
		}
	}
}

void triangle_scalar_stress(struct vec *restrict c, struct element *restrict element)
{
	struct triangle *triangle = container_of(element, struct triangle, element);
	float sxx = 0;
	float sxy = 0;
	float syy = 0;

	for (int m = 0; m < 3; m++) {
		if (!element->vertices[m]->enabled) {
			continue;
		}

		int i = element->vertices[m]->id;

		sxx += c->x[2*i] * triangle->dof_grad[m].x[0];
		syy += c->x[2*i + 1] * triangle->dof_grad[m].x[1];
		sxy += 0.5f * (c->x[2*i] * triangle->dof_grad[m].x[1] + c->x[2*i + 1] * triangle->dof_grad[m].x[0]);
	}

	sxx *= element->elasticity;
	sxy *= element->elasticity;
	syy *= element->elasticity;

	float pressure = -0.5f * (sxx + syy);
	sxx += pressure;
	syy += pressure;

	element->scalar_stress = sqrtf(SQR(sxx) + 2*SQR(sxy) + SQR(syy));
}

void mesh_construct_problem(struct mesh *restrict mesh, struct sparse *restrict A, struct vec *b)
{
	for (int i = 0; i < mesh->nelements; i++) {
		mesh->elements[i]->vtable->stiffness_add(A, mesh->elements[i]);
		mesh->elements[i]->vtable->forces_add(b, mesh->elements[i]);
	}
}

void mesh_scalar_stress(struct mesh *restrict mesh, struct vec *restrict c)
{
#ifdef _OPENMP
#pragma omp parallel for num_threads(OMP_NTHREAD)
#endif
	for (int i = 0; i < mesh->nelements; i++) {
		mesh->elements[i]->vtable->scalar_stress(c, mesh->elements[i]);
	}
}
