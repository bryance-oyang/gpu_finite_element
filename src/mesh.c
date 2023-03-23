/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <stdlib.h>
#include <math.h>
#include "mesh.h"

int mesh_init(struct mesh *restrict mesh)
{
	mesh->nvertices = 0;
	mesh->vertices = NULL;
	mesh->ntriangles = 0;
	mesh->triangles = NULL;
	mesh->nenabled = 0;
	return 0;
}

void mesh_destroy(struct mesh *restrict mesh)
{
	if (mesh->vertices != NULL) {
		free(mesh->vertices);
		mesh->vertices = NULL;
	}
	if (mesh->triangles != NULL) {
		free(mesh->triangles);
		mesh->triangles = NULL;
	}
}

int mesh_add_vertex(struct mesh *restrict mesh, number x, number y, bool enabled)
{
	mesh->nvertices++;
	mesh->vertices = reallocarray(mesh->vertices, mesh->nvertices, sizeof(*mesh->vertices));
	if (mesh->vertices == NULL) {
		return -1;
	}
	int idx = mesh->nvertices - 1;

	mesh->vertices[idx].pos.x[0] = x;
	mesh->vertices[idx].pos.x[1] = y;
	mesh->vertices[idx].enabled = enabled;
	if (enabled) {
		mesh->nenabled++;
	}

	return idx;
}

int mesh_add_triangle(struct mesh *restrict mesh, int v0, int v1, int v2, number density, number elasticity)
{
	mesh->ntriangles++;
	mesh->triangles = reallocarray(mesh->triangles, mesh->ntriangles, sizeof(*mesh->triangles));
	if (mesh->triangles == NULL) {
		return -1;
	}
	int idx = mesh->ntriangles - 1;

	mesh->triangles[idx].mesh = mesh;

	mesh->triangles[idx].vertices[0] = v0;
	mesh->triangles[idx].vertices[1] = v1;
	mesh->triangles[idx].vertices[2] = v2;

	mesh->triangles[idx].density = density;
	mesh->triangles[idx].elasticity = elasticity;

	return idx;
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

struct triangle triangle_make(int v0, int v1, int v2, number density, number elasticity)
{
	struct triangle triangle;
	triangle.vertices[0] = v0;
	triangle.vertices[1] = v1;
	triangle.vertices[2] = v2;
	triangle.density = density;
	triangle.elasticity = elasticity;
	return triangle;
}

void triangle_compute_area(struct triangle *restrict triangle)
{
	struct vec2 v10, v20;
	struct vec2 *v0 = &triangle->mesh->vertices[triangle->vertices[0]].pos;
	struct vec2 *v1 = &triangle->mesh->vertices[triangle->vertices[1]].pos;
	struct vec2 *v2 = &triangle->mesh->vertices[triangle->vertices[2]].pos;

	vec2_sub(v1, v0, &v10);
	vec2_sub(v2, v0, &v20);
	triangle->area = fabs(0.5 * (v10.x[0]*v20.x[1] - v10.x[1]*v20.x[0]));
}

void triangle_compute_dof(struct triangle *restrict triangle)
{
	struct vec2 v01, v21;

	for (int i = 0; i < 3; i++) {
		struct vec2 *v0 = &triangle->mesh->vertices[triangle->vertices[(i+0)%3]].pos;
		struct vec2 *v1 = &triangle->mesh->vertices[triangle->vertices[(i+1)%3]].pos;
		struct vec2 *v2 = &triangle->mesh->vertices[triangle->vertices[(i+2)%3]].pos;

		// gram schmidt
		vec2_sub(v0, v1, &v01);
		vec2_sub(v2, v1, &v21);
		vec2_scale(vec2_dot(&v01, &v21) / vec2_dot(&v21, &v21), &v21);
		vec2_sub(&v01, &v21, &triangle->dof_grad[(i+0)%3]);

		number s = 1.0 / vec2_dot(&v01, &triangle->dof_grad[(i+0)%3]);
		vec2_scale(s, &triangle->dof_grad[(i+0)%3]);
	}
}

number triangle_dof(struct triangle *restrict triangle, int vertex, number x, number y)
{
	struct vec2 vxy;
	vxy.x[0] = x;
	vxy.x[1] = y;
	return 1 + vec2_dot(&triangle->dof_grad[vertex], &vxy) - vec2_dot(&triangle->dof_grad[vertex], &triangle->mesh->vertices[triangle->vertices[vertex]].pos);
}

/*
 * integral of (E/2) * (D_i u_j D_i u_j + D_i u_i + D_j u_j)
 */
void stiffness_add_triangle(struct sparse *restrict A, struct triangle *restrict triangle)
{
	/* triangle vertex loops */
	for (int m = 0; m < 3; m++) {
		struct vertex *vm = &triangle->mesh->vertices[triangle->vertices[m]];
		if (!vm->enabled) {
			continue;
		}

		for (int n = 0; n < 3; n++) {
			struct vertex *vn = &triangle->mesh->vertices[triangle->vertices[n]];
			if (!vn->enabled) {
				continue;
			}

			int i = vm->id;
			int j = vn->id;

			/* xy loops */
			for (int r = 0; r < DIM; r++) {
				for (int s = 0; s < DIM; s++) {
					number entry = 0;
					if (r == s) {
						entry += vec2_dot(&triangle->dof_grad[m], &triangle->dof_grad[n]);
					}
					entry += triangle->dof_grad[m].x[r] * triangle->dof_grad[n].x[s];

					entry *= triangle->area * triangle->elasticity / 2;
					sparse_add(A, DIM*i + r, DIM*j + s, entry);
				}
			}
		}
	}
}

void forces_add_triangle(struct vec *restrict b, struct triangle *restrict triangle)
{
	number third_weight = triangle->area * triangle->density / 3;
	for (int n = 0; n < 3; n++) {
		struct vertex *v = &triangle->mesh->vertices[triangle->vertices[n]];
		if (v->enabled) {
			int i = v->id;
			b->x[DIM*i + 1] += -third_weight;
		}
	}
}

void mesh_construct_problem(struct mesh *restrict mesh, struct sparse *restrict A, struct vec *b)
{
	for (int i = 0; i < mesh->ntriangles; i++) {
		triangle_compute_area(&mesh->triangles[i]);
		triangle_compute_dof(&mesh->triangles[i]);
		stiffness_add_triangle(A, &mesh->triangles[i]);
		forces_add_triangle(b, &mesh->triangles[i]);
	}
}
