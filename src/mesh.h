/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef MESH_H
#define MESH_H

#include <stdbool.h>
#include "linear_algebra.h"

#define ELEMENT_MAX_VERTICES 8

struct vertex {
	struct vec2 pos;
	int id;
	bool enabled;
};

struct element {
	struct element_vtable *vtable;
	struct mesh *mesh;

	int nvertices;
	struct vertex *vertices[ELEMENT_MAX_VERTICES];

	float density;
	float elasticity;
	float scalar_stress;
};

struct element_vtable {
	void (*stiffness_add)(struct sparse *restrict A, struct element *restrict element);
	void (*forces_add)(struct vec *restrict b, struct element *restrict element);
	void (*scalar_stress)(struct vec *restrict c, struct element *restrict element);
};

struct triangle {
	struct element element;
	float area;
	struct vec2 dof_grad[3];
};

struct mesh {
	int nvertices;
	int nenabled;
	struct vertex *vertices;

	int nelements;
	struct element **elements;
};

int mesh_init(struct mesh *restrict mesh);
void mesh_destroy(struct mesh *restrict mesh);
struct vertex *mesh_add_vertex(struct mesh *restrict mesh, float x, float y, bool enabled);
struct triangle *mesh_add_triangle(struct mesh *restrict mesh, struct vertex *v0,
	struct vertex *v1, struct vertex *v2, float density, float elasticity);
void mesh_assign_vertex_ids(struct mesh *restrict mesh);

void triangle_compute_area(struct triangle *restrict triangle);
void triangle_compute_dof(struct triangle *restrict triangle);
float triangle_dof(struct triangle *restrict triangle, int vertex, float x, float y);

void stiffness_add_triangle(struct sparse *restrict A, struct element *restrict element);
void forces_add_triangle(struct vec *restrict b, struct element *restrict element);
void triangle_scalar_stress(struct vec *restrict c, struct element *restrict element);

void mesh_construct_problem(struct mesh *restrict mesh, struct sparse *restrict A, struct vec *b);

#endif /* MESH_H */
