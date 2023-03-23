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

#define DIM 2

struct vertex {
	struct vec2 pos;
	int id;
	bool enabled;
};

struct triangle {
	struct mesh *mesh;
	int vertices[3];
	number density;
	number elasticity;

	number area;
	struct vec2 dof_grad[3];
};

struct mesh {
	int nvertices;
	int nenabled;
	struct vertex *vertices;

	int ntriangles;
	struct triangle *triangles;
};

int mesh_init(struct mesh *restrict mesh);
void mesh_destroy(struct mesh *restrict mesh);
int mesh_add_vertex(struct mesh *restrict mesh, number x, number y, bool enabled);
int mesh_add_triangle(struct mesh *restrict mesh, int v0, int v1, int v2, number density, number elasticity);
void mesh_assign_vertex_ids(struct mesh *restrict mesh);

void triangle_compute_area(struct triangle *restrict triangle);
void triangle_compute_dof(struct triangle *restrict triangle);
number triangle_dof(struct triangle *restrict triangle, int vertex, number x, number y);

void stiffness_add_triangle(struct sparse *restrict A, struct triangle *restrict triangle);
void forces_add_triangle(struct vec *restrict b, struct triangle *restrict triangle);

void mesh_construct_problem(struct mesh *restrict mesh, struct sparse *restrict A, struct vec *b);

#endif /* MESH_H */
