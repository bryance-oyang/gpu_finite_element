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

#ifndef MESH_H
#define MESH_H

#include <stdbool.h>
#include "hash_table.h"
#include "linear_algebra.h"

#define ELEMENT_MAX_VERTICES 8
#define ELEMENT_MAX_EDGES 12
#define ELEMENT_MAX_FACES 8

struct vertex {
	struct vec2 pos;
	int id;
	bool enabled;
};

struct edge {
	struct ht_node node;

	bool is_boundary;
	/* 0,1 are defining vertices, but also can contain midpoint */
	int vertices[3];
};

struct face {
	struct ht_node node;

	bool is_boundary;
	int vertices[4];
};

struct element {
	struct element_vtable *vtable;
	struct mesh *mesh;

	int nvertices;
	/* indices of vertices */
	int vertices[ELEMENT_MAX_VERTICES];
	int nedges;
	struct edge *edges[ELEMENT_MAX_EDGES];
	int nfaces;
	struct face *faces[ELEMENT_MAX_FACES];

	float density;
	float elasticity;
	float scalar_stress;
};

struct element_vtable {
	void (*stiffness_add)(struct sparse *restrict A, struct element *restrict element);
	void (*forces_add)(struct vec *restrict b, struct element *restrict element);
	void (*scalar_stress)(struct vec *restrict c, struct element *restrict element);
};

struct mesh {
	int nvertices;
	int vertices_size;
	int nenabled;
	struct vertex *vertices;

	struct ht edges_table;
	struct ht faces_table;

	int nelements;
	int elements_size;
	struct element **elements;
};

struct triangle {
	struct element element;

	float area;
	struct vec2 dof_grad[3];
};

/**
 *    2
 *    |.
 *    | .
 *   4|  .3
 *    |   .
 *    |    .
 *   0|-----1
 *       5
 *
 * The dof's are pulling each vertex in the x and y direction by 1 unit while
 * keeping all others fixed. For v: vertex, d: direction (xy), i: vector index,
 * u_{vdi} = f_v delta_{di} where f_v is the quadratic function 1 at v and 0 for
 * other vertices, determined by coefficients in canonical triangle
 */
struct triangle2 {
	struct element element;

	float jacob;
	/* J[i][j] = dx^i / dr^j */
	float J[2][2];
	/* inv_J[i][j] = dr^i / dx^j */
	float inv_J[2][2];
	/* normed J: sum_i dx^i/dr^j dx^i/dr^j = 1: make so that dof in xy frame is normed to 1 as opposed to coordinate basis */
	float N[2][2];
	float inv_N[2][2];
};

int mesh_init(struct mesh *restrict mesh);
void mesh_destroy(struct mesh *restrict mesh);
int mesh_add_vertex(struct mesh *restrict mesh, float x, float y, bool enabled);
struct vertex *get_vert(struct element *restrict element, int vidx);
struct edge *mesh_add_edge(struct mesh *restrict mesh, struct vertex *v0, struct vertex *v1);
void mesh_assign_vertex_ids(struct mesh *restrict mesh);
void mesh_construct_problem(struct mesh *restrict mesh, struct sparse *restrict A, struct vec *b);
void mesh_scalar_stress(struct mesh *restrict mesh, struct vec *restrict c);

struct triangle *mesh_add_triangle(struct mesh *restrict mesh, int v0,
	int v1, int v2, float density, float elasticity);
void triangle_stiffness_add(struct sparse *restrict A, struct element *restrict element);
void triangle_forces_add(struct vec *restrict b, struct element *restrict element);
void triangle_scalar_stress(struct vec *restrict c, struct element *restrict element);

struct triangle2 *mesh_add_triangle2(struct mesh *restrict mesh, int v0,
	int v1, int v2, float density, float elasticity);
void triangle2_stiffness_add(struct sparse *restrict A, struct element *restrict element);
void triangle2_forces_add(struct vec *restrict b, struct element *restrict element);
void triangle2_scalar_stress(struct vec *restrict c, struct element *restrict element);

#endif /* MESH_H */
