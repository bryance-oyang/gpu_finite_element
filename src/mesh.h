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

#define ELEMENT_MAX_VERTICES 10
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
	/* 0,1 are defining vertices, but also can contain middle points */
	int vertices[4];
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

	double density;
	double elasticity;
};

struct element_vtable {
	void (*stiffness_add)(struct sparse *restrict A, struct element *restrict element);
	void (*forces_add)(struct vec *restrict b, struct element *restrict element);
	double (*scalar_stress)(struct vec *restrict c, struct element *restrict element, double *x);
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

	double area;
	struct vec2 dof_grad[3];

	double scalar_stress;
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

	double jacob;
	/* inv_J[2*i + j] = dr^i / dx^j */
	double inv_J[4];
};

/**
 *    2
 *    |.
 *   7| .6
 *    |  .
 *   8| 9 .5
 *    |    .
 *   0|-----1
 *      3 4
 */
struct triangle3 {
	struct element element;

	double jacob;
	/* inv_J[2*i + j] = dr^i / dx^j */
	double inv_J[4];
};

int mesh_init(struct mesh *restrict mesh);
void mesh_destroy(struct mesh *restrict mesh);
int mesh_add_vertex(struct mesh *restrict mesh, double x, double y, bool enabled);
struct vertex *get_vert(struct element *restrict element, int vidx);
struct edge *mesh_add_edge(struct mesh *restrict mesh, struct vertex *v0, struct vertex *v1);
void mesh_assign_vertex_ids(struct mesh *restrict mesh);
void mesh_construct_problem(struct mesh *restrict mesh, struct sparse *restrict A, struct vec *b);

struct triangle *mesh_add_triangle(struct mesh *restrict mesh, int v0,
	int v1, int v2, double density, double elasticity);
void triangle_stiffness_add(struct sparse *restrict A, struct element *restrict element);
void triangle_forces_add(struct vec *restrict b, struct element *restrict element);
double triangle_scalar_stress(struct vec *restrict c, struct element *restrict element, double *x);

struct triangle2 *mesh_add_triangle2(struct mesh *restrict mesh, int v0,
	int v1, int v2, double density, double elasticity);
void triangle2_stiffness_add(struct sparse *restrict A, struct element *restrict element);
void triangle2_forces_add(struct vec *restrict b, struct element *restrict element);
double triangle2_scalar_stress(struct vec *restrict c, struct element *restrict element, double *x);

struct triangle3 *mesh_add_triangle3(struct mesh *restrict mesh, int v0,
	int v1, int v2, double density, double elasticity);
void triangle3_stiffness_add(struct sparse *restrict A, struct element *restrict element);
void triangle3_forces_add(struct vec *restrict b, struct element *restrict element);
double triangle3_scalar_stress(struct vec *restrict c, struct element *restrict element, double *x);

#endif /* MESH_H */
