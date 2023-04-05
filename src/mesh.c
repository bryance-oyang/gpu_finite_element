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
#include <stdio.h>
#include <math.h>

#include "mesh.h"
#include "container_of.h"

#define HBUF_LEN 256
char hash_buf[HBUF_LEN];

static struct element_vtable triangle_vtable = {
	.stiffness_add = stiffness_add_triangle,
	.forces_add = forces_add_triangle,
	.scalar_stress = triangle_scalar_stress,
};

int mesh_init(struct mesh *restrict mesh)
{
	mesh->nvertices = 0;
	mesh->vertices_size = 128;
	mesh->vertices = malloc(mesh->vertices_size * sizeof(*mesh->vertices));
	if (mesh->vertices == NULL) {
		goto err_vertices;
	}

	mesh->nelements = 0;
	mesh->elements_size = 128;
	mesh->elements = malloc(mesh->elements_size * sizeof(*mesh->elements));
	if (mesh->elements == NULL) {
		goto err_elements;
	}
	mesh->nenabled = 0;

	if (ht_init(&mesh->edges_table) != 0) {
		goto err_edges_table;
	}

	if (ht_init(&mesh->faces_table) != 0) {
		goto err_faces_table;
	}

	return 0;

err_faces_table:
	ht_destroy(&mesh->edges_table);
err_edges_table:
	free(mesh->elements);
err_elements:
	free(mesh->vertices);
err_vertices:
	return -1;
}

void mesh_destroy(struct mesh *restrict mesh)
{
	free(mesh->vertices);

	for (int i = 0; i < mesh->nelements; i++) {
		free(mesh->elements[i]);
	}
	free(mesh->elements);

	ht_destroy(&mesh->edges_table);
	ht_destroy(&mesh->faces_table);
}

static inline int get_vertex_idx(struct mesh *mesh, struct vertex *v)
{
	return v - mesh->vertices;
}

struct vertex *mesh_add_vertex(struct mesh *restrict mesh, float x, float y, bool enabled)
{
	if (mesh->nvertices == mesh->vertices_size) {
		int size = mesh->vertices_size * 2;
		void *tmp = realloc(mesh->vertices, size * sizeof(*mesh->vertices));
		if (tmp == NULL) {
			return NULL;
		}

		mesh->vertices = tmp;
		mesh->vertices_size = size;
	}

	mesh->nvertices++;
	int idx = mesh->nvertices - 1;

	mesh->vertices[idx].pos.x[0] = x;
	mesh->vertices[idx].pos.x[1] = y;
	mesh->vertices[idx].enabled = enabled;
	if (enabled) {
		mesh->nenabled++;
	}

	return &mesh->vertices[idx];
}

static struct edge *alloc_new_edge()
{
	struct edge *edge = malloc(sizeof(*edge));
	if (edge == NULL) {
		return NULL;
	}
	ht_node_init(&edge->node);
	return edge;
}

/**
 * find edge given v0, v1. v0 and v1 will be modified so they are sorted by
 * increasing index (v0_idx < v1_idx). if pprev is not NULL, set to hash table
 * node so that (*pprev)->next is the position to insert edge if not found
 */
static struct edge *find_edge(struct mesh *restrict mesh,
	struct vertex **v0, struct vertex **v1, struct ht_node **pprev)
{
	/* order v0, v1 by index */
	int v0_idx = get_vertex_idx(mesh, *v0);
	int v1_idx = get_vertex_idx(mesh, *v1);
	if (v0_idx > v1_idx) {
		int tmp = v0_idx;
		v0_idx = v1_idx;
		v1_idx = tmp;

		struct vertex *tmp2 = *v0;
		*v0 = *v1;
		*v1 = tmp2;
	}

	/* key for hash table */
	if (snprintf(hash_buf, HBUF_LEN, "%d_%d", v0_idx, v1_idx) >= HBUF_LEN) {
		return NULL;
	}
	unsigned int key = knuth_hash(hash_buf);

	/* find edge in hash table: traverse linked list until equality match */
	struct ht_node *cur, *prev;
	struct edge *edge = NULL;
	ht_for_each(&mesh->edges_table, cur, prev, key) {
		struct edge *tmp = container_of(cur, struct edge, node);
		if (tmp->vertices[0] == *v0 && tmp->vertices[1] == *v1) {
			edge = tmp;
			break;
		}
	}

	if (pprev != NULL) {
		*pprev = prev;
	}
	return edge;
}

struct edge *mesh_add_edge(struct mesh *restrict mesh, struct vertex *v0, struct vertex *v1)
{
	struct ht_node *prev;
	struct edge *edge = find_edge(mesh, &v0, &v1, &prev);

	if (edge == NULL) {
		/* not found in ht */
		if ((edge = alloc_new_edge()) == NULL) {
			return NULL;
		}
		prev->next = &edge->node;

		edge->is_boundary = true;
		edge->vertices[0] = v0;
		edge->vertices[1] = v1;
	} else {
		/* cannot be boundary edge if another triangle shares edge */
		edge->is_boundary = false;
	}

	return edge;
}

/* assign matrix indices: only enabled vertices will enter the matrix */
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

	/*
	triangle->element.nedges = 3;
	triangle->element.edges[0] = mesh_add_edge(mesh, v0, v1);
	triangle->element.edges[1] = mesh_add_edge(mesh, v1, v2);
	triangle->element.edges[2] = mesh_add_edge(mesh, v2, v0);
	*/

	triangle_compute_area(triangle);
	triangle_compute_dof(triangle);

	return triangle;
}

void triangle_compute_area(struct triangle *restrict triangle)
{
	struct vec2 v10, v20;

	struct vec2 *v0 = &triangle->element.vertices[0]->pos;
	struct vec2 *v1 = &triangle->element.vertices[1]->pos;
	struct vec2 *v2 = &triangle->element.vertices[2]->pos;

	/* area = 0.5 * cross product of edges */
	vec2_sub(v1, v0, &v10);
	vec2_sub(v2, v0, &v20);
	triangle->area = fabsf(0.5f * (v10.x[0]*v20.x[1] - v10.x[1]*v20.x[0]));
}

/**
 * computes the gradients of functions that are 1 at one vertex and 0 at the
 * other 2; do this for all 3 vertices
 */
void triangle_compute_dof(struct triangle *restrict triangle)
{
	struct vec2 v01, v21;

	for (int i = 0; i < 3; i++) {
		/* v0 is the vertex where the function will be 1 */
		struct vec2 *v0 = &triangle->element.vertices[(i+0)%3]->pos;
		struct vec2 *v1 = &triangle->element.vertices[(i+1)%3]->pos;
		struct vec2 *v2 = &triangle->element.vertices[(i+2)%3]->pos;

		/* gram schmidt to get perp vector to opposite edge of v0 */
		vec2_sub(v0, v1, &v01);
		vec2_sub(v2, v1, &v21);
		vec2_scale(vec2_dot(&v01, &v21) / vec2_dot(&v21, &v21), &v21);
		vec2_sub(&v01, &v21, &triangle->dof_grad[(i+0)%3]);

		float s = 1.0f / vec2_dot(&v01, &triangle->dof_grad[(i+0)%3]);
		vec2_scale(s, &triangle->dof_grad[(i+0)%3]);
	}
}

/* function that is 1 at one vertex and 0 at other 2 */
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

	/* mn: triangle vertex loops */
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

			/* rs: xy loops */
			for (int r = 0; r < DIM; r++) {
				for (int s = 0; s < DIM; s++) {
					float entry = 0;
					if (r == s) {
						/* (D_i u_j) (D^i u^j) */
						entry += vec2_dot(&triangle->dof_grad[m], &triangle->dof_grad[n]);
					}
					/* D_i u_j D^j u_i */
					entry += triangle->dof_grad[m].x[s] * triangle->dof_grad[n].x[r];

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

	/* gravity acting on vertex */
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
