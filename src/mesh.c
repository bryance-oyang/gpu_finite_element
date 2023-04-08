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

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mesh.h"
#include "container_of.h"

#define HBUF_LEN 256
char hash_buf[HBUF_LEN];

static struct element_vtable triangle_vtable = {
	.stiffness_add = triangle_stiffness_add,
	.forces_add = triangle_forces_add,
	.scalar_stress = triangle_scalar_stress,
};

/**
 * each degree of freedom is a0 r^2 + a1 r s + a2 s^2 + a3 r + a4 s + a5
 * index ordering: vertex, xy, vecind, partial, coeff#
 */

#define TRIANGLE2_NVERTEX 6
#define TRIANGLE2_NXY 2
#define TRIANGLE2_NDERIV 2
#define TRIANGLE2_NCOEFF 6
#define TRIANGLE2_NDCOEFF 3

/**
 * f_v = function 1 on vertex v and 0 on others
 * f_v = a_0 r^2 + a_1 rs + a_2 s^2 + a_3 r + a_4 s + a_5
 */
static struct canon_triangle2 {
	int is_computed;
	/* r^2 + rs + s^2 + r + s + 1 */
	float a[TRIANGLE2_NVERTEX][TRIANGLE2_NCOEFF];
	/* A_0 r + A_1 s + A_2 */
	float Da[TRIANGLE2_NVERTEX][TRIANGLE2_NDERIV][TRIANGLE2_NDCOEFF];
	/* integrals of r^2, rs, s^2, r, s, 1 */
	float I0[TRIANGLE2_NCOEFF];
	/* integral of D_dr0 f_v0 D_dr1 f_v1 */
	float I1[TRIANGLE2_NVERTEX][TRIANGLE2_NVERTEX][TRIANGLE2_NDERIV][TRIANGLE2_NDERIV];
	/* u^j_m  */
	float I3[TRIANGLE2_NVERTEX];
} canon_triangle2;

static struct element_vtable triangle2_vtable = {
	.stiffness_add = triangle2_stiffness_add,
	.forces_add = triangle2_forces_add,
	.scalar_stress = triangle2_scalar_stress,
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

int mesh_add_vertex(struct mesh *restrict mesh, float x, float y, bool enabled)
{
	if (mesh->nvertices == mesh->vertices_size) {
		int size = mesh->vertices_size * 2;
		void *tmp = realloc(mesh->vertices, size * sizeof(*mesh->vertices));
		if (tmp == NULL) {
			return -1;
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

	return idx;
}

struct vertex *get_vert(struct element *restrict element, int vidx)
{
	return &element->mesh->vertices[element->vertices[vidx]];
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
	int *v0, int *v1, struct ht_node **pprev)
{
	/* order v0, v1 by index */
	if (v0 > v1) {
		int tmp = *v0;
		*v0 = *v1;
		*v1 = tmp;
	}

	/* key for hash table */
	if (snprintf(hash_buf, HBUF_LEN, "%d_%d", *v0, *v1) >= HBUF_LEN) {
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

static struct edge *add_get_edge(struct mesh *restrict mesh, int v0, int v1)
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

static int mesh_add_element(struct mesh *restrict mesh, struct element *restrict element)
{
	int nelements = mesh->nelements + 1;
	void *elements = realloc(mesh->elements, nelements * sizeof(*mesh->elements));
	if (elements == NULL) {
		return -1;
	}
	mesh->elements = elements;
	mesh->nelements = nelements;
	mesh->elements[nelements - 1] = element;
	return 0;
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

static void triangle_compute_area(struct triangle *restrict triangle)
{
	struct vec2 v10, v20;

	struct vec2 *v0 = &get_vert(&triangle->element, 0)->pos;
	struct vec2 *v1 = &get_vert(&triangle->element, 1)->pos;
	struct vec2 *v2 = &get_vert(&triangle->element, 2)->pos;

	/* area = 0.5 * cross product of edges */
	vec2_sub(v1, v0, &v10);
	vec2_sub(v2, v0, &v20);
	triangle->area = fabsf(0.5f * (v10.x[0]*v20.x[1] - v10.x[1]*v20.x[0]));
}

/**
 * computes the gradients of functions that are 1 at one vertex and 0 at the
 * other 2; do this for all 3 vertices
 */
static void triangle_compute_dof(struct triangle *restrict triangle)
{
	struct vec2 v01, v21;

	for (int i = 0; i < 3; i++) {
		/* v0 is the vertex where the function will be 1 */
		struct vec2 *v0 = &get_vert(&triangle->element, (i+0)%3)->pos;
		struct vec2 *v1 = &get_vert(&triangle->element, (i+1)%3)->pos;
		struct vec2 *v2 = &get_vert(&triangle->element, (i+2)%3)->pos;

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
static float triangle_dof(struct triangle *restrict triangle, int vertex, float x, float y)
{
	struct vec2 vxy;
	vxy.x[0] = x;
	vxy.x[1] = y;
	return 1 + vec2_dot(&triangle->dof_grad[vertex], &vxy) - vec2_dot(&triangle->dof_grad[vertex], &get_vert(&triangle->element, vertex)->pos);
}

struct triangle *mesh_add_triangle(struct mesh *restrict mesh, int v0,
	int v1, int v2, float density, float elasticity)
{
	struct triangle *triangle = malloc(sizeof(*triangle));
	if (triangle == NULL || mesh_add_element(mesh, &triangle->element) != 0) {
		return NULL;
	}

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

/*
 * integral of (E/2) * (D_i u_j D_i u_j + D_i u_i + D_j u_j)
 */
void triangle_stiffness_add(struct sparse *restrict A, struct element *restrict element)
{
	struct triangle *triangle = container_of(element, struct triangle, element);

	/* mn: triangle vertex loops */
	for (int m = 0; m < 3; m++) {
		struct vertex *vm = get_vert(element, m);
		if (!vm->enabled) {
			continue;
		}

		for (int n = 0; n < 3; n++) {
			struct vertex *vn = get_vert(element, n);
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

void triangle_forces_add(struct vec *restrict b, struct element *restrict element)
{
	struct triangle *triangle = container_of(element, struct triangle, element);

	/* gravity acting on vertex */
	float third_weight = triangle->area * element->density / 3;
	for (int n = 0; n < 3; n++) {
		struct vertex *v = get_vert(element, n);
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
		struct vertex *vert = get_vert(element, m);
		if (!vert->enabled) {
			continue;
		}

		int i = vert->id;

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

	element->scalar_stress = sqrtf(1.5 * (SQR(sxx) + 2*SQR(sxy) + SQR(syy)));
}

/**
 * integrates (A[0] r + A[1] s + A[2]) * (D[0] r + D[1] s + D[2])
 */
static float canon_triangle2_integral_grad(float *A, float *D)
{
	return
	canon_triangle2.I0[0] * (A[0]*D[0] + A[1]*D[1]) +
	canon_triangle2.I0[5] * (A[2]*D[2]) +
	canon_triangle2.I0[1] * (A[0]*D[1] + A[1]*D[0]) +
	canon_triangle2.I0[3] * ((A[0] + A[1])*D[2] + (D[0] + D[1])*A[2]);
}

static void canon_triangle2_acoeff()
{
	int dim = 6;
	float M[36] = {
	/*	r^2	rs	s^2	r	s	1 */
		0,	0,	0,	0,	0,	1,	/* vertex 0 */
		1,	0,	0,	1,	0,	1,	/* vertex 1 */
		0,	0,	1,	0,	1,	1,	/* vertex 2 */
		0.25,	0.25,	0.25,	0.5,	0.5,	1,	/* vertex 3 */
		0,	0,	0.25,	0,	0.5,	1,	/* vertex 4 */
		0.25,	0,	0,	0.5,	0,	1	/* vertex 5 */
	};

	float inv_M[36];
	get_inverse(M, dim, inv_M);
	for (int i = 0; i < dim; i++) {
	for (int j = 0; j < dim; j++) {
		canon_triangle2.a[j][i] = inv_M[i*dim + j];
	}}
}

static void canon_triangle2_Dacoeff()
{
	for (int v = 0; v < TRIANGLE2_NVERTEX; v++) {
		canon_triangle2.Da[v][0][0] = 2 * canon_triangle2.a[v][0];
		canon_triangle2.Da[v][0][1] = canon_triangle2.a[v][1];
		canon_triangle2.Da[v][0][2] = canon_triangle2.a[v][3];

		canon_triangle2.Da[v][1][0] = canon_triangle2.a[v][1];
		canon_triangle2.Da[v][1][1] = 2 * canon_triangle2.a[v][2];
		canon_triangle2.Da[v][1][2] = canon_triangle2.a[v][4];
	}
}

static void canon_triangle2_I1()
{
	for (int v0 = 0; v0 < TRIANGLE2_NVERTEX; v0++) {
	for (int v1 = 0; v1 < TRIANGLE2_NVERTEX; v1++) {
	for (int dr0 = 0; dr0 < TRIANGLE2_NDERIV; dr0++) {
	for (int dr1 = 0; dr1 < TRIANGLE2_NDERIV; dr1++) {
		canon_triangle2.I1[v0][v1][dr0][dr1] = canon_triangle2_integral_grad(
			canon_triangle2.Da[v0][dr0],
			canon_triangle2.Da[v1][dr1]
		);
	}}}}
}

static void canon_triangle2_I3()
{
	for (int v = 0; v < TRIANGLE2_NVERTEX; v++) {
	for (int c = 0; c < TRIANGLE2_NCOEFF; c++) {
		canon_triangle2.I3[v] += canon_triangle2.a[v][c] * canon_triangle2.I0[c];
	}}
}

static void canon_triangle2_compute_all()
{
	canon_triangle2.I0[0] = 1.0f/12.0f; /* r^2 */
	canon_triangle2.I0[1] = 1.0f/24.0f; /* rs */
	canon_triangle2.I0[2] = 1.0f/12.0f; /* s^2 */
	canon_triangle2.I0[3] = 1.0f/6.0f; /* r */
	canon_triangle2.I0[4] = 1.0f/6.0f; /* s */
	canon_triangle2.I0[5] = 0.5f; /* 1 */

	canon_triangle2_acoeff();
	canon_triangle2_Dacoeff();
	canon_triangle2_I1();
	canon_triangle2_I3();

	canon_triangle2.is_computed = 1;
}

/** midpoints for edges */
static void triangle2_add_edge_midpoints(struct triangle2 *restrict triangle2, int v0, int v1, int v2)
{
	int vidx[3] = {v0, v1, v2};
	int midvidx[3] = {5, 3, 4}; // indices of midpoint in order of v0v1, v1v2, v2v0
	struct vec2 midpoint;
	int midv;

	struct element *restrict element = &triangle2->element;
	struct mesh *restrict mesh = element->mesh;

	for (int i = 0; i < 3; i++) {
		/* is_boundary true means edge is not shared yet, hence midpoint vertex not created yet */
		if (element->edges[i]->is_boundary) {
			/* midpoint of v0, v1, cyclical */
			vec2_midpoint(&mesh->vertices[vidx[(i+0)%3]].pos, &mesh->vertices[vidx[(i+1)%3]].pos, &midpoint);
			if ((midv = mesh_add_vertex(mesh, midpoint.x[0], midpoint.x[1], true)) < 0) {
				raise(SIGSEGV);
			}
			element->edges[i]->vertices[2] = midv;
			element->vertices[midvidx[i]] = midv;
		} else {
			element->vertices[midvidx[i]] = element->edges[i]->vertices[2];
		}
	}
}

static void triangle2_compute_geometry(struct triangle2 *restrict triangle2)
{
	float J[2][2];
	struct vec2 v01, v02;
	struct vec2 *v0 = &get_vert(&triangle2->element, 0)->pos;
	struct vec2 *v1 = &get_vert(&triangle2->element, 1)->pos;
	struct vec2 *v2 = &get_vert(&triangle2->element, 2)->pos;

	vec2_sub(v1, v0, &v01);
	vec2_sub(v2, v0, &v02);
	for (int i = 0; i < 2; i++) {
		J[i][0] = v01.x[i];
		J[i][1] = v02.x[i];
	}

	float det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	triangle2->inv_J[0][0] = J[1][1] / det_J;
	triangle2->inv_J[0][1] = -J[0][1] / det_J;
	triangle2->inv_J[1][0] = -J[1][0] / det_J;
	triangle2->inv_J[1][1] = J[0][0] / det_J;
	triangle2->jacob = fabsf(det_J);
}

struct triangle2 *mesh_add_triangle2(struct mesh *restrict mesh, int v0,
	int v1, int v2, float density, float elasticity)
{
	struct triangle2 *triangle2 = malloc(sizeof(*triangle2));
	if (triangle2 == NULL || mesh_add_element(mesh, &triangle2->element) != 0) {
		return NULL;
	}
	struct element *element = &triangle2->element;

	element->vtable = &triangle2_vtable;
	element->mesh = mesh;

	element->nvertices = 6; /* includes midpoints of edges */
	element->vertices[0] = v0;
	element->vertices[1] = v1;
	element->vertices[2] = v2;

	element->density = density;
	element->elasticity = elasticity;

	element->nedges = 3;
	element->edges[0] = add_get_edge(mesh, v0, v1);
	element->edges[1] = add_get_edge(mesh, v1, v2);
	element->edges[2] = add_get_edge(mesh, v2, v0);

	triangle2_add_edge_midpoints(triangle2, v0, v1, v2);
	triangle2_compute_geometry(triangle2);

	return triangle2;
}

void triangle2_stiffness_add(struct sparse *restrict A, struct element *restrict element)
{
	if (!canon_triangle2.is_computed) {
		canon_triangle2_compute_all();
	}

	struct triangle2 *triangle2 = container_of(element, struct triangle2, element);

	for (int v0 = 0; v0 < TRIANGLE2_NVERTEX; v0++) {
	for (int v1 = 0; v1 < TRIANGLE2_NVERTEX; v1++) {
		struct vertex *vert0 = get_vert(element, v0);
		if (!vert0->enabled) {
			continue;
		}
		int id0 = vert0->id;

		struct vertex *vert1 = get_vert(element, v1);
		if (!vert1->enabled) {
			continue;
		}
		int id1 = vert1->id;

		for (int xy0 = 0; xy0 < TRIANGLE2_NXY; xy0++) {
		for (int xy1 = 0; xy1 < TRIANGLE2_NXY; xy1++) {
			float entry = 0;

			for (int dr0 = 0; dr0 < TRIANGLE2_NDERIV; dr0++) {
			for (int dr1 = 0; dr1 < TRIANGLE2_NDERIV; dr1++) {
				/* D_i u_j D^i u^j */
				if (xy0 == xy1) {
					for (int i = 0; i < 2; i++) {
						entry += canon_triangle2.I1[v0][v1][dr0][dr1]
							* triangle2->inv_J[dr0][i] * triangle2->inv_J[dr1][i];
					}
				}

				/* D_i u^j D_j u^i */
				entry += canon_triangle2.I1[v0][v1][dr0][dr1]
					* triangle2->inv_J[dr0][xy1] * triangle2->inv_J[dr1][xy0];
			}}

			entry *= triangle2->jacob * element->elasticity / 2;
			sparse_add(A, 2*id0 + xy0, 2*id1 + xy1, entry);
		}}
	}}
}

void triangle2_forces_add(struct vec *restrict b, struct element *restrict element)
{
	struct triangle2 *triangle2 = container_of(element, struct triangle2, element);

	float factor = triangle2->jacob * element->density;

	for (int v = 0; v < TRIANGLE2_NVERTEX; v++) {
		struct vertex *vert = get_vert(element, v);
		if (!vert->enabled) {
			continue;
		}
		int id = vert->id;

		b->x[2*id + 1] -= factor * canon_triangle2.I3[v];
	}
}

void triangle2_scalar_stress(struct vec *restrict c, struct element *restrict element)
{

}
