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
#include <signal.h>
#include <math.h>

#include "elements.h"

#define HBUF_LEN 256
char hash_buf[HBUF_LEN];

/** first-order triangle */
static struct element_vtable triangle_vtable = {
	.stiffness_add = triangle_stiffness_add,
	.forces_add = triangle_forces_add,
	.scalar_stress = triangle_scalar_stress,
};

/**
 * second-order triangle:
 * each degree of freedom is a0 r^2 + a1 r s + a2 s^2 + a3 r + a4 s + a5
 * index ordering: vertex, vecind, partial, coeff#
 */

#define TRIANGLE2_NVERTEX 6
#define TRIANGLE2_NDERIV 2
#define TRIANGLE2_NCOEFF 6
#define TRIANGLE2_NDCOEFF 3
#define TRIANGLE2_NXY 2

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
 * f_v = function 1 on vertex v and 0 on others
 * f_v = a_0 r^2 + a_1 rs + a_2 s^2 + a_3 r + a_4 s + a_5
 */
static struct canon_triangle2 {
	int is_computed;
	/* a_0 r^2 + a_1 rs + a_2 s^2 + a_3 r + a_4 s + a_5 */
	double a[TRIANGLE2_NVERTEX][TRIANGLE2_NCOEFF];
	/* A_0 r + A_1 s + A_2 */
	double Da[TRIANGLE2_NVERTEX][TRIANGLE2_NDERIV][TRIANGLE2_NDCOEFF];
	/* integrals of r^2, rs, s^2, r, s, 1 */
	double I0[TRIANGLE2_NCOEFF];
	/* integral of D_dr0 f_v0 D_dr1 f_v1 */
	double I1[TRIANGLE2_NVERTEX][TRIANGLE2_NVERTEX][TRIANGLE2_NDERIV][TRIANGLE2_NDERIV];
	/* integral of f_v  */
	double I2[TRIANGLE2_NVERTEX];
} canon_triangle2;

static struct element_vtable triangle2_vtable = {
	.stiffness_add = triangle2_stiffness_add,
	.forces_add = triangle2_forces_add,
	.scalar_stress = triangle2_scalar_stress,
};

/**
 * third-order triangle:
 * each degree of freedom is a0 r^3 + a1 r^2 s + a2 r s^2 + a3 s^3 + a4 r^2 + a5 r s + a6 s^2 + a7 r + a8 s + a9
 * index ordering: vertex, vecind, partial, coeff#
 */

#define TRIANGLE3_NVERTEX 10
#define TRIANGLE3_NDERIV 2
#define TRIANGLE3_NCOEFF 10
#define TRIANGLE3_NDCOEFF 6
#define TRIANGLE3_NXY 2
#define TRIANGLE3_NINTEGRAL 15

/**
 *    2
 *    |.
 *   7| .6
 *    |  .
 *   8| 9 .5
 *    |    .
 *   0|-----1
 *      3 4
 *
 * f_v = function 1 on vertex v and 0 on others
 * f_v = a0 r^3 + a1 r^2 s + a2 r s^2 + a3 s^3 + a4 r^2 + a5 r s + a6 s^2 + a7 r + a8 s + a9
 */
static struct canon_triangle3 {
	int is_computed;
	/* a0 r^3 + a1 r^2 s + a2 r s^2 + a3 s^3 + a4 r^2 + a5 r s + a6 s^2 + a7 r + a8 s + a9 */
	double a[TRIANGLE3_NVERTEX][TRIANGLE3_NCOEFF];
	/* A_0 r^2 + A_1 rs + A_2 s^2 + A_3 r + A_4 s + A_5 */
	double Da[TRIANGLE3_NVERTEX][TRIANGLE3_NDERIV][TRIANGLE3_NDCOEFF];
	/* integrals of r^3, r^2 s, ... and r^4, r^3 s, ... */
	double I0[TRIANGLE3_NINTEGRAL];
	/* integral of D_dr0 f_v0 D_dr1 f_v1 */
	double I1[TRIANGLE3_NVERTEX][TRIANGLE3_NVERTEX][TRIANGLE3_NDERIV][TRIANGLE3_NDERIV];
	/* integral of f_v  */
	double I2[TRIANGLE3_NVERTEX];
} canon_triangle3;

static struct element_vtable triangle3_vtable = {
	.stiffness_add = triangle3_stiffness_add,
	.forces_add = triangle3_forces_add,
	.scalar_stress = triangle3_scalar_stress,
};

/* helper function for computing stresses in triangle by mapping xy to rs coords */
static void canon_triangle_coord(struct vec *restrict c,
	struct element *restrict element, double *restrict xy, double *restrict rs)
{
	struct vec2 vxy = {.x = {xy[0], xy[1]}};
	struct vec2 v[3];
	struct vec2 v01, v02, xy0;
	double J[4];
	double inv_J[4];

	for (int i = 0; i < 3; i++) {
		struct vertex *restrict vert = get_vert(element, i);
		v[i] = vert->pos;
		if (vert->enabled) {
			int id = vert->id;
			for (int k = 0; k < 2; k++) {
				v[i].x[k] += c->x[2*id + k];
			}
		}
	}

	vec2_sub(&v[1], &v[0], &v01);
	vec2_sub(&v[2], &v[0], &v02);
	vec2_sub(&vxy, &v[0], &xy0);

	for (int i = 0; i < 2; i++) {
		J[2*i + 0] = v01.x[i];
		J[2*i + 1] = v02.x[i];
	}
	inverse_matrix2(J, inv_J);

	rs[0] = 0;
	rs[1] = 0;
	for (int j = 0; j < 2; j++) {
		rs[0] += inv_J[2*0 + j] * xy0.x[j];
		rs[1] += inv_J[2*1 + j] * xy0.x[j];
	}
}

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

int mesh_add_vertex(struct mesh *restrict mesh, double x, double y, bool enabled)
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
	if (*v0 > *v1) {
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

static void triangle_compute_area(struct triangle *restrict triangle)
{
	struct vec2 v10, v20;

	struct vec2 *v0 = &get_vert(&triangle->element, 0)->pos;
	struct vec2 *v1 = &get_vert(&triangle->element, 1)->pos;
	struct vec2 *v2 = &get_vert(&triangle->element, 2)->pos;

	/* area = 0.5 * cross product of edges */
	vec2_sub(v1, v0, &v10);
	vec2_sub(v2, v0, &v20);
	triangle->area = fabs(0.5 * (v10.x[0]*v20.x[1] - v10.x[1]*v20.x[0]));
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

		double s = 1.0 / vec2_dot(&v01, &triangle->dof_grad[(i+0)%3]);
		vec2_scale(s, &triangle->dof_grad[(i+0)%3]);
	}
}

struct triangle *mesh_add_triangle(struct mesh *restrict mesh, int v0,
	int v1, int v2, double density, double elasticity)
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
					double entry = 0;
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
	double third_weight = triangle->area * element->density / 3;
	for (int n = 0; n < 3; n++) {
		struct vertex *v = get_vert(element, n);
		if (v->enabled) {
			int i = v->id;
			b->x[DIM*i + 1] += -third_weight;
		}
	}
}

double triangle_scalar_stress(struct vec *restrict c, struct element *restrict element, double *x)
{
	double rs[2];
	canon_triangle_coord(c, element, x, rs);
	if (!(rs[0] >= 0 && rs[1] >= 0 && rs[0] + rs[1] <= 1)) {
		return NAN;
	}

	struct triangle *triangle = container_of(element, struct triangle, element);

	double sxx = 0;
	double sxy = 0;
	double syy = 0;

	for (int m = 0; m < 3; m++) {
		struct vertex *vert = get_vert(element, m);
		if (!vert->enabled) {
			continue;
		}

		int i = vert->id;

		sxx += c->x[2*i] * triangle->dof_grad[m].x[0];
		syy += c->x[2*i + 1] * triangle->dof_grad[m].x[1];
		sxy += 0.5 * (c->x[2*i] * triangle->dof_grad[m].x[1] + c->x[2*i + 1] * triangle->dof_grad[m].x[0]);
	}

	sxx *= element->elasticity;
	sxy *= element->elasticity;
	syy *= element->elasticity;

	double pressure = -0.5 * (sxx + syy);
	sxx += pressure;
	syy += pressure;

	triangle->scalar_stress = sqrt(1.5 * (SQR(sxx) + 2*SQR(sxy) + SQR(syy)));

	return triangle->scalar_stress;
}

/**
 * integrates (A[0] r + A[1] s + A[2]) * (D[0] r + D[1] s + D[2])
 */
static double canon_triangle2_integral_grad(double *A, double *D)
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
	double M[36] = {
	/*	r^2	rs	s^2	r	s	1 */
		0,	0,	0,	0,	0,	1,	/* vertex 0 */
		1,	0,	0,	1,	0,	1,	/* vertex 1 */
		0,	0,	1,	0,	1,	1,	/* vertex 2 */
		0.25,	0.25,	0.25,	0.5,	0.5,	1,	/* vertex 3 */
		0,	0,	0.25,	0,	0.5,	1,	/* vertex 4 */
		0.25,	0,	0,	0.5,	0,	1	/* vertex 5 */
	};

	double inv_M[36];
	inverse_matrix(M, dim, inv_M);
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

static void canon_triangle2_I2()
{
	for (int v = 0; v < TRIANGLE2_NVERTEX; v++) {
	for (int c = 0; c < TRIANGLE2_NCOEFF; c++) {
		canon_triangle2.I2[v] += canon_triangle2.a[v][c] * canon_triangle2.I0[c];
	}}
}

static void canon_triangle2_compute_all()
{
	canon_triangle2.I0[0] = 1.0/12.0; /* r^2 */
	canon_triangle2.I0[1] = 1.0/24.0; /* rs */
	canon_triangle2.I0[2] = 1.0/12.0; /* s^2 */
	canon_triangle2.I0[3] = 1.0/6.0; /* r */
	canon_triangle2.I0[4] = 1.0/6.0; /* s */
	canon_triangle2.I0[5] = 0.5; /* 1 */

	canon_triangle2_acoeff();
	canon_triangle2_Dacoeff();
	canon_triangle2_I1();
	canon_triangle2_I2();

	canon_triangle2.is_computed = 1;
}

/** midpoints for edges */
static void triangle2_add_edge_midpoints(struct triangle2 *restrict triangle2, int v0, int v1, int v2)
{
	int vidx[3] = {v0, v1, v2};
	int midvidx[3] = {5, 3, 4}; /* indices of midpoint in order of v0v1, v1v2, v2v0 */
	struct vec2 midpoint;
	int midv;

	struct element *restrict element = &triangle2->element;
	struct mesh *restrict mesh = element->mesh;

	for (int i = 0; i < 3; i++) {
		/* is_boundary true means edge is not shared yet, hence midpoint vertex not created yet */
		if (element->edges[i]->is_boundary) {
			/* midpoint of v0, v1, cyclical */
			struct vertex *v0 = &mesh->vertices[vidx[(i+0)%3]];
			struct vertex *v1 = &mesh->vertices[vidx[(i+1)%3]];
			vec2_midpoint(&v0->pos, &v1->pos, &midpoint);

			bool enabled = v0->enabled || v1->enabled;
			if ((midv = mesh_add_vertex(mesh, midpoint.x[0], midpoint.x[1], enabled)) < 0) {
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
	double J[4];
	struct vec2 v01, v02;
	struct vec2 *v0 = &get_vert(&triangle2->element, 0)->pos;
	struct vec2 *v1 = &get_vert(&triangle2->element, 1)->pos;
	struct vec2 *v2 = &get_vert(&triangle2->element, 2)->pos;

	vec2_sub(v1, v0, &v01);
	vec2_sub(v2, v0, &v02);
	for (int i = 0; i < 2; i++) {
		J[2*i + 0] = v01.x[i];
		J[2*i + 1] = v02.x[i];
	}

	triangle2->jacob = fabs(matrix_det2(J));
	inverse_matrix2(J, triangle2->inv_J);
}

struct triangle2 *mesh_add_triangle2(struct mesh *restrict mesh, int v0,
	int v1, int v2, double density, double elasticity)
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
			double entry = 0;

			for (int dr0 = 0; dr0 < TRIANGLE2_NDERIV; dr0++) {
			for (int dr1 = 0; dr1 < TRIANGLE2_NDERIV; dr1++) {
				/* D_i u_j D^i u^j */
				if (xy0 == xy1) {
					for (int i = 0; i < 2; i++) {
						entry += canon_triangle2.I1[v0][v1][dr0][dr1]
							* triangle2->inv_J[2*dr0 + i] * triangle2->inv_J[2*dr1 + i];
					}
				}

				/* D_i u^j D_j u^i */
				entry += canon_triangle2.I1[v0][v1][dr0][dr1]
					* triangle2->inv_J[2*dr0 + xy1] * triangle2->inv_J[2*dr1 + xy0];
			}}

			entry *= triangle2->jacob * element->elasticity / 2;
			sparse_add(A, 2*id0 + xy0, 2*id1 + xy1, entry);
		}}
	}}
}

void triangle2_forces_add(struct vec *restrict b, struct element *restrict element)
{
	struct triangle2 *triangle2 = container_of(element, struct triangle2, element);

	double factor = triangle2->jacob * element->density;

	for (int v = 0; v < TRIANGLE2_NVERTEX; v++) {
		struct vertex *vert = get_vert(element, v);
		if (!vert->enabled) {
			continue;
		}
		int id = vert->id;

		b->x[2*id + 1] -= factor * canon_triangle2.I2[v];
	}
}

double triangle2_scalar_stress(struct vec *restrict c, struct element *restrict element, double *x)
{
	/* canonical coordinates and additional 1 for easier index contraction */
	double r[3] = {0, 0, 1};
	canon_triangle_coord(c, element, x, r);
	if (!(r[0] >= 0 && r[1] >= 0 && r[0] + r[1] <= 1)) {
		return NAN;
	}

	struct triangle2 *restrict triangle2 = container_of(element, struct triangle2, element);

	double s[2][2] = {{0, 0}, {0, 0}};
	for (int v = 0; v < TRIANGLE2_NVERTEX; v++) {
	struct vertex *vert = get_vert(element, v);
	if (!vert->enabled) {
		continue;
	}
	int id = vert->id;

	for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
	for (int dr = 0; dr < TRIANGLE2_NDERIV; dr++) {
	for (int d = 0; d < TRIANGLE2_NDCOEFF; d++) {
		s[i][j] += triangle2->inv_J[2*dr + i]
			* canon_triangle2.Da[v][dr][d]
			* r[d]
			* c->x[2*id + j];
	}}}}}

	/* symmetrize */
	s[0][0] *= 2;
	s[1][1] *= 2;
	s[0][1] += s[1][0];
	s[1][0] = s[0][1];

	for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
		s[i][j] *= 0.5 * element->elasticity;
	}}

	double pressure = -0.5 * (s[0][0] + s[1][1]);
	s[0][0] += pressure;
	s[1][1] += pressure;

	return sqrt(1.5 * (SQR(s[0][0]) + 2*SQR(s[0][1]) + SQR(s[1][1])));
}

/**
 * integrates (A_0 r^2 + A_1 rs + A_2 s^2 + A_3 r + A_4 s + A_5)
 * * (B_0 r^2 + B_1 rs + B_2 s^2 + B_3 r + B_4 s + B_5)
 */
static double canon_triangle3_integral_grad(double *A, double *B)
{
	return
	canon_triangle3.I0[10] * (A[0]*B[0]) +
	canon_triangle3.I0[11] * (A[0]*B[1] + A[1]*B[0]) +
	canon_triangle3.I0[12] * (A[0]*B[2] + A[1]*B[1] + A[2]*B[0]) +
	canon_triangle3.I0[13] * (A[1]*B[2] + A[2]*B[1]) +
	canon_triangle3.I0[14] * (A[2]*B[2]) +

	canon_triangle3.I0[0] * (A[0]*B[3] + A[3]*B[0]) +
	canon_triangle3.I0[1] * (A[0]*B[4] + A[4]*B[0] + A[1]*B[3] + A[3]*B[1]) +
	canon_triangle3.I0[2] * (A[3]*B[2] + A[2]*B[3] + A[1]*B[4] + A[4]*B[1]) +
	canon_triangle3.I0[3] * (A[2]*B[4] + A[4]*B[2]) +

	canon_triangle3.I0[4] * (A[3]*B[3] + A[0]*B[5] + A[5]*B[0]) +
	canon_triangle3.I0[5] * (A[1]*B[5] + A[5]*B[1] + A[3]*B[4] + A[4]*B[3]) +
	canon_triangle3.I0[6] * (A[4]*B[4] + A[2]*B[5] + A[5]*B[2]) +

	canon_triangle3.I0[7] * (A[3]*B[5] + A[5]*B[3]) +
	canon_triangle3.I0[8] * (A[4]*B[5] + A[5]*B[4]) +

	canon_triangle3.I0[9] * (A[5]*B[5]);
}

static void canon_triangle3_acoeff()
{
	int dim = 10;
	double M[100] = {
	/*	r^3	r^2 s	r s^2	s^3	r^2	rs	s^2	r	s	1 */
		0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	/* vertex 0 */
		1,	0,	0,	0,	1,	0,	0,	1,	0,	1,	/* vertex 1 */
		0,	0,	0,	1,	0,	0,	1,	0,	1,	1,	/* vertex 2 */
		1.0/27,	0,	0,	0,	1.0/9,	0,	0,	1.0/3,	0,	1,	/* vertex 3 */
		8.0/27,	0,	0,	0,	4.0/9,	0,	0,	2.0/3,	0,	1,	/* vertex 4 */
		8.0/27,	4.0/27,	2.0/27,	1.0/27,	4.0/9,	2.0/9,	1.0/9,	2.0/3,	1.0/3,	1,	/* vertex 5 */
		1.0/27,	2.0/27,	4.0/27,	8.0/27,	1.0/9,	2.0/9,	4.0/9,	1.0/3,	2.0/3,	1,	/* vertex 6 */
		0,	0,	0,	8.0/27,	0,	0,	4.0/9,	0,	2.0/3,	1,	/* vertex 7 */
		0,	0,	0,	1.0/27,	0,	0,	1.0/9,	0,	1.0/3,	1,	/* vertex 8 */
		1.0/27,	1.0/27,	1.0/27,	1.0/27,	1.0/9,	1.0/9,	1.0/9,	1.0/3,	1.0/3,	1,	/* vertex 9 */
	};

	double inv_M[100];
	inverse_matrix(M, dim, inv_M);
	for (int i = 0; i < dim; i++) {
	for (int j = 0; j < dim; j++) {
		canon_triangle3.a[j][i] = inv_M[i*dim + j];
	}}
}

static void canon_triangle3_Dacoeff()
{
	for (int v = 0; v < TRIANGLE3_NVERTEX; v++) {
		canon_triangle3.Da[v][0][0] = 3*canon_triangle3.a[v][0];
		canon_triangle3.Da[v][0][1] = 2*canon_triangle3.a[v][1];
		canon_triangle3.Da[v][0][2] = canon_triangle3.a[v][2];
		canon_triangle3.Da[v][0][3] = 2*canon_triangle3.a[v][4];
		canon_triangle3.Da[v][0][4] = canon_triangle3.a[v][5];
		canon_triangle3.Da[v][0][5] = canon_triangle3.a[v][7];

		canon_triangle3.Da[v][1][0] = canon_triangle3.a[v][1];
		canon_triangle3.Da[v][1][1] = 2*canon_triangle3.a[v][2];
		canon_triangle3.Da[v][1][2] = 3*canon_triangle3.a[v][3];
		canon_triangle3.Da[v][1][3] = canon_triangle3.a[v][5];
		canon_triangle3.Da[v][1][4] = 2*canon_triangle3.a[v][6];
		canon_triangle3.Da[v][1][5] = canon_triangle3.a[v][8];
	}
}

static void canon_triangle3_I1()
{
	for (int v0 = 0; v0 < TRIANGLE3_NVERTEX; v0++) {
	for (int v1 = 0; v1 < TRIANGLE3_NVERTEX; v1++) {
	for (int dr0 = 0; dr0 < TRIANGLE3_NDERIV; dr0++) {
	for (int dr1 = 0; dr1 < TRIANGLE3_NDERIV; dr1++) {
		canon_triangle3.I1[v0][v1][dr0][dr1] = canon_triangle3_integral_grad(
			canon_triangle3.Da[v0][dr0],
			canon_triangle3.Da[v1][dr1]
		);
	}}}}
}

static void canon_triangle3_I2()
{
	for (int v = 0; v < TRIANGLE3_NVERTEX; v++) {
	for (int c = 0; c < TRIANGLE3_NCOEFF; c++) {
		canon_triangle3.I2[v] += canon_triangle3.a[v][c] * canon_triangle3.I0[c];
	}}
}

static void canon_triangle3_compute_all()
{
	canon_triangle3.I0[0] = 1.0/20.0; /* r^3 */
	canon_triangle3.I0[1] = 1.0/60.0; /* r^2 s */
	canon_triangle3.I0[2] = 1.0/60.0; /* r s^2 */
	canon_triangle3.I0[3] = 1.0/20.0; /* s^3 */

	canon_triangle3.I0[4] = 1.0/12.0; /* r^2 */
	canon_triangle3.I0[5] = 1.0/24.0; /* r s */
	canon_triangle3.I0[6] = 1.0/12.0; /* s^2 */

	canon_triangle3.I0[7] = 1.0/6.0; /* r */
	canon_triangle3.I0[8] = 1.0/6.0; /* s */

	canon_triangle3.I0[9] = 0.5; /* 1 */

	canon_triangle3.I0[10] = 1.0/30; /* r^4 */
	canon_triangle3.I0[11] = 1.0/120; /* r^3 s */
	canon_triangle3.I0[12] = 1.0/180; /* r^2 s^2 */
	canon_triangle3.I0[13] = 1.0/120; /* r s^3 */
	canon_triangle3.I0[14] = 1.0/30; /* s^4 */

	canon_triangle3_acoeff();
	canon_triangle3_Dacoeff();
	canon_triangle3_I1();
	canon_triangle3_I2();

	canon_triangle3.is_computed = 1;
}

/** points for edges */
static void triangle3_add_edge_points(struct triangle3 *restrict triangle3, int v0, int v1, int v2)
{
	int vidx[3] = {v0, v1, v2};
	struct vec2 v01, midp;
	int midv;

	struct element *restrict element = &triangle3->element;
	struct mesh *restrict mesh = element->mesh;

	for (int i = 0; i < 3; i++) {
		/* is_boundary true means edge is not shared yet, hence midpoint vertex not created yet */
		if (element->edges[i]->is_boundary) {
			/* midpoint of v0, v1, cyclical */
			struct vertex *vert0 = &mesh->vertices[vidx[(i+0)%3]];
			struct vertex *vert1 = &mesh->vertices[vidx[(i+1)%3]];

			/* v0 + (v1 - v0) / 3 */
			vec2_sub(&vert1->pos, &vert0->pos, &v01);
			vec2_scale(1.0/3, &v01);
			vec2_add(&vert0->pos, &v01, &midp);
			bool enabled = vert0->enabled || vert1->enabled;
			if ((midv = mesh_add_vertex(mesh, midp.x[0], midp.x[1], enabled)) < 0) {
				raise(SIGSEGV);
			}
			element->edges[i]->vertices[2] = midv;
			element->vertices[2*i + 3] = midv;

			/* v0 + 2*(v1 - v0) / 3 */
			vec2_add(&midp, &v01, &midp);
			if ((midv = mesh_add_vertex(mesh, midp.x[0], midp.x[1], enabled)) < 0) {
				raise(SIGSEGV);
			}
			element->edges[i]->vertices[3] = midv;
			element->vertices[2*i + 4] = midv;
		} else {
			/* the ordering of the vertices matters: e.g. the 2*i + 3 one should be closer to v0 */
			/* v0 --- 3 --- 4 --- v1 and NOT other way */
			struct vertex *vert0 = &mesh->vertices[vidx[(i+0)%3]];
			int ev2 = element->edges[i]->vertices[2];
			int ev3 = element->edges[i]->vertices[3];
			struct vertex *evert2 = &mesh->vertices[ev2];
			struct vertex *evert3 = &mesh->vertices[ev3];
			struct vec2 v0_ev2, v0_ev3;
			vec2_sub(&evert2->pos, &vert0->pos, &v0_ev2);
			vec2_sub(&evert3->pos, &vert0->pos, &v0_ev3);

			/* determine closeness */
			if (vec2_dot(&v0_ev2, &v0_ev2) < vec2_dot(&v0_ev3, &v0_ev3)) {
				element->vertices[2*i + 3] = ev2;
				element->vertices[2*i + 4] = ev3;
			} else {
				element->vertices[2*i + 3] = ev3;
				element->vertices[2*i + 4] = ev2;
			}
		}
	}

	/* (v0 + v1 + v2) / 3 */
	struct vertex *vert0 = &mesh->vertices[v0];
	struct vertex *vert1 = &mesh->vertices[v1];
	struct vertex *vert2 = &mesh->vertices[v2];
	vec2_add(&vert0->pos, &vert1->pos, &midp);
	vec2_add(&vert2->pos, &midp, &midp);
	vec2_scale(1.0/3, &midp);
	bool enabled = vert0->enabled || vert1->enabled || vert2->enabled;
	if ((midv = mesh_add_vertex(mesh, midp.x[0], midp.x[1], enabled)) < 0) {
		raise(SIGSEGV);
	}
	element->vertices[9] = midv;
}

static void triangle3_compute_geometry(struct triangle3 *restrict triangle3)
{
	double J[4];
	struct vec2 v01, v02;
	struct vec2 *v0 = &get_vert(&triangle3->element, 0)->pos;
	struct vec2 *v1 = &get_vert(&triangle3->element, 1)->pos;
	struct vec2 *v2 = &get_vert(&triangle3->element, 2)->pos;

	vec2_sub(v1, v0, &v01);
	vec2_sub(v2, v0, &v02);
	for (int i = 0; i < 2; i++) {
		J[2*i + 0] = v01.x[i];
		J[2*i + 1] = v02.x[i];
	}

	triangle3->jacob = fabs(matrix_det2(J));
	inverse_matrix2(J, triangle3->inv_J);
}

struct triangle3 *mesh_add_triangle3(struct mesh *restrict mesh, int v0,
	int v1, int v2, double density, double elasticity)
{
	struct triangle3 *triangle3 = malloc(sizeof(*triangle3));
	if (triangle3 == NULL || mesh_add_element(mesh, &triangle3->element) != 0) {
		return NULL;
	}
	struct element *element = &triangle3->element;

	element->vtable = &triangle3_vtable;
	element->mesh = mesh;

	element->nvertices = 10; /* includes edge points and center point */
	element->vertices[0] = v0;
	element->vertices[1] = v1;
	element->vertices[2] = v2;

	element->density = density;
	element->elasticity = elasticity;

	element->nedges = 3;
	element->edges[0] = add_get_edge(mesh, v0, v1);
	element->edges[1] = add_get_edge(mesh, v1, v2);
	element->edges[2] = add_get_edge(mesh, v2, v0);

	triangle3_add_edge_points(triangle3, v0, v1, v2);
	triangle3_compute_geometry(triangle3);

	return triangle3;
}

void triangle3_stiffness_add(struct sparse *restrict A, struct element *restrict element)
{
	if (!canon_triangle3.is_computed) {
		canon_triangle3_compute_all();
	}

	struct triangle3 *triangle3 = container_of(element, struct triangle3, element);

	for (int v0 = 0; v0 < TRIANGLE3_NVERTEX; v0++) {
	for (int v1 = 0; v1 < TRIANGLE3_NVERTEX; v1++) {
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

		for (int xy0 = 0; xy0 < TRIANGLE3_NXY; xy0++) {
		for (int xy1 = 0; xy1 < TRIANGLE3_NXY; xy1++) {
			double entry = 0;

			for (int dr0 = 0; dr0 < TRIANGLE3_NDERIV; dr0++) {
			for (int dr1 = 0; dr1 < TRIANGLE3_NDERIV; dr1++) {
				/* D_i u_j D^i u^j */
				if (xy0 == xy1) {
					for (int i = 0; i < 2; i++) {
						entry += canon_triangle3.I1[v0][v1][dr0][dr1]
							* triangle3->inv_J[2*dr0 + i] * triangle3->inv_J[2*dr1 + i];
					}
				}

				/* D_i u^j D_j u^i */
				entry += canon_triangle3.I1[v0][v1][dr0][dr1]
					* triangle3->inv_J[2*dr0 + xy1] * triangle3->inv_J[2*dr1 + xy0];
			}}

			entry *= triangle3->jacob * element->elasticity / 2;
			sparse_add(A, 2*id0 + xy0, 2*id1 + xy1, entry);
		}}
	}}
}

void triangle3_forces_add(struct vec *restrict b, struct element *restrict element)
{
	struct triangle3 *triangle3 = container_of(element, struct triangle3, element);

	double factor = triangle3->jacob * element->density;

	for (int v = 0; v < TRIANGLE3_NVERTEX; v++) {
		struct vertex *vert = get_vert(element, v);
		if (!vert->enabled) {
			continue;
		}
		int id = vert->id;

		b->x[2*id + 1] -= factor * canon_triangle3.I2[v];
	}
}

double triangle3_scalar_stress(struct vec *restrict c, struct element *restrict element, double *x)
{
	/* canonical coordinates */
	double r[2] = {0, 0};
	canon_triangle_coord(c, element, x, r);
	if (!(r[0] >= 0 && r[1] >= 0 && r[0] + r[1] <= 1)) {
		return NAN;
	}

	double r_poly[6];
	r_poly[0] = SQR(r[0]); /* r^2 */
	r_poly[1] = r[0] * r[1]; /* rs */
	r_poly[2] = SQR(r[1]); /* s^2 */
	r_poly[3] = r[0]; /* r */
	r_poly[4] = r[1]; /* s */
	r_poly[5] = 1.0; /* 1 */

	struct triangle3 *restrict triangle3 = container_of(element, struct triangle3, element);

	double s[2][2] = {{0, 0}, {0, 0}};
	for (int v = 0; v < TRIANGLE3_NVERTEX; v++) {
	struct vertex *vert = get_vert(element, v);
	if (!vert->enabled) {
		continue;
	}
	int id = vert->id;

	for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
	for (int dr = 0; dr < TRIANGLE3_NDERIV; dr++) {
	for (int d = 0; d < TRIANGLE3_NDCOEFF; d++) {
		s[i][j] += triangle3->inv_J[2*dr + i]
			* canon_triangle3.Da[v][dr][d]
			* r_poly[d]
			* c->x[2*id + j];
	}}}}}

	/* symmetrize */
	s[0][0] *= 2;
	s[1][1] *= 2;
	s[0][1] += s[1][0];
	s[1][0] = s[0][1];

	for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
		s[i][j] *= 0.5 * element->elasticity;
	}}

	double pressure = -0.5 * (s[0][0] + s[1][1]);
	s[0][0] += pressure;
	s[1][1] += pressure;

	return sqrt(1.5 * (SQR(s[0][0]) + 2*SQR(s[0][1]) + SQR(s[1][1])));
}
