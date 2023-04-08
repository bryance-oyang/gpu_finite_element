/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @brief sends mesh to webpage via websocket_ctube
 */

#include <math.h>

#include "linear_algebra.h"
#include "visualize.h"
#include "mesh.h"
#include "ws_ctube.h"

/* linearly map x in lo, hi to Lo, Hi */
static int32_t lin_scale(float x, float lo, float hi, int32_t Lo, int32_t Hi)
{
	if (hi - lo < EPSILON) {
		return 0;
	} else if (x > hi) {
		return Hi;
	} else if (x < lo) {
		return Lo;
	} else {
		return (Hi - Lo) * (x - lo) / (hi - lo) + Lo;
	}
}

static void get_scaling(struct vis *restrict vis, struct mesh *restrict mesh)
{
	vis->min_x = FLT_MAX;
	vis->max_x = -FLT_MAX;

	vis->min_y = FLT_MAX;
	vis->max_y = -FLT_MAX;

	for (int i = 0; i < mesh->nelements; i++) {
		struct element *restrict element = mesh->elements[i];
		for (int v = 0; v < element->nvertices; v++) {
			struct vertex *restrict vert = get_vert(element, v);
			float x = vert->pos.x[0];
			float y = vert->pos.x[1];
			vis->min_x = fminf(vis->min_x, x);
			vis->min_y = fminf(vis->min_y, y);
			vis->max_x = fmaxf(vis->max_x, x);
			vis->max_y = fmaxf(vis->max_y, y);
		}
	}

	vis->mid_x = 0.5f * (vis->min_x + vis->max_x);
	vis->mid_y = 0.5f * (vis->min_y + vis->max_y);

	float dx = vis->max_x - vis->min_x;
	float dy = vis->max_y - vis->min_y;
	vis->max_x += 0.1*dx;
	vis->min_x -= 0.1*dx;
	vis->max_y += 0.1*dy;
	vis->min_y -= 0.1*dy;

	vis->slope = fmaxf((vis->max_x - vis->min_x) / IMAGE_WIDTH, (vis->max_y - vis->min_y) / IMAGE_HEIGHT);
}

static void get_xy(struct vis *restrict vis, int i, int j, float *x, float *y)
{
	*x = vis->slope * (j - IMAGE_WIDTH / 2) + vis->mid_x;
	*y = -vis->slope * (i - IMAGE_HEIGHT / 2) + vis->mid_y;
}

static bool xy_in_triangle(struct element *restrict element, struct vec *restrict c, float x, float y)
{
	struct vec2 v[3];

	for (int i = 0; i < 3; i++) {
		struct vertex *vert = get_vert(element, i);
		int id = vert->id;
		v[i] = vert->pos;

		if (vert->enabled) {
			/* position after being strained */
			v[i].x[0] += c->x[2*id + 0];
			v[i].x[1] += c->x[2*id + 1];
		}
	}

	struct vec2 w = {.x = {x, y}};
	struct vec2 w0, v01;
	float crosses[3];
	for (int i = 0; i < 3; i++) {
		vec2_sub(&v[(i+1)%3], &v[i], &v01);
		vec2_sub(&w, &v[i], &w0);
		crosses[i] = v01.x[0] * w0.x[1] - v01.x[1] * w0.x[0];
	}

	return (crosses[0] >= 0 && crosses[1] >= 0 && crosses[2] >= 0)
		|| (crosses[0] <= 0 && crosses[1] <= 0 && crosses[2] <= 0);
}

static struct element *xy_to_element(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c, float x, float y)
{
	if (vis->cached_element != NULL && xy_in_triangle(vis->cached_element, c, x, y)) {
		return vis->cached_element;
	}

	for (int i = 0; i < mesh->nelements; i++) {
		struct element *element = mesh->elements[i];
		if (xy_in_triangle(element, c, x, y)) {
			vis->cached_element = element;
			return element;
		}
	}

	return NULL;
}

int vis_init(struct vis *restrict vis, struct mesh *restrict mesh)
{
	vis->data_bytes = IMAGE_HEIGHT * IMAGE_WIDTH * sizeof(*vis->data);
	vis->data = malloc(vis->data_bytes);
	if (vis->data == NULL) {
		goto err_nodata;
	}

	vis->stresses = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof(*vis->stresses));
	if (vis->stresses == NULL) {
		goto err_nostresses;
	}

	vis->sorted_stresses = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof(*vis->sorted_stresses));
	if (vis->sorted_stresses == NULL) {
		goto err_nosorted_stresses;
	}

	vis->ctube = ws_ctube_open(9743, 2, 0, 24);
	if (vis->ctube == NULL) {
		goto err_noctube;
	}

	vis->cached_element = NULL;
	get_scaling(vis, mesh);
	return 0;

err_noctube:
	free(vis->sorted_stresses);
err_nosorted_stresses:
	free(vis->stresses);
err_nostresses:
	free(vis->data);
err_nodata:
	return -1;
}

void vis_destroy(struct vis *restrict vis)
{
	free(vis->sorted_stresses);
	free(vis->stresses);
	free(vis->data);
	ws_ctube_close(vis->ctube);
}

static int number_cmp(const void *a, const void *b)
{
	float aa = *(float *)a;
	float bb = *(float *)b;
	return (aa > bb) - (aa < bb);
}

static void fill_stresses(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c)
{
	float x, y;
	vis->nsorted_stresses = 0;

	for (int i = 0; i < IMAGE_HEIGHT; i++) {
		for (int j = 0; j < IMAGE_WIDTH; j++) {
			get_xy(vis, i, j, &x, &y);
			struct element *element = xy_to_element(vis, mesh, c, x, y);
			if (element == NULL) {
				vis->stresses[i*IMAGE_WIDTH + j] = NAN;
			} else {
				float stress = element->vtable->scalar_stress(c, element, x, y);
				vis->stresses[i*IMAGE_WIDTH + j] = stress;

				vis->sorted_stresses[vis->nsorted_stresses] = stress;
				vis->nsorted_stresses++;
			}
		}
	}
}

void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c)
{
	fill_stresses(vis, mesh, c);

	/* sort and get percentile of stresses by index */
	qsort(vis->sorted_stresses, vis->nsorted_stresses, sizeof(*vis->sorted_stresses), number_cmp);
	float min_stress = vis->sorted_stresses[(int)(0.01 * (vis->nsorted_stresses - 1))];
	float max_stress = vis->sorted_stresses[(int)(0.99 * (vis->nsorted_stresses - 1))];

	for (int i = 0; i < IMAGE_HEIGHT; i++) {
		for (int j = 0; j < IMAGE_WIDTH; j++) {
			float stress = vis->stresses[i*IMAGE_WIDTH + j];
			if (isnanf(stress)) {
				vis->data[i*IMAGE_WIDTH + j] = -1;
			} else {
				vis->data[i*IMAGE_WIDTH + j] = lin_scale(stress, min_stress, max_stress, 236, 0);
			}
		}
	}
}

void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(vis->ctube, vis->data, vis->data_bytes);
}
