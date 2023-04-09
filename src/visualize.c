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

struct bounding_box {
	/* [x/y][min/max] */
	float xy[2][2];
	int ij[2][2];
};

/* linearly map x in lo, hi to Lo, Hi */
static float lin_scale(float x, float lo, float hi, float Lo, float Hi)
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

static void get_ij(struct vis *restrict vis, float x, float y, int *i, int *j)
{
	*i = -(y - vis->mid_y) / vis->slope + IMAGE_HEIGHT / 2;
	*j = (x - vis->mid_x) / vis->slope + IMAGE_WIDTH / 2;

	if (*i < 0) {
		*i = 0;
	} else if (*i >= IMAGE_HEIGHT) {
		*i = IMAGE_HEIGHT - 1;
	}
	if (*j < 0) {
		*j = 0;
	} else if (*j >= IMAGE_WIDTH) {
		*j = IMAGE_WIDTH - 1;
	}
}

static struct bounding_box triangle_bounding_box(struct vis *restrict vis, struct element *restrict element, struct vec *restrict c)
{
	struct bounding_box bb = {
		.xy = {{FLT_MAX, -FLT_MAX}, {FLT_MAX, -FLT_MAX}}
	};

	for (int v = 0; v < 3; v++) {
		struct vertex *restrict vert = get_vert(element, v);
		int id = vert->id;
		for (int k = 0; k < 2; k++) {
			float coord = vert->pos.x[k];
			if (vert->enabled) {
				coord += c->x[2*id + k];
			}
			bb.xy[k][0] = fminf(bb.xy[k][0], coord);
			bb.xy[k][1] = fmaxf(bb.xy[k][1], coord);
		}
	}

	get_ij(vis, bb.xy[0][0], bb.xy[1][1], &bb.ij[0][0], &bb.ij[1][0]);
	get_ij(vis, bb.xy[0][1], bb.xy[1][0], &bb.ij[0][1], &bb.ij[1][1]);

	return bb;
}

int vis_init(struct vis *restrict vis, struct mesh *restrict mesh)
{
	vis->data_bytes = (3*IMAGE_HEIGHT*IMAGE_WIDTH + 6*mesh->nelements) * sizeof(*vis->data);
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
	float xy[2];
	vis->nsorted_stresses = 0;

	for (int i = 0; i < IMAGE_HEIGHT * IMAGE_WIDTH; i++) {
		vis->stresses[i] = NAN;
	}

	for (int e = 0; e < mesh->nelements; e++) {
		struct element *element = mesh->elements[e];
		struct bounding_box bb = triangle_bounding_box(vis, element, c);
		for (int i = bb.ij[0][0]; i <= bb.ij[0][1]; i++) {
			for (int j = bb.ij[1][0]; j <= bb.ij[1][1]; j++) {
				get_xy(vis, i, j, &xy[0], &xy[1]);
				float stress = element->vtable->scalar_stress(c, element, xy);
				if (isnan(vis->stresses[i*IMAGE_WIDTH + j])) {
					vis->stresses[i*IMAGE_WIDTH + j] = stress;
				}

				if (!isnan(stress)) {
					vis->sorted_stresses[vis->nsorted_stresses] = stress;
					vis->nsorted_stresses++;
				}
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
			if (isnan(stress)) {
				for (int k = 0; k < 3; k++) {
					vis->data[(i*IMAGE_WIDTH + j)*3 + k] = 0;
				}
			} else {
				for (int k = 0; k < 3; k++) {
					vis->data[(i*IMAGE_WIDTH + j)*3 + k] = 0;
				}
				float value = lin_scale(stress, min_stress, max_stress, 0, 1);
				vis->data[(i*IMAGE_WIDTH + j)*3 + 0] = 194*value;
				vis->data[(i*IMAGE_WIDTH + j)*3 + 1] = 63 + 192*value;
				vis->data[(i*IMAGE_WIDTH + j)*3 + 2] = 37 + 218*value;
			}
		}
	}

	/* triangles x0, y0, x1, y1, x2, y2 */
	for (int e = 0; e < mesh->nelements; e++) {
		for (int v = 0; v < 3; v++) {
			struct vertex *restrict vert = get_vert(mesh->elements[e], v);
			float x, y;
			int i, j;

			if (vert->enabled) {
				int id = vert->id;
				x = vert->pos.x[0] + c->x[DIM*id + 0];
				y = vert->pos.x[1] + c->x[DIM*id + 1];
			} else {
				x = vert->pos.x[0];
				y = vert->pos.x[1];
			}
			get_ij(vis, x, y, &i, &j);
			vis->data[3*IMAGE_HEIGHT*IMAGE_WIDTH + (e*3 + v)*DIM + 0] = j;
			vis->data[3*IMAGE_HEIGHT*IMAGE_WIDTH + (e*3 + v)*DIM + 1] = i;
		}
	}
}

void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(vis->ctube, vis->data, vis->data_bytes);
}
