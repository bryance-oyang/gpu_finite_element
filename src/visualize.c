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

#include "visualize.h"
#include "mesh.h"
#include "ws_ctube.h"

struct bounding_box {
	/* [x/y][min/max] */
	double xy[2][2];
	int ij[2][2];
};

/* linearly map x in lo, hi to Lo, Hi */
static double lin_scale(double x, double lo, double hi, double Lo, double Hi)
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
			double x = vert->pos.x[0];
			double y = vert->pos.x[1];
			vis->min_x = fmin(vis->min_x, x);
			vis->min_y = fmin(vis->min_y, y);
			vis->max_x = fmax(vis->max_x, x);
			vis->max_y = fmax(vis->max_y, y);
		}
	}

	vis->mid_x = 0.5 * (vis->min_x + vis->max_x);
	vis->mid_y = 0.5 * (vis->min_y + vis->max_y);

	double dx = vis->max_x - vis->min_x;
	double dy = vis->max_y - vis->min_y;
	vis->max_x += 0.1*dx;
	vis->min_x -= 0.1*dx;
	vis->max_y += 0.1*dy;
	vis->min_y -= 0.1*dy;

	vis->slope = fmax((vis->max_x - vis->min_x) / IMAGE_WIDTH, (vis->max_y - vis->min_y) / IMAGE_HEIGHT);
}

static void get_xy(struct vis *restrict vis, int i, int j, double *x, double *y)
{
	*x = vis->slope * (j - IMAGE_WIDTH / 2) + vis->mid_x;
	*y = -vis->slope * (i - IMAGE_HEIGHT / 2) + vis->mid_y;
}

static void get_ij(struct vis *restrict vis, double x, double y, int *i, int *j)
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
			double coord = vert->pos.x[k];
			if (vert->enabled) {
				coord += c->x[2*id + k];
			}
			bb.xy[k][0] = fmin(bb.xy[k][0], coord);
			bb.xy[k][1] = fmax(bb.xy[k][1], coord);
		}
	}

	get_ij(vis, bb.xy[0][0], bb.xy[1][1], &bb.ij[0][0], &bb.ij[1][0]);
	get_ij(vis, bb.xy[0][1], bb.xy[1][0], &bb.ij[0][1], &bb.ij[1][1]);

	return bb;
}

static inline void cathedral(double value, uint16_t *rgb)
{
	rgb[0] = 194*value;
	rgb[1] = 36 + 219*value;
	rgb[2] = 18 + 237*value;
}

static inline void cathedral_discrete(double value, uint16_t *rgb)
{
	value = ((int)(value * 16)) / 16.0;
	rgb[0] = 194*value;
	rgb[1] = 36 + 219*value;
	rgb[2] = 18 + 237*value;
}

static inline void rainbow(double value, uint16_t *rgb)
{
	double H = lin_scale(value, 0, 1, 248, 0);
	double S = 1;
	double L = -0.28 * SQR(SQR(1 - value)) + 0.5;
	double C = (1 - fabs(2*L - 1)) * S;
	double m = L - C / 2;

	double Hp = H / 60;
	double Hp_mod2 = Hp - 2 * ((int)Hp / 2);
	double X = C * (1 - fabs(Hp_mod2 - 1));

	double R, G, B;
	if (0 <= Hp  && Hp < 1) {
		R = C;
		G = X;
		B = 0;
	} else if (1 <= Hp  && Hp < 2) {
		R = X;
		G = C;
		B = 0;
	} else if (2 <= Hp && Hp  < 3) {
		R = 0;
		G = C;
		B = X;
	} else if (3 <= Hp && Hp  < 4) {
		R = 0;
		G = X;
		B = C;
	} else if (4 <= Hp && Hp  < 5) {
		R = X;
		G = 0;
		B = C;
	} else if (5 <= Hp && Hp  < 6) {
		R = C;
		G = 0;
		B = X;
	}

	rgb[0] = 255*(R + m);
	rgb[1] = 255*(G + m);
	rgb[2] = 255*(B + m);
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
	double aa = *(double *)a;
	double bb = *(double *)b;
	return (aa > bb) - (aa < bb);
}

static void fill_stresses(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c)
{
	double xy[2];

	for (int i = 0; i < IMAGE_HEIGHT * IMAGE_WIDTH; i++) {
		vis->stresses[i] = NAN;
	}

	for (int e = 0; e < mesh->nelements; e++) {
		struct element *element = mesh->elements[e];
		struct bounding_box bb = triangle_bounding_box(vis, element, c);
		for (int i = bb.ij[0][0]; i <= bb.ij[0][1]; i++) {
			for (int j = bb.ij[1][0]; j <= bb.ij[1][1]; j++) {
				get_xy(vis, i, j, &xy[0], &xy[1]);
				double stress = element->vtable->scalar_stress(c, element, xy);
				if (isnan(vis->stresses[i*IMAGE_WIDTH + j])) {
					vis->stresses[i*IMAGE_WIDTH + j] = stress;
				}
			}
		}
	}

	vis->nsorted_stresses = 0;
	for (int i = 0; i < IMAGE_HEIGHT * IMAGE_WIDTH; i++) {
		double stress = vis->stresses[i];
		if (!isnan(stress)) {
			vis->sorted_stresses[vis->nsorted_stresses] = stress;
			vis->nsorted_stresses++;
		}
	}
	qsort(vis->sorted_stresses, vis->nsorted_stresses, sizeof(*vis->sorted_stresses), number_cmp);
}

void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c)
{
	fill_stresses(vis, mesh, c);

	/* get percentile of stresses by index */
	double min_stress = vis->sorted_stresses[(int)(0.01 * (vis->nsorted_stresses - 1))];
	double max_stress = vis->sorted_stresses[(int)(0.99 * (vis->nsorted_stresses - 1))];

	printf("99%% percentile stress: %lf\n", max_stress);

	for (int i = 0; i < IMAGE_HEIGHT; i++) {
		for (int j = 0; j < IMAGE_WIDTH; j++) {
			double stress = vis->stresses[i*IMAGE_WIDTH + j];
			if (isnan(stress)) {
				for (int k = 0; k < 3; k++) {
					vis->data[(i*IMAGE_WIDTH + j)*3 + k] = 0;
				}
			} else {
				double value = lin_scale(stress, min_stress, max_stress, 0, 1);
				rainbow(value, &vis->data[(i*IMAGE_WIDTH + j)*3]);
			}
		}
	}

	double min_y = FLT_MAX;

	/* triangles x0, y0, x1, y1, x2, y2 */
	for (int e = 0; e < mesh->nelements; e++) {
		for (int v = 0; v < 3; v++) {
			struct vertex *restrict vert = get_vert(mesh->elements[e], v);
			double x, y;
			int i, j;

			if (vert->enabled) {
				int id = vert->id;
				x = vert->pos.x[0] + c->x[DIM*id + 0];
				y = vert->pos.x[1] + c->x[DIM*id + 1];
			} else {
				x = vert->pos.x[0];
				y = vert->pos.x[1];
			}
			min_y = fmin(min_y, y);
			get_ij(vis, x, y, &i, &j);
			vis->data[3*IMAGE_HEIGHT*IMAGE_WIDTH + (e*3 + v)*DIM + 0] = j;
			vis->data[3*IMAGE_HEIGHT*IMAGE_WIDTH + (e*3 + v)*DIM + 1] = i;
		}
	}

	printf("min_y: %lf\n", min_y);
}

void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(vis->ctube, vis->data, vis->data_bytes);
}
