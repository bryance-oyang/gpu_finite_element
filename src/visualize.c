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

int vis_init(struct vis *restrict vis, struct mesh *restrict mesh)
{
	vis->data_bytes = (3*DIM + 1) * mesh->nelements * sizeof(*vis->data);
	vis->data = malloc(vis->data_bytes);
	if (vis->data == NULL) {
		goto err_nodata;
	}

	vis->ctube = ws_ctube_open(9743, 2, 0, 24);
	if (vis->ctube == NULL) {
		goto err_noctube;
	}
	return 0;

err_noctube:
	free(vis->data);
err_nodata:
	return -1;
}

void vis_destroy(struct vis *restrict vis)
{
	free(vis->data);
	ws_ctube_close(vis->ctube);
}

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

static int number_cmp(const void *a, const void *b)
{
	float aa = *(float *)a;
	float bb = *(float *)b;
	return (aa > bb) - (aa < bb);
}

void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c)
{
	float min_xy = FLT_MAX;
	float max_xy = -FLT_MAX;

	/* determine min/max of coordinates */
	for (int i = 0; i < mesh->nelements; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = mesh->elements[i]->vertices[j];
			int idx = v->id;
			for (int k = 0; k < DIM; k++) {
				float coord = v->pos.x[k] + c->x[DIM*idx + k];
				if (coord < min_xy) {
					min_xy = coord;
				}
				if (coord > max_xy) {
					max_xy = coord;
				}
			}
		}
	}

	float *stresses = malloc(mesh->nelements * sizeof(*stresses));
	for (int i = 0; i < mesh->nelements; i++) {
		stresses[i] = mesh->elements[i]->scalar_stress;
	}
	/* sort and get percentile of stresses by index */
	qsort(stresses, mesh->nelements, sizeof(*stresses), number_cmp);
	float min_stress = stresses[(int)(0.01 * (mesh->nelements - 1))];
	float max_stress = stresses[(int)(0.99 * (mesh->nelements - 1))];

	/* store coordinates and stresses of all triangles in the format
	 * x0,y0,x1,y1,x2,y2,stress */
	for (int i = 0; i < mesh->nelements; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = mesh->elements[i]->vertices[j];
			int idx = v->id;
			for (int k = 0; k < DIM; k++) {
				int32_t Lo, Hi;
				if (k == 0) {
					/* x */
					Lo = 0;
					Hi = IMAGE_WIDTH;
				} else {
					/* y */
					Lo = IMAGE_HEIGHT;
					Hi = 0;
				}

				float coord = v->pos.x[k] + c->x[DIM*idx + k];
				vis->data[i*(3*DIM + 1) + j*DIM + k] = lin_scale(coord, min_xy - 0.1*fabsf(min_xy), max_xy + 0.1*fabsf(max_xy), Lo, Hi);
			}
		}

		float stress = mesh->elements[i]->scalar_stress;
		vis->data[(i+1)*(3*DIM + 1) - 1] = lin_scale(stress, min_stress, max_stress, 236, 0);
	}
}

void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(vis->ctube, vis->data, vis->data_bytes);
}
