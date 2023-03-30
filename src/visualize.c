/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <math.h>
#include "visualize.h"
#include "ws_ctube.h"
#include "mesh.h"

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
static int32_t lin_scale(number x, number lo, number hi, int32_t Lo, int32_t Hi)
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
	number aa = *(number *)a;
	number bb = *(number *)b;
	return (aa > bb) - (aa < bb);
}

void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh)
{
	number min_xy = NUMBER_MAX;
	number max_xy = -NUMBER_MAX;

	number *stresses = malloc(mesh->nelements * sizeof(*stresses));

	/* determine min/max */
	for (int i = 0; i < mesh->nelements; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = mesh->elements[i]->vertices[j];
			for (int k = 0; k < DIM; k++) {
				number coord = v->pos.x[k];
				if (coord < min_xy) {
					min_xy = coord;
				}
				if (coord > max_xy) {
					max_xy = coord;
				}
			}
		}
		stresses[i] = mesh->elements[i]->scalar_stress;
	}

	/* percentile of stresses */
	qsort(stresses, mesh->nelements, sizeof(*stresses), number_cmp);
	number min_stress = stresses[(int)(0.05 * (mesh->nelements - 1))];
	number max_stress = stresses[(int)(0.95 * (mesh->nelements - 1))];

	for (int i = 0; i < mesh->nelements; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = mesh->elements[i]->vertices[j];
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

				number coord = v->pos.x[k];
				vis->data[i*(3*DIM + 1) + j*DIM + k] = lin_scale(coord, min_xy - 0.1*fabsf(min_xy), max_xy + 0.1*fabsf(max_xy), Lo, Hi);
			}
		}

		number stress = mesh->elements[i]->scalar_stress;
		vis->data[(i+1)*(3*DIM + 1) - 1] = lin_scale(stress, min_stress, max_stress, 0, 255);
	}
}

void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(vis->ctube, vis->data, vis->data_bytes);
}
