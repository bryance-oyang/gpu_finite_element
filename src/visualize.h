#ifndef VISUALIZE_H
#define VISUALIZE_H

#include <stdint.h>
#include <math.h>
#include "ws_ctube.h"
#include "mesh.h"

struct vis {
	struct ws_ctube *ctube;
	size_t data_bytes;
	int32_t *data;
};

static int vis_init(struct vis *restrict vis, struct mesh *restrict mesh)
{
	vis->data_bytes = (3*DIM + 1) * mesh->ntriangles * sizeof(*vis->data);
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

static void vis_destroy(struct vis *restrict vis)
{
	free(vis->data);
	ws_ctube_close(vis->ctube);
}

/* linearly map x in lo, hi to Lo, Hi */
static int32_t lin_scale(number x, number lo, number hi, int32_t Lo, int32_t Hi)
{
	if (hi - lo < EPSILON) {
		return 0;
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

static void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh)
{
	number min_xy = NUMBER_MAX;
	number max_xy = -NUMBER_MAX;

	number *stresses = malloc(mesh->ntriangles * sizeof(*stresses));

	/* determine min/max */
	for (int i = 0; i < mesh->ntriangles; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = &mesh->vertices[mesh->triangles[i].vertices[j]];
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
		stresses[i] = mesh->triangles[i].scalar_stress;
	}

	/* percentile of stresses */
	qsort(stresses, mesh->ntriangles, sizeof(*stresses), number_cmp);
	number min_stress = stresses[(int)(0.05 * mesh->ntriangles)];
	number max_stress = stresses[(int)(0.95 * mesh->ntriangles)];

	for (int i = 0; i < mesh->ntriangles; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = &mesh->vertices[mesh->triangles[i].vertices[j]];
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
				vis->data[i*(3*DIM + 1) + j*DIM + k] = lin_scale(coord, min_xy - 0.1*fabs(min_xy), max_xy + 0.1*fabs(max_xy), Lo, Hi);
			}
		}

		number stress = mesh->triangles[i].scalar_stress;
		vis->data[(i+1)*(3*DIM + 1) - 1] = lin_scale(stress, min_stress, max_stress, 0, 255);
	}
}

static void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(vis->ctube, vis->data, vis->data_bytes);
}

#endif /* VISUALIZE_H */
