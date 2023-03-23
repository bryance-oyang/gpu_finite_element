#ifndef VISUALIZE_H
#define VISUALIZE_H

#include "ws_ctube.h"
#include "mesh.h"

struct vis {
	struct ws_ctube ctube;
	size_t data_bytes;
	float *data;
};

static int vis_init(struct vis *restrict vis, struct mesh *restrict mesh)
{
	vis->data_bytes = (3*DIM + 1) * mesh->ntriangles * sizeof(*vis->data);
	vis->data = malloc(vis->data_bytes);
	if (vis->data == NULL) {
		goto err_nodata;
	}
	if (ws_ctube_init(&vis->ctube, 9743, 2, 0, 24) != 0) {
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
	ws_ctube_destroy(&vis->ctube);
}

static void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh)
{
	for (int i = 0; i < mesh->ntriangles; i++) {
		for (int j = 0; j < 3; j++) {
			struct vertex *v = &mesh->vertices[mesh->triangles[i].vertices[j]];
			for (int k = 0; k < DIM; k++) {
				vis->data[i*(3*DIM + 1) + j*DIM + k] = v->pos.x[k];
			}
		}
		vis->data[(i+1)*(3*DIM + 1) - 1] = mesh->triangles[i].scalar_stress;
	}
}

static void vis_send(struct vis *restrict vis)
{
	ws_ctube_broadcast(&vis->ctube, vis->data, vis->data_bytes);
}

#endif /* VISUALIZE_H */
