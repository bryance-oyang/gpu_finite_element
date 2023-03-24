#ifndef VISUALIZE_H
#define VISUALIZE_H

#include <stdlib.h>
#include <stdint.h>
#include "mesh.h"

struct vis {
	struct ws_ctube *ctube;
	size_t data_bytes;
	int32_t *data;
};

int vis_init(struct vis *restrict vis, struct mesh *restrict mesh);
void vis_destroy(struct vis *restrict vis);
void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh);
void vis_send(struct vis *restrict vis);

#endif /* VISUALIZE_H */
