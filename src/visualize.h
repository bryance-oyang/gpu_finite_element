/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 */

#ifndef VISUALIZE_H
#define VISUALIZE_H

#include <stdlib.h>

struct mesh;
struct vec;
struct element;

struct vis {
	struct ws_ctube *ctube;
	size_t data_bytes;
	int32_t *data;

	float *stresses;

	int nsorted_stresses;
	float *sorted_stresses;
	struct element *cached_element;

	float min_x;
	float max_x;
	float mid_x;

	float min_y;
	float max_y;
	float mid_y;

	float slope;
};

int vis_init(struct vis *restrict vis, struct mesh *restrict mesh);
void vis_destroy(struct vis *restrict vis);
void vis_fill(struct vis *restrict vis, struct mesh *restrict mesh, struct vec *restrict c);
void vis_send(struct vis *restrict vis);

#endif /* VISUALIZE_H */
