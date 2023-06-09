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

#ifndef LOAD_OBJ_H
#define LOAD_OBJ_H

#include <stdlib.h>
#include <stdio.h>

#include "elements.h"

int load_obj(char *filename, struct mesh *mesh, int elasticity)
{
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		return -1;
	}

	char *line = NULL;
	size_t line_size;

	while (getline(&line, &line_size, file) > 0) {
		double n0, n1, n2;
		int v0, v1, v2, _v0, _v1, _v2, v0_, v1_, v2_;

		if (sscanf(line, "v %lf %lf %lf", &n0, &n1, &n2) == 3) {
			bool enabled = true;
			double x = n0;
			double y = -n2;

			if (x < 0) {
				enabled = false;
			}
			mesh_add_vertex(mesh, x, y, enabled);
		} else if (sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d", &v0, &_v0, &v0_, &v1, &_v1, &v1_, &v2, &_v2, &v2_) == 9) {
			mesh_add_triangle3(mesh, v0-1, v1-1, v2-1, 1, elasticity);
		}
	}

	if (line != NULL) {
		free(line);
	}
	fclose(file);
	return 0;
}

#endif /* LOAD_OBJ_H */
