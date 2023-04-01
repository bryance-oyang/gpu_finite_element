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

#include "mesh.h"

int load_obj(char *filename, struct mesh *mesh)
{
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		return -1;
	}

	char *line = NULL;
	size_t line_size;

	while (getline(&line, &line_size, file) > 0) {
		float n0, n1, n2;
		int v0, v1, v2, _v0, _v1, _v2, v0_, v1_, v2_;

		if (sscanf(line, "v %f %f %f", &n0, &n1, &n2) == 3) {
			mesh_add_vertex(mesh, n0, -n2, true);
		} else if (sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d", &v0, &_v0, &v0_, &v1, &_v1, &v1_, &v2, &_v2, &v2_) == 9) {
			mesh_add_triangle(mesh, &mesh->vertices[v0-1],
				&mesh->vertices[v1-1], &mesh->vertices[v2-1],
				1, 1);
		}
	}

	if (line != NULL) {
		free(line);
	}
	fclose(file);
	return 0;
}

#endif /* LOAD_OBJ_H */
