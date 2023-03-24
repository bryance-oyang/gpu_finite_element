/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifdef DEBUG
#define _GNU_SOURCE 1
#include <fenv.h>
#endif /* DEBUG */

#include <stdio.h>
#include <unistd.h>
#include "finite_element.h"
#include "visualize.h"

int main()
{
#ifdef DEBUG
	feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
#endif /* DEBUG */

	struct mesh mesh;
	mesh_init(&mesh);

	int v0, v1, v2, v3, v4, v5, v6, v7;
	v0 = mesh_add_vertex(&mesh, 0, 1, false);
	v1 = mesh_add_vertex(&mesh, 0, 0, false);
	v2 = mesh_add_vertex(&mesh, 1, 1, true);
	v3 = mesh_add_vertex(&mesh, 1, 0, true);
	v4 = mesh_add_vertex(&mesh, 2, 1, true);
	v5 = mesh_add_vertex(&mesh, 2, 0, true);
	v6 = mesh_add_vertex(&mesh, 3, 1, true);
	v7 = mesh_add_vertex(&mesh, 3, 0, true);

	mesh_add_triangle(&mesh, v0, v1, v2, 1, 100);
	mesh_add_triangle(&mesh, v1, v2, v3, 1, 100);
	mesh_add_triangle(&mesh, v2, v3, v5, 1, 100);
	mesh_add_triangle(&mesh, v2, v4, v5, 1, 100);
	mesh_add_triangle(&mesh, v4, v5, v6, 1, 100);
	mesh_add_triangle(&mesh, v5, v6, v7, 1, 100);

	struct finite_element_problem problem;
	fep_init(&problem, &mesh);
	fep_solve(&problem, 0.00);
	fep_scalar_stress(&problem);

	for (int i = 0; i < DIM*mesh.nenabled; i++) {
		printf("%f\n", problem.c.x[i]);
	}

	printf("scalar_stresses:\n");
	for (int i = 0; i < mesh.ntriangles; i++) {
		printf("%f\n", mesh.triangles[i].scalar_stress);
	}

	struct vis vis;
	vis_init(&vis, &mesh);
	vis_fill(&vis, &mesh);
	for (;;) {
		vis_send(&vis);
		usleep(10000);
	}

	vis_destroy(&vis);
	fep_destroy(&problem);
	mesh_destroy(&mesh);
	return 0;
}
