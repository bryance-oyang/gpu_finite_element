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
#include "load_obj.h"
#include "finite_element.h"
#include "visualize.h"

int main()
{
#ifdef DEBUG
	feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
#endif /* DEBUG */

	struct mesh mesh;
	mesh_init(&mesh);

	printf("loading...\n");
	fflush(stdout);
	load_obj("../obj/beam.obj", &mesh);
	for (int i = 0; i < mesh.nvertices; i++) {
		float x = mesh.vertices[i].pos.x[0];
		float y = mesh.vertices[i].pos.x[1];
		if (x < 0) {
			mesh.vertices[i].enabled = false;
		}
	}

	struct vis vis;
	vis_init(&vis, &mesh);

	struct finite_element_problem problem;
	printf("initing...\n");
	fflush(stdout);
	fep_init(&problem, &mesh);
	printf("solving...\n");
	fflush(stdout);
	if (fep_solve(&problem, 0.001, &vis) != 0) {
		printf("error in solve!!!!!!!!!!!!!!\n");
		return -1;
	}
	printf("calculating stresses...\n");
	fflush(stdout);
	fep_scalar_stress(&problem);
	printf("done\n");
	fflush(stdout);

/*
	vis_fill(&vis, &mesh);
	for (;;) {
		vis_send(&vis);
		sleep(1);
	}
	*/
	vis_destroy(&vis);

	fep_destroy(&problem);
	mesh_destroy(&mesh);
	return 0;
}
