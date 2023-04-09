/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @brief main
 */


#ifdef DEBUG
#define _GNU_SOURCE 1
#include <fenv.h>
#endif /* DEBUG */

#include <unistd.h>
#include <time.h>

#include "load_obj.h"
#include "visualize.h"
#include "finite_element.h"

static struct timespec start, end;

static void benchmark_start()
{
#ifdef CLOCK_MONOTONIC
		if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
			clock_gettime(CLOCK_REALTIME, &start);
		}
#else
		clock_gettime(CLOCK_REALTIME, &start);
#endif /* CLOCK_MONOTONIC */
}

static void benchmark_end(const char *msg)
{
#ifdef CLOCK_MONOTONIC
		if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
			clock_gettime(CLOCK_REALTIME, &end);
		}
#else
		clock_gettime(CLOCK_REALTIME, &end);
#endif /* CLOCK_MONOTONIC */

	double msec = (end.tv_sec - start.tv_sec)*1e3 + (end.tv_nsec - start.tv_nsec)*1e-6;
	printf(msg, msec);
}


int main()
{
#ifdef DEBUG
	feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
#endif /* DEBUG */

#ifndef GPU_COMPUTE
	printf("=== Running cpu only version (for gpu acceleration, compile with `make gpu`) ===\n");
#endif /* GPU_COMPUTE */

	struct mesh mesh;
	mesh_init(&mesh);

	printf("loading...\n");
	fflush(stdout);
	load_obj("../obj/beam.obj", &mesh, 16743);

	struct vis vis;
	vis_init(&vis, &mesh);

	struct finite_element_problem problem;
	printf("initing...\n");
	fflush(stdout);
	benchmark_start();
	if (fep_init(&problem, &mesh) != 0) {
		printf("error in fep_init!!!!!!!!!!!!!!\n");
		return -1;
	}
	benchmark_end("benchmarked build time: %g ms\n");

	printf("solving...\n");
	fflush(stdout);
	benchmark_start();
	if (fep_solve(&problem, 0.001, &vis) != 0) {
		printf("error in fep_solve!!!!!!!!!!!!!!\n");
		return -1;
	}
	benchmark_end("benchmarked solve time: %g ms\n");

#ifndef GPU_COMPUTE
	printf("=== Running cpu only version (for gpu acceleration, compile with `make gpu`) ===\n");
#endif /* GPU_COMPUTE */
	fflush(stdout);

	benchmark_start();
	vis_fill(&vis, &mesh, &problem.c);
	benchmark_end("benchmarked vis time: %g ms\n");
	for (;;) {
		vis_send(&vis);
		sleep(1);
	}
	vis_destroy(&vis);

	fep_destroy(&problem);
	mesh_destroy(&mesh);
	return 0;
}
