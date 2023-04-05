/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @brief setup finite element problem
 */

#include <stdio.h>
#include <time.h>
#include <math.h>

#include "finite_element.h"
#include "cuda.h"

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

	float msec = (end.tv_sec - start.tv_sec)*1e3 + (end.tv_nsec - start.tv_nsec)*1e-6;
	printf(msg, msec);
}

int fep_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh)
{
	if (sparse_init(&p->A) != 0) {
		goto err_noA;
	}
	if (vec_init(&p->b, DIM * mesh->nenabled) != 0) {
		goto err_nob;
	}
	if (vec_init(&p->c, DIM * mesh->nenabled) != 0) {
		goto err_noc;
	}
	p->mesh = mesh;


	benchmark_start();
	mesh_assign_vertex_ids(mesh);
	mesh_construct_problem(mesh, &p->A, &p->b);
	benchmark_end("benchmarked mesh build time: %g ms\n");

	benchmark_start();
	sparse_sort(&p->A);
	sparse_consolidate(&p->A);
	benchmark_end("benchmarked sparse sort time: %g ms\n");

#ifdef GPU_COMPUTE
	if (cuda_init(p) != 0) {
		goto err_nocuda;
	}
#endif /* GPU_COMPUTE */

	return 0;

#ifdef GPU_COMPUTE
err_nocuda:
	vec_destroy(&p->c);
#endif /* GPU_COMPUTE */

err_noc:
	vec_destroy(&p->b);
err_nob:
	sparse_destroy(&p->A);
err_noA:
	return -1;
}

void fep_destroy(struct finite_element_problem *restrict p)
{
#ifdef GPU_COMPUTE
	cuda_destroy(p);
#endif
	vec_destroy(&p->c);
	vec_destroy(&p->b);
	sparse_destroy(&p->A);
}

int fep_solve(struct finite_element_problem *restrict p, float tolerance, struct vis *restrict vis)
{
	int retval;

	benchmark_start();

#ifdef GPU_COMPUTE
	retval = gpu_conj_gradient(p, tolerance);
#else /* GPU_COMPUTE */
	retval = sparse_conj_grad(&p->A, &p->b, &p->c, tolerance, vis, p->mesh);
#endif /* GPU_COMPUTE */

	benchmark_end("benchmarked solve time: %g ms\n");
	return retval;
}
