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

#include "solver.h"
#include "cuda.h"

int solver_init(struct solver *restrict solver, struct mesh *restrict mesh)
{
	if (sparse_init(&solver->A) != 0) {
		goto err_noA;
	}
	if (vec_init(&solver->b, DIM * mesh->nenabled) != 0) {
		goto err_nob;
	}
	if (vec_init(&solver->c, DIM * mesh->nenabled) != 0) {
		goto err_noc;
	}
	solver->mesh = mesh;


	mesh_assign_vertex_ids(mesh);
	mesh_construct_problem(mesh, &solver->A, &solver->b);

	sparse_sort(&solver->A);
	sparse_consolidate(&solver->A);

#ifdef GPU_COMPUTE
	if (cuda_init(solver) != 0) {
		goto err_nocuda;
	}
#endif /* GPU_COMPUTE */

	return 0;

#ifdef GPU_COMPUTE
err_nocuda:
	vec_destroy(&p->c);
#endif /* GPU_COMPUTE */

err_noc:
	vec_destroy(&solver->b);
err_nob:
	sparse_destroy(&solver->A);
err_noA:
	return -1;
}

void solver_destroy(struct solver *restrict solver)
{
#ifdef GPU_COMPUTE
	cuda_destroy(solver);
#endif
	vec_destroy(&solver->c);
	vec_destroy(&solver->b);
	sparse_destroy(&solver->A);
}

int solver_solve(struct solver *restrict solver, double tolerance, struct vis *restrict vis)
{
	int retval;

#ifdef GPU_COMPUTE
	retval = gpu_conj_gradient(solver, tolerance);
#else /* GPU_COMPUTE */
	retval = sparse_conj_grad(&solver->A, &solver->b, &solver->c, tolerance, vis, solver->mesh);
#endif /* GPU_COMPUTE */

	return retval;
}
