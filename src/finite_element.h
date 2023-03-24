/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include "mesh.h"

/* Ac = b */
struct finite_element_problem {
	struct mesh *mesh;
	struct sparse A;
	struct vec c;
	struct vec b;

#ifdef GPU_COMPUTE

#endif /* GPU_COMPUTE */
};

int fep_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh);
void fep_destroy(struct finite_element_problem *restrict p);
void fep_solve(struct finite_element_problem *restrict p, number tolerance);

void fep_scalar_stress(struct finite_element_problem *restrict p);

#endif /* FINITE_ELEMENT_H */
