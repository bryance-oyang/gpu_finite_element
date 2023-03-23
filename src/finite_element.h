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
	struct sparse A;
	struct vec c;
	struct vec b;
};

int finite_element_problem_init(struct finite_element_problem *restrict p, struct mesh *restrict mesh);
void finite_element_problem_destroy(struct finite_element_problem *restrict p);
void finite_element_problem_solve(struct finite_element_problem *restrict p, number tolerance);

#endif /* FINITE_ELEMENT_H */
