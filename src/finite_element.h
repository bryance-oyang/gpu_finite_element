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
