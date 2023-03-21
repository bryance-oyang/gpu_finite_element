/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <stdio.h>
#include "element.h"

//     0     1     2     3
// --- * --- * --- * --- *

int main()
{
	struct sparse A;
	struct vec gravity;
	struct vec answer;

	sparse_init(&A);
	vec_init(&gravity, 4);
	vec_init(&answer, 4);

	gravity.x[0] = -1;
	gravity.x[1] = -1;
	gravity.x[2] = -1;
	gravity.x[3] = -0.5;

	insert_boundary(&A, 0, 1, 1);
	insert_line(&A, 0, 1, 1, 1);
	insert_line(&A, 1, 2, 1, 1);
	insert_line(&A, 2, 3, 1, 1);

	sparse_conj_grad(0.08, &A, &gravity, &answer);

	for (int i = 0; i < 4; i++) {
		printf("%f\n", answer.x[i]);
	}

	vec_destroy(&answer);
	vec_destroy(&gravity);
	sparse_destroy(&A);
	return 0;
}
