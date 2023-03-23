/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include "linear_algebra.h"

static inline int insert_line(struct sparse *restrict S, int i, int j, number len, number E)
{
	int retval = 0;
	number k = E / len;
	retval += sparse_add(S, i, i, k);
	retval += sparse_add(S, i, j, -k);
	retval += sparse_add(S, j, i, -k);
	retval += sparse_add(S, j, j, k);
	return retval;
}

static inline int insert_boundary(struct sparse *restrict S, int i, number len, number E)
{
	number k = E / len;
	return sparse_add(S, i, i, k);
}

#endif /* ELEMENT_H */
