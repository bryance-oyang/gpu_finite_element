/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef CUDA_H
#define CUDA_H

int cuda_init(struct finite_element_problem *restrict p);
void cuda_destroy(struct finite_element_problem *restrict p);

#endif /* CUDA_H */
