/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 */

#ifndef CUDA_H
#define CUDA_H

#ifdef GPU_COMPUTE

int cuda_init(struct solver *restrict solver);
void cuda_destroy(struct solver *restrict solver);

int gpu_conj_gradient(struct solver *restrict solver, double tolerance);

#endif /* GPU_COMPUTE */

#endif /* CUDA_H */
