/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "finite_element.h"

#ifdef GPU_COMPUTE
#include <cuda_runtime.h>
#include <cusparse.h>
#endif /* GPU_COMPUTE */

int cuda_init(struct finite_element_problem *restrict p)
{
#ifdef GPU_COMPUTE
#else /* GPU_COMPUTE */
	(void)p;
#endif /* GPU_COMPUTE */

	return 0;
}

void cuda_destroy(struct finite_element_problem *restrict p)
{
#ifdef GPU_COMPUTE
#else /* GPU_COMPUTE */
	(void)p;
#endif /* GPU_COMPUTE */
}
