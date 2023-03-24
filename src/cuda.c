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
