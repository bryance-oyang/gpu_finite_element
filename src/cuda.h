#ifndef CUDA_H
#define CUDA_H

int cuda_init(struct finite_element_problem *restrict p);
void cuda_destroy(struct finite_element_problem *restrict p);

#endif /* CUDA_H */
