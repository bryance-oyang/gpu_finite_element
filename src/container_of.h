/**
 * @file
 */

#ifndef CONTAINER_OF_H
#define CONTAINER_OF_H

#include <stddef.h>

/**
 * get the container of ptr which is a member of a struct
 *
 * Example:
 *
 * struct a {
 * 	struct b bmember;
 * 	struct c cmember;
 * };
 *
 * struct a *a_ptr;
 * struct c *c_ptr = &a_ptr->cmember;
 *
 * assert(a_ptr == container_of(c_ptr, struct a, cmember))
 */
#define container_of(ptr, type, member) ({ \
	typeof(((type *)0)->member) *_container_of_ptr = (ptr); \
	((type *)((char *)(_container_of_ptr) - offsetof(type, member)));})

#endif /* CONTAINER_OF_H */
