#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "container_of.h"

#define HT_BITS 16
#define HT_SIZE (1 << (HT_BITS))

struct ht_node {
	struct ht_node *next;
};

struct ht {
	struct ht_node table[HT_SIZE];
};

#define ht_for_each(ht, cur, prev, key) \
	for (prev = &(ht)->table[(key)], cur = (ht)->table[(key)].next; cur != NULL; prev = cur, cur = cur->next)

/**
 * try to find ptr's parent in ht. if found, cur will point to it. if not found,
 * ptr can be insert as prev->next
 */
#define ht_find(ht, cur, prev, key, ptr, type, member, eq_fn, found) \
	do { \
		found = 0; \
		ht_for_each(ht, cur, prev, key) {\
			if (eq_fn(container_of(ptr, type, member), container_of(cur, type, member))) { \
				found = 1; \
				break; \
			} \
		} \
	} while (0);

int ht_init(struct ht *ht);
void ht_destroy(struct ht *ht);
void ht_node_init(struct ht_node *node);
uint32_t knuth_hash(const char *str);

#endif /* HASH_TABLE_H */
