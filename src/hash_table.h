#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <stdlib.h>
#include <string.h>
#include "container_of.h"

#define HT_SIZE (1 << 16)

struct ht_node {
	struct ht_node *next;
};

struct ht {
	struct ht_node table[HT_SIZE];
};

int ht_init(struct ht *ht)
{
	return 0;
}

void ht_destroy(struct ht *ht)
{
	struct ht_node *prev, *cur;
	for (int i = 0; i < HT_SIZE; i++) {
		if (ht->table[i].next == NULL) {
			continue;
		}
		cur = ht->table[i].next;
		while(cur != NULL) {
			prev = cur;
			cur = cur->next;
			free(prev);
		}
	}
}

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

unsigned int knuth_hash(const char *str)
{
	unsigned int len = strlen(str);
	unsigned int hash = len;
	for (unsigned int i = 0; i < len; str++, i++) {
		hash = ((hash << 5) ^ (hash >> 27)) ^ (*str);
	}
	return hash;
}

#endif /* HASH_TABLE_H */
