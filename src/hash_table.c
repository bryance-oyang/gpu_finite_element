/**
 * @file
 * @brief hash table
 */

#include <stdlib.h>
#include <string.h>

#include "hash_table.h"

int ht_init(struct ht *ht)
{
	for (int i = 0; i < HT_SIZE; i++) {
		ht->table[i].next = NULL;
	}
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

void ht_node_init(struct ht_node *node)
{
	node->next = NULL;
}

uint32_t knuth_hash(const char *str)
{
	uint32_t mask = 0;
	mask = ~mask;
	mask <<= HT_BITS;
	mask = ~mask;

	uint32_t len = strlen(str);
	uint32_t hash = len;
	for (uint32_t i = 0; i < len; str++, i++) {
		hash = ((hash << 5) ^ (hash >> 27)) ^ (*str);
	}
	return hash & mask;
}
