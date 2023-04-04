#include "hash_table.h"

int ht_init(struct ht *ht)
{
	(void)ht;
	return 0;
}

/*
void ht_destroy(struct ht *ht, bool free_entries)
{
	if (!free_entries) {
		return;
	}

	struct ht_node *prev, *cur;
	for (int i = 0; i < HT_SIZE; i++) {
		if (ht->table[i].next == NULL) {
			continue;
		}
		cur = ht->table[i].next;
		while(cur != NULL) {
			prev = cur;
			cur = cur->next;
			free(container_of(...));
		}
	}
}
*/

unsigned int knuth_hash(const char *str)
{
	unsigned int len = strlen(str);
	unsigned int hash = len;
	for (unsigned int i = 0; i < len; str++, i++) {
		hash = ((hash << 5) ^ (hash >> 27)) ^ (*str);
	}
	return hash;
}
