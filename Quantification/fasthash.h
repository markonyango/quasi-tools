/*
 * fasthash.h
 *
 *  Created on: 25.07.2011
 *      Author: mark
 */

#ifndef FASTHASH_H_
#define FASTHASH_H_

#define NHASH 29989
#define MULT 31
#define NODEGROUP 1000

typedef struct nodes *nodeptr;

typedef struct nodes {
	char *word;
	uint16_t count[256];
	nodeptr next;
} nodes;

unsigned int fhash(char *p);
nodeptr nmalloc();
char *smalloc(int n);
void incword(char *s, uint8_t files);

nodeptr bin[NHASH];
nodeptr freenode;


#endif /* FASTHASH_H_ */
