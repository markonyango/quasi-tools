/* Copyright (C) 1999 Lucent Technologies */
/* From 'Programming Pearls' by Jon Bentley */

/* wordfreq.c -- list of words in file, with counts */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "fasthash.h"

/* A few global variables to be seen from everywhere */
int nodesleft = 0;
nodeptr freenode = NULL;
nodeptr bin[NHASH];

unsigned int fhash(char *p)
{	unsigned int h = 0;
	for ( ; *p; p++)
		h = MULT * h + *p;
	return h % NHASH;
}

nodeptr nmalloc()
{	if (nodesleft == 0) {
		freenode = malloc(NODEGROUP*sizeof(nodes));
		nodesleft = NODEGROUP;
	}
	nodesleft--;
	return freenode++;
}

#define CHARGROUP 10000
int charsleft = 0;
char *freechar;

char *smalloc(int n)
{	if (charsleft < n) {
		freechar = malloc(n+CHARGROUP);
		charsleft = n+CHARGROUP;
	}
	charsleft -= n;
	freechar += n;
	return freechar - n;
}

void incword(char *s, uint8_t files)
{	nodeptr p;
	int h = fhash(s);
	for (p = bin[h]; p != NULL; p = p->next)
		if (strcmp(s, p->word) == 0) {
			(p->count[files])++;
			return;
		}

	p = nmalloc();
	p->count[files] = 0;
	p->word = smalloc(strlen(s)+1);
	strcpy(p->word, s);
	p->next = bin[h];
	bin[h] = p;
}
