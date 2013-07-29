/*
 * output.c
 *
 *  Created on: 26.07.2011
 *      Author: mark
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "fasthash.h"
#include "output.h"


void output_DESEQ(uint8_t filectr, uint8_t argc, char *argv[]){

	nodeptr p;
	uint16_t k = 0;
	uint16_t j = 0;

	FILE *output = fopen("matrix.txt","w");
	fprintf(output,"name\t");
	for(k=0; k < argc-2; k++)
		fprintf(output,"%s\t",argv[k+1]);
	fprintf(output,"%s\n",argv[k+1]);

	for(k=0; k < NHASH; k++){
		for(p = bin[k]; p != NULL; p = p->next){
			if(strcmp(p->word,"") != 0){											/* sanity check: no emtpy ids */
				fprintf(output,"%s\t",p->word);
				for(j = 0; j < filectr-1; j++)
					fprintf(output,"%d\t",p->count[j] > 0 ? p->count[j] : 0);		/* sanity check: no negative values */
				fprintf(output,"%d\n",p->count[j] > 0 ? p->count[j] : 0);			/* sanity check: no negative values */
			}
		}
	}
	fclose(output);
}
