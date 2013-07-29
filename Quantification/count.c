/*
 * count.c
 *
 *  Created on: 22.11.2011
 *      Author: markonyango
 */

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fasthash.h"
#include "output.h"


clock_t BeginTimer()
{
    /* timer declaration */
    clock_t Begin; /* initialize Begin */

    Begin = clock() ; /* start the timer */

    return Begin;
}

clock_t EndTimer(clock_t begin)
{
    clock_t End;
    End = clock() ;   /* stop the timer */
    return End;
}

int main(int argc, char *argv[]) {

	double begin = 0;

	uint8_t tokenctr = 0;
	uint8_t filectr = 0;
	uint8_t mapped = 0;		/* Boolean: is read mapped or not? */
	uint32_t readcount_all = 0;
	uint32_t readcount_unmapped = 0;
	uint32_t readcount_mapped = 0;
	unsigned total = 0;
	long k=0;
	char *ptr;
	char line[1024];
	char *temp_string;
	char delimiter[] ="\t";

	float elapTicks;
	float elapMilli, elapSeconds;
	nodeptr p;

	FILE* file;

	if(argc == 1){
		printf("Usage:\n");
		exit(1);
	}


	for(filectr = 0; filectr < (argc-1); filectr++){

		begin  = BeginTimer();

		printf("%d. Extracting IDs from %s (this may take some time)\n",filectr+1, argv[filectr+1]);

		file = fopen(argv[filectr+1], "r");
		if (file==NULL) {fputs ("File error!\n",stderr); exit (1);}
		while(fgets(line,sizeof(line),file) != NULL){
			/* We found a record in the SAM file
			 * lets check that line token by token
			 */

			if(line[1] != '@' && line[0] != '@'){

				readcount_all++;

				ptr = strtok(line,delimiter);
				tokenctr++;
				while(ptr != NULL && tokenctr < 4){
					if(strcmp(ptr,"*") != 0 && (tokenctr == 3)){
						readcount_mapped++;
						mapped = 1;
						incword(ptr, filectr);
					}
					else if(strcmp(ptr,"*") == 0 && (tokenctr ==3)){
						readcount_unmapped++;
					}
					ptr = strtok(NULL,delimiter);
					tokenctr++;
				}
				tokenctr = 0;
				mapped = 0;
			}
			/* We found a SAM header line
			 * lets index the id
			 */
			else{
				if(filectr == 0){
					ptr = strtok(line,delimiter);
					if(ptr[1] == 'S'){
						ptr = strtok(NULL,delimiter);

						temp_string = strndup(ptr+3,strlen(ptr));
						incword(temp_string, filectr);
						free(temp_string); /* fixes memory leak */
					}
				}
			}
		}
		fclose(file);


		for (k = 0; k < NHASH; k++){
				for (p = bin[k]; p != NULL; p = p->next){
					total += p->count[filectr] ? 1 : 0;
				}
		}

		printf("Distinct IDs: %d\n", total);
		printf("Total Reads found: %d\n", readcount_all);
		printf("Unmapped Reads found: %d (%.2f%)\n", readcount_unmapped, ((float)readcount_unmapped/readcount_all)*100);
		printf("Mapped Reads found: %d (%.2f%)\n", readcount_mapped, ((float)readcount_mapped/readcount_all)*100);

		readcount_all = 0;
		readcount_unmapped = 0;
		readcount_mapped = 0;
		total = 0;

		printf("==================================================\n");
	}

	/* variable definitions to calculate time taken */
	elapTicks = EndTimer(begin);    /* stop the timer, and calculate the time taken */
	elapMilli = elapTicks/1000;     /* milliseconds from Begin to End */
	elapSeconds = elapMilli/1000;   /* seconds from Begin to End */
	/* elapMinutes = elapSeconds/60;    minutes from Begin to End */

	printf ("Seconds passed: %.2f\n\n", elapSeconds);

	printf("Calculating gene expression levels\n");
	output_DESEQ(filectr,argc,argv);

	return 0;
}
