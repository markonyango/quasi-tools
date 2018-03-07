/*
 ============================================================================
 Name        : qa.c
 Author      : Mark Onyango
 Version     :
 Copyright   : 
 Description : 
 ============================================================================
 */

#include "qa.h"

int main(int argc, char *argv[])
{

	/* Each output file is give it's unique filepointer
	 * to avoid confusion during output.
	 */
	FILE *fileptr;
	FILE *boxplot_fileptr;
	FILE *phred_dist_fileptr;
	FILE *length_dist_fileptr;
	FILE *bases_fileptr;

	/* We need to keep track of the current working directory
	 * so that output gets written to the directory this tool
	 * was called from.
	 */
	char cwd[1024];

	/* The size of the line array should account for very long
	 * reads in addition to the readname etc. 1024 characters
	 * should suffice for most FASTQ files.
	 */
	char line[1024];

	/* For the assessment of tile quality it is important to know
	 * the Illumina FASTQ version since they show differences in
	 * the Name section (included whitespaces, FC identifier, etc)
	 */
	/* uint8_t illumina_fastq_version; // never used */

	uint32_t linecount;
	uint16_t i, j;
	uint32_t readcount;
	uint64_t *dist_phred;
	uint64_t *dist_phred_ATGC;
	uint32_t *dist_seqlength;
	uint32_t *boxplot_bins;

	float elapTicks;
	float elapMilli, elapSeconds;
	double begin = 0;

	char tempFileName[2048];
	char *phred_sequence;

	char buff[PATH_MAX];

	/* Count how many times each base was sequenced */
	uint64_t sum_A = 0;
	uint64_t sum_T = 0;
	uint64_t sum_G = 0;
	uint64_t sum_C = 0;
	uint64_t sum_N = 0;

	base *bases;
	bases = NULL;

	dist_phred = calloc(200, sizeof(uint64_t));
	dist_phred_ATGC = calloc(5000, sizeof(uint64_t));

	linecount = 0;
	i = 0;
	max_seqlength = 0;

	/* Let's get the current working directory first or exit with 1 */
	if (getcwd(cwd, sizeof(cwd)) == NULL)
	{
		exit(1);
	}
	else
	{
		printf("Current working directory: %s\n", cwd);
		printf("Basename of file: %s\n", basename(argv[1]));
		chdir("..");
		printf("Basename: %s\n", basename(getcwd(buff, sizeof(buff))));
	}

	/* Open the input file in read mode */
	fileptr = fopen(argv[1], "r");
	if (fileptr == NULL)
	{
		fputs("You neglected to specify a input file!\n", stderr);
		exit(1);
	} else {
		printf("Processing file %s...\n", argv[1]);
	}

	begin = BeginTimer();

	phred_offset = detect_phred(fileptr);

	while (fgets(line, sizeof(line), fileptr) != NULL)
	{
		linecount++;

		if ((linecount % 2 == 0) && (linecount % 4 != 0))
		{
			i = strlen(line) - 1;
			max_seqlength = max_seqlength < i ? i : max_seqlength;
		}
	}

	linecount = 0;
	rewind(fileptr);

	boxplot_bins = calloc(max_seqlength * 50, sizeof(uint32_t));
	dist_seqlength = calloc(max_seqlength, sizeof(uint32_t));
	bases = malloc(max_seqlength * sizeof(base));
	phred_sequence = malloc((max_seqlength * sizeof(char)) + 1);
	number_of_bases_per_cycle = calloc(max_seqlength, sizeof(uint64_t));

	while (fgets(line, sizeof(line), fileptr) != NULL)
	{

		linecount++;

		if ((linecount % 2 == 0) && (linecount % 4 != 0))
		{
			/* Since i has not been used before we can use it as
			 * a temporary variable as it will be re-initialized
			 * with 0 further down.
			 */
			i = strlen(line) - 1;
			dist_seqlength[i] = dist_seqlength[i] + 1;

			base_content(bases, line, i);
			/* Save the read for further use when calculating the phred
             * score distribution per base.
             */
			strcpy(phred_sequence, line);
		}
		else if (linecount % 4 == 0)
		{
			distribution_phred(dist_phred, dist_phred_ATGC, line, i, phred_sequence);
			boxplot(boxplot_bins, line, i);

			/* Put code for progress bar here as this line is
			 * always the last line of a read block.*/
		}
	}
	fclose(fileptr);

	/* Number of reads that were just read */
	readcount = linecount / 4;

	/* Write the distribution of bases per cycle to file.
	 * Also write the GC-content per cycle to file.
	 */
	i = sprintf(tempFileName, "%s/%s-base_dist.txt", cwd, basename(argv[1]));
	bases_fileptr = fopen(tempFileName, "w");
	if (bases_fileptr == NULL)
	{
		printf("Warning: Could not allocate space for _bases_dist.txt!\n");
	}
	else
	{
		for (i = 0; i < max_seqlength; i++)
		{
			fprintf(bases_fileptr, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
					(float)bases[i].A / (readcount - (readcount - number_of_bases_per_cycle[i])) * 100,
					(float)bases[i].T / (readcount - (readcount - number_of_bases_per_cycle[i])) * 100,
					(float)bases[i].G / (readcount - (readcount - number_of_bases_per_cycle[i])) * 100,
					(float)bases[i].C / (readcount - (readcount - number_of_bases_per_cycle[i])) * 100,
					(float)bases[i].N / (readcount - (readcount - number_of_bases_per_cycle[i])) * 100,
					((float)bases[i].G + bases[i].C) / (readcount - (readcount - number_of_bases_per_cycle[i])) * 100);
		}
	}
	fclose(bases_fileptr);

	/* Write the distribution of phred scores to file.
	 * Phred scores range from 0 to 60 so only those will be output 
     */
	i = sprintf(tempFileName, "%s/%s-phred_dist.txt", cwd, basename(argv[1]));
	phred_dist_fileptr = fopen(tempFileName, "w");
	if (phred_dist_fileptr == NULL)
	{
		printf("Warning: Could not allocate space for phred_dist.txt!\n");
	}
	else
	{
		for (i = 0; i < max_seqlength; i++)
		{
			sum_A += bases[i].A;
			sum_T += bases[i].T;
			sum_G += bases[i].G;
			sum_C += bases[i].C;
			sum_N += bases[i].N;
		}

		/* Write the distribution for all bases' phred scores to file
        for(i = 0; i < 61; i++){
			double temp = ((double)dist_phred[i]/(readcount*seqlength))*100;
			fprintf(phred_dist_fileptr,"%f ",temp);
		}*/

		/* Write the distribution for A phred scores to file */
		for (i = 0; i < 61; i++)
		{
			double temp = ((double)dist_phred_ATGC[(4 * i) + 0] / sum_A) * 100;
			fprintf(phred_dist_fileptr, "%f ", temp);
		}

		/* Write the distribution for T phred scores to file */
		for (i = 0; i < 61; i++)
		{
			double temp = ((double)dist_phred_ATGC[(4 * i) + 1] / sum_T) * 100;
			fprintf(phred_dist_fileptr, "%f ", temp);
		}

		/* Write the distribution for G phred scores to file */
		for (i = 0; i < 61; i++)
		{
			double temp = ((double)dist_phred_ATGC[(4 * i) + 2] / sum_G) * 100;
			fprintf(phred_dist_fileptr, "%f ", temp);
		}

		/* Write the distribution for C phred scores to file */
		for (i = 0; i < 61; i++)
		{
			double temp = ((double)dist_phred_ATGC[(4 * i) + 3] / sum_C) * 100;
			fprintf(phred_dist_fileptr, "%f ", temp);
		}
	}
	fclose(phred_dist_fileptr);

	/* Write the distribution of sequence lengths to file.
	 * The for loop ranges from 0 to max_seqlength.
	 */
	//strcpy(tempFileName, cwd);
	i = sprintf(tempFileName, "%s/%s-length_dist.txt", cwd, basename(argv[1]));
	length_dist_fileptr = fopen(tempFileName, "w");
	if (length_dist_fileptr == NULL)
	{
		printf("Warning: Could not allocate space for length_dist.txt!\n");
	}
	else
	{
		for (i = 1; i <= max_seqlength; i++)
		{
			/*fwrite((const void*)&dist_seqlength[i], sizeof(uint32_t), 1, length_dist_fileptr);*/
			fprintf(length_dist_fileptr, "%d ", dist_seqlength[i]);
		}
	}
	fclose(length_dist_fileptr);

	/* Write the bins for the boxplot to file.
     * Each array subset corresponds to a cycle and ranges from 1 - 50 containing the number of times
     * each quality score was counted at the corresponding cycle.
     */
	i = sprintf(tempFileName, "%s/%s-boxplotdata.txt", cwd, basename(argv[1]));
	boxplot_fileptr = fopen(tempFileName, "w");
	if (boxplot_fileptr == NULL)
	{
		printf("Boxplot File could not be created!\n");
	}
	for (i = 0; i < 50; i++)
	{
		for (j = 0; j < max_seqlength; j++)
		{
			fprintf(boxplot_fileptr, "%d ", boxplot_bins[i * max_seqlength + j]);
		}
		fprintf(boxplot_fileptr, "\n");
	}
	fclose(boxplot_fileptr);

	printf("Reads: %d\n", readcount);
	printf("Max. Length: %d\n", max_seqlength);

	/* variable definitions to calculate time taken */
	elapTicks = EndTimer(begin);	/* stop the timer, and calculate the time taken */
	elapMilli = elapTicks / 1000;   /* milliseconds from Begin to End */
	elapSeconds = elapMilli / 1000; /* seconds from Begin to End */
	/* elapMinutes = elapSeconds/60;    minutes from Begin to End */

	printf("Seconds passed: %.2f\n", elapSeconds);
	/* printf("Milliseconds passed: %.2f\n\n",elapMilli); */

	/* free(bases);	24.07.2012 - This caused invalid calls to free memory under Ubuntu */
	return 0;
}

uint32_t countReads(FILE *fileptr)
{
	char line[1024];
	uint32_t readCounts = 0;

	while ((fgets(line, sizeof(line), fileptr) != NULL))
	{
		readCounts++;
	}
	rewind(fileptr);
	return (readCounts / 4);
}

void boxplot(uint32_t *boxplot_bins, char *line, uint16_t length)
{

	uint16_t cycle = 0;
	/*uint16_t length = strlen(line)-1;*/
	uint8_t phred = 0;

	for (cycle = 0; cycle < length; cycle++)
	{
		phred = line[cycle] - phred_offset; /* Get the phred score from the array */
											/* Arrays are accessed in row-major order: array[x][y] = row * NUMCOLS + col */
		boxplot_bins[(phred * max_seqlength) + cycle]++;
	}
}

uint8_t detect_phred(FILE *fileptr)
{

	uint16_t i = 0;
	uint16_t length;
	uint8_t lowest_char = 255;
	char line[1024];
	uint16_t linecount = 0, reads = 0;

	while ((fgets(line, sizeof(line), fileptr) != NULL) && reads <= 10000)
	{
		linecount++;
		if (linecount % 4 == 0)
		{
			reads++;
			length = strlen(line) - 1;
			for (i = 0; i < length; i++)
			{
				lowest_char = line[i] < lowest_char ? line[i] : lowest_char;
			}
		}
	}

	rewind(fileptr);

	/* Detect phred-score Version */
	if (lowest_char < 59)
	{
		printf("Sanger / Illumina 1.9 encoding detected\n");
		return (33); /* Quality is Sanger encoded */
	}
	else if (lowest_char < 64)
	{
		printf("Solexa encoding detected\n");
		return (59); /* Quality is Solexa encoded */
	}
	else if (lowest_char < 66)
	{
		printf("Illumina 1.3+ encoding detected\n");
		return (64); /* Quality is Illumina 1.3+ encoded */
	}
	else if (lowest_char <= 126)
	{
		printf("Illumina 1.5+ encoding detected\n");
		return (64); /* Quality is Illumina 1.5+ encoded */
	}

	printf("Quality encoding unknown! Setting encoding to Sanger per default!\n");
	return (33); /* Quality encoding could not be determined. Set to Sanger per default */
}

void distribution_phred(uint64_t *dist_phred, uint64_t *dist_phred_ATGC, char *cigar, uint16_t length, char *phred_sequence)
{
	uint16_t i = 0;
	uint8_t score = 0;

	for (i = 0; i < length; i++)
	{
		score = (uint8_t)cigar[i] - phred_offset;
		dist_phred[score]++;
		switch (phred_sequence[i])
		{
		case 'A':
			dist_phred_ATGC[(4 * score) + 0]++;
			break;
		case 'T':
			dist_phred_ATGC[(4 * score) + 1]++;
			break;
		case 'G':
			dist_phred_ATGC[(4 * score) + 2]++;
			break;
		case 'C':
			dist_phred_ATGC[(4 * score) + 3]++;
			break;
		default:
			break;
		}
	}
}

void base_content(base *bases, char *read, uint16_t length)
{
	/*uint16_t length = strlen(read)-1;*/
	uint16_t i = 0;

	for (i = 0; i < length; i++)
	{
		switch (read[i])
		{
		case 'A':
			bases[i].A++;
			break;
		case 'T':
			bases[i].T++;
			break;
		case 'G':
			bases[i].G++;
			break;
		case 'C':
			bases[i].C++;
			break;
		case 'N':
			bases[i].N++;
			break;
		case '.': /* 24.11.2011: Recognizing . as N */
			bases[i].N++;
			break;
		default:
			printf("Base couldn't be identified! (%d)\n", i);
			break;
		}
		number_of_bases_per_cycle[i] = number_of_bases_per_cycle[i] + 1;
	}
}

uint8_t detect_fastq_version(char *line)
{
	uint16_t i = 0;
	size_t len = strlen(line);

	for (i = 0; i < len; i++)
		if (line[i] == ' ')
			return (1);

	return (0);
}

clock_t BeginTimer()
{
	/* timer declaration */
	clock_t Begin; /* initialize Begin */

	Begin = clock(); /* start the timer */

	return Begin;
}

clock_t EndTimer(clock_t begin)
{
	clock_t End;
	End = clock(); /* stop the timer */
	return End;
}

/*char * strndup(const char *s, size_t n){
	size_t len = strnlen(s,n);
	char *new = malloc(len+1);

	if(new==NULL)
		return NULL;

	new[len]='\0';
	return memcpy(new,s,len);
}*/
