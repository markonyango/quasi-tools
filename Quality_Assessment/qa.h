#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <inttypes.h>
#include <math.h>

typedef struct base
{
    uint32_t A;
    uint32_t T;
    uint32_t G;
    uint32_t C;
    uint32_t N;
} base;

/* This global variable holds the offset for the true Phred score */
uint8_t phred_offset;

/* This global array holds the number of sequenced bases per cycle */
uint64_t *number_of_bases_per_cycle;

/* This global variable hold the maximum sequence length */
uint16_t max_seqlength;

/* Function declarations */
void base_content(base bases[], char *read, uint16_t length);
void distribution_phred(uint64_t *dist_phred, uint64_t *dist_phred_ATGC, char *cigar, uint16_t length, char *phred_sequence);
uint8_t detect_phred(FILE *fileptr);
void boxplot(uint32_t *boxplot_bins, char *line, uint16_t length);
clock_t BeginTimer();
clock_t EndTimer(clock_t begin);
uint32_t countReads(FILE *fileptr);