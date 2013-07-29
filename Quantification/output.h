/*
 * output.h
 *
 *  Created on: 26.07.2011
 *      Author: mark
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

typedef float elem_type;

elem_type median(elem_type m[], long int n);
void output_DESEQ(uint8_t filectr, uint8_t argc, char *argv[]);
void output_EDGER(uint8_t filectr, char *argv[]); /* August 15th 2011: now obsolete */
void estimateSizeFactors(uint8_t filectr, uint16_t loc_max_libsize, uint8_t argc, char *argv[]);

#endif /* OUTPUT_H_ */
