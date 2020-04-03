/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : <1010124>
 *   Name        : <Yuxuan Liang>
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char* argv[]) {

	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	char* q2_file = NULL;
	char* q4_file = NULL;
	char* q5_file = NULL;
	double xo;
	char* q6_file = NULL;

	q2_file = argv[1];
	q4_file = argv[2];
	q5_file = argv[3];
	xo = atof(argv[4]);
	q6_file = argv[5];
	/* TODO: Add timing for each task and output running time in ms */
	struct timeval start;
	struct timeval stop;

	/* Question 2 */
	gettimeofday(&start, NULL);
	shockwave(q2_file);
	gettimeofday(&stop, NULL);
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("Question 2:  %.2f milliseconds\n", elapsed_ms);

	/* Question 4 */
	gettimeofday(&start, NULL);
	linalgbsys(q4_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("Question 4:  %.2f milliseconds\n", elapsed_ms);

	/* Question 5 */
	gettimeofday(&start, NULL);
	interp(q5_file, xo);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("Question 5:  %.2f milliseconds\n", elapsed_ms);

	/* Question 6 */
	gettimeofday(&start, NULL);
	waveeqn(q6_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("Question 6:  %.2f milliseconds\n", elapsed_ms);

	return (EXIT_SUCCESS);
}
