/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	char* q2_file = NULL;
	char* q4_file = NULL;
	char* q5_file = NULL;
	double xo;
	char* q6_file = NULL;

	/* TODO: Add timing for each task and output running time in ms */
    
	/* Question 2 */
	shockwave(q2_file);
	
	/* Question 4 */
	linalgbsys(q4_file);
	
	/* Question 5 */
	interp(q5_file,xo);
	
	/* Question 6 */
	waveeqn(q6_file);
    
	return (EXIT_SUCCESS);
}
