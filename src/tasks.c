/***************************************************************************
 *
 *   File        : tasks.c
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

#define MAX_BUF_LEN 2048
#define pi 3.1415926535
#define CRITERION 1e-9
#define dx 1e-12
#define MAX_ITERATION 1000
#define INITIAL_SIZE 5
#define BETA_L_INITIAL_GUESS 30
#define BETA_U_INITIAL_GUESS 85
#define COLUMNS_OF_INFILE_Q2 1
#define COLUMNS_OF_INFILE_Q4 4
#define COLUMNS_OF_INFILE_Q5 2
#define COLUMNS_OF_INFILE_Q6 4

// define the struct for linear system in Question 4
typedef struct {
	double a;
	double b;
	double c;
	double Q;
	double x;
} linear_t;

// define the struct for data points in Question 5
typedef struct {
	double x;
	double y;
} data_t;

// define the struct for interpolation function in Question 5
typedef struct {
	double a;
	double b;
	double c;
	double d;
} interp_t;

// define the struct for the tri_diagonal matrix in Question 5
typedef struct {
	double a;
	double b;
	double c;
} tri_diagonal_t;

// output error if the file does not exist
FILE * safe_fopen(const char* path, const char* mode)
{
	FILE* fp = fopen(path, mode);
	if (fp == NULL) {
		perror("file open error.");
		exit(EXIT_FAILURE);
	}
	return fp;
}

// output error if the pointer is NULL
void* safe_malloc(size_t num_bytes)
{
	void* ptr = malloc(num_bytes);
	if (ptr == NULL) {
		printf("ERROR: malloc(%lu)\n", num_bytes);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

// output error if the new pointer is NULL
void* safe_realloc(void* ptr, size_t num_bytes)
{
	void* new_ptr = realloc(ptr, num_bytes);
	if (new_ptr == NULL) {
		printf("ERROR: realloc(%lu)\n", num_bytes);
		exit(EXIT_FAILURE);
	}
	return new_ptr;
}

// function for changing degree into radian
double degree_to_radian(double degree)
{
	return degree * pi / 180;
}

// function for changing radian into degree
double radian_to_degree(double radian)
{
	return radian / pi * 180;
}

// f(beta) in Question 2
double f_beta(double M, double theta, double beta, double gamma)
{
	return 2 * cos(beta) / sin(beta) * (pow(M, 2) * pow(sin(beta), 2) - 1) / (pow(M, 2) * (gamma + cos(2 * beta)) + 2) - tan(theta);
}

// calculate the derivative of f(beta)
double fprime_beta(double M, double theta, double beta, double gamma)
{
	return (f_beta(M, theta, beta + dx / 2, gamma) - f_beta(M, theta, beta - dx / 2, gamma)) / dx;
}

// Newton Raphson method
double newtonraph(double M, double theta, double beta, double gamma)
{
	for (int i = 0; i < MAX_ITERATION; i++)
	{
		if (fabs(f_beta(M, theta, beta, gamma)) < CRITERION)
		{
			break;
		}
		// find and return the root when iteration reaches the maximum iteration or f(beta) is smaller than the criterion
		beta = beta - f_beta(M, theta, beta, gamma) / fprime_beta(M, theta, beta, gamma);
	}
	return beta;
}

// RHS function for 1st-order upwind method in Question 6
double RHS_1U(double** func, int i, int j, int c, double deltax)
{
	if (i == 0)
		// for i = 0, use the boundary stencil delta f/delta x = fni+1-fni/delta x
		// RHS = -c * delta f/delta x
		return -c * (func[i + 1][j] - func[i][j]) / deltax;
	// first-order accurate upwind scheme: delta f/delta x = fni-fni-1/delta x for i=1,...,Nx
	// RHS = -c * delta f/delta x
	return -c * (func[i][j] - func[i - 1][j]) / deltax;
}

// RHS function for 2nd-order central method in Question 6
double RHS_2C(double** func, int i, int j, int c, double deltax, int Nx)
{
	if (i == 0)
		// for i = 0, use the boundary stencil delta f/delta x = fni+1 - fni/delta x
		// RHS = -c * delta f/delta x
		return -c * (func[i + 1][j] - func[i][j]) / deltax;
	else if (i == Nx)
		// for i = Nx, use the boundary stencil delta f/delta x = fni - fni-1/delta x
		// RHS = -c * delta f/delta x
		return -c * (func[i][j] - func[i - 1][j]) / deltax;
	// second-order accurate central scheme: delta f/delta x = fni+1 - fni-1/2 * delta x for i=1,...,Nx-1
	// RHS = -c * delta f/delta x
	return -c * (func[i + 1][j] - func[i - 1][j]) / (2 * deltax);
}

void shockwave(const char* q2_file)
{
	double M, theta, beta_l, beta_u, gamma;
	FILE* fp1 = safe_fopen(q2_file, "r");
	char* buf = (char*)safe_malloc(MAX_BUF_LEN * sizeof(char));

	// question(a)
	fgets(buf, MAX_BUF_LEN, fp1); // skip the first line of the file
	fscanf(fp1, "%lf,%lf,%lf,%lf,%lf\n", &M, &theta, &beta_l, &beta_u, &gamma); // read M from the second line
	// set the initial guess for beta_l and beta_u
	beta_l = BETA_L_INITIAL_GUESS;
	beta_u = BETA_U_INITIAL_GUESS;
	// change the angle from degrees to radians
	theta = degree_to_radian(theta);
	beta_l = degree_to_radian(beta_l);
	beta_u = degree_to_radian(beta_u);
	// find beta_l and beta_u using Newton-Raphson method
	beta_l = newtonraph(M, theta, beta_l, gamma);
	beta_u = newtonraph(M, theta, beta_u, gamma);
//	printf("(a)\nbeta_l = %.6f, beta_u = %.6f\n", radian_to_degree(beta_l), radian_to_degree(beta_u));

	// question(b)
	fgets(buf, MAX_BUF_LEN, fp1); // skip the third line
	fscanf(fp1, "%lf\n", &M); // read M from the fourth line
	beta_l = asin(1 / M);
	beta_u = pi / 2;
//	printf("(b)\n%.5f,%.5f\n", radian_to_degree(beta_l), radian_to_degree(beta_u));
	for (theta = 1; beta_l >= 0 && beta_l <= pi / 2 && beta_u >= 0 && beta_u <= pi / 2; theta++)
	{
		beta_l = newtonraph(M, degree_to_radian(theta), beta_l, gamma);
		beta_u = newtonraph(M, degree_to_radian(theta), beta_u, gamma);
		if (beta_l >= 0 && beta_l <= pi / 2 && beta_u >= 0 && beta_u <= pi / 2)
		{
		//	printf("%d,%.6f,%.6f\n", (int)theta, radian_to_degree(beta_l), radian_to_degree(beta_u));
		}
	}
	FILE* fp2 = safe_fopen("out_shock.csv", "a");
	fprintf(fp2, "M,theta,beta_lower,beta_upper\n");
	fgets(buf, MAX_BUF_LEN, fp1); // skip the fifth line
	while (fscanf(fp1, "%lf\n", &M) == COLUMNS_OF_INFILE_Q2) // keep reading M until there is no next line
	{
		// when theta = 0, beta_lower = asin(1 / M) and beta_upper = 90бу
		theta = 0;
		beta_l = asin(1 / M);
		beta_u = pi / 2;
		fprintf(fp2, "%.6f,%.0f,%.6f,%.6f\n", M, radian_to_degree(theta), radian_to_degree(beta_l), radian_to_degree(beta_u));

		// begin with theta = 1бу, stop when beta is not in the range of 0бу- 90бу
		for (theta = 1; beta_l >= 0 && beta_l <= pi / 2 && beta_u >= 0 && beta_u <= pi / 2; theta++)
		{
			beta_l = newtonraph(M, degree_to_radian(theta), beta_l, gamma);
			beta_u = newtonraph(M, degree_to_radian(theta), beta_u, gamma);
			// if beta is in the right range, print it
			if (beta_l >= 0 && beta_l <= pi / 2 && beta_u >= 0 && beta_u <= pi / 2)
			{
				fprintf(fp2, "%.6f,%d,%.6f,%.6f\n", M, (int)theta, radian_to_degree(beta_l), radian_to_degree(beta_u));
			}
		}
	}
	// close the input and output file
	fclose(fp1);
	fclose(fp2);
	// free the memory of the buffer
	free(buf);
}

void linalgbsys(const char* q4_file)
{
	int N = INITIAL_SIZE, N_read = 0;
	// allocate memories for the linear system
	linear_t* linear = (linear_t*)safe_malloc(N * sizeof(linear_t));
	FILE* fp = safe_fopen(q4_file, "r");
	char* buf = (char*)safe_malloc(MAX_BUF_LEN * sizeof(char));
	fgets(buf, MAX_BUF_LEN, fp); // skip the first line of the file
	while (fscanf(fp, "%lf,%lf,%lf,%lf\n", &linear[N_read].a, &linear[N_read].b, &linear[N_read].c, &linear[N_read].Q) == COLUMNS_OF_INFILE_Q4)
	{
		N_read++;
		// when the lines reach the maximum size of the array, double the size
		if (N == N_read)
		{
			N = 2 * N;
			linear = (linear_t*)safe_realloc(linear, N * sizeof(linear_t));
		}
	}
	fclose(fp);

	// realloc the memory to the final size (free the extra memory)
	N = N_read;
	linear = (linear_t*)safe_realloc(linear, N * sizeof(linear_t));
	
	// rewrite the matrix, ai* and Qi* need to be calculated for i = 2,3,...,N, so the index of the array should start with 1
	for (int i = 1; i < N; i++)
	{
		// ai* = ai - cibi-1/ai-1*
		linear[i].a = linear[i].a - linear[i].c * linear[i - 1].b / linear[i - 1].a;
		// Qi* = Qi - ciQi-1*/ai-1*
		linear[i].Q = linear[i].Q - linear[i].c * linear[i - 1].Q / linear[i - 1].a;
	}
	// then calculate the solution to the original tri-diagonal matrix
	// xi = Qi*/ai* for i = N, i.e. index of array = N - 1
	linear[N - 1].x = linear[N - 1].Q / linear[N - 1].a;

	for (int i = N - 2; i >= 0; i--)
	{	
		// xi = (Qi* - bixi+1)/ai* for i = N-1, N-2,...,1, i.e. index starts from N-2
		linear[i].x = (linear[i].Q - linear[i].b * linear[i + 1].x) / linear[i].a;
	}

	fp = safe_fopen("out_linalsys.csv", "a");
	fprintf(fp, "x\n");
	for (int i = 0; i < N; i++)
	{
		// output x1, x2, x3... xN
		fprintf(fp, "%.6f\n", linear[i].x);
	}
	fclose(fp);
	free(buf);
	free(linear);
}

void interp(const char* q5_file, const double xo)
{
	int N = INITIAL_SIZE, N_read = 0;
	// allocate memories for the data points
	data_t* data = (data_t*)safe_malloc(N * sizeof(data_t));
	FILE* fp = safe_fopen(q5_file, "r");
	char* buf = (char*)safe_malloc(MAX_BUF_LEN * sizeof(char));
	fgets(buf, MAX_BUF_LEN, fp); // skip the first line of the file
	while (fscanf(fp, "%lf,%lf\n", &data[N_read].x, &data[N_read].y) == COLUMNS_OF_INFILE_Q5)
	{
		N_read++;
		// when the lines reach the maximum size of the array, double the size
		if (N == N_read)
		{
			N = 2 * N;
			data = (data_t*)safe_realloc(data, N * sizeof(data_t));
		}
	}
	fclose(fp);
	// realloc the memory to the final size (free the extra memory)
	N = N_read;
	data = (data_t*)safe_realloc(data, N * sizeof(data_t));

	// allocate memories for coefficients of the interpolated function
	interp_t* interp = (interp_t*)safe_malloc(N * sizeof(interp_t));
	memset(interp, 0, N * sizeof(interp_t));

	// obtain the hi's
	int n = N - 1; // n is number of intervals, which equals to N - 1

	// allocate memories for hi's
	double* h = (double*)safe_malloc(n * sizeof(double));
	memset(h, 0, n * sizeof(double));

	for (int i = 0; i < n; i++)
	{
		// calculate hi = xi+1 - xi
		h[i] = data[i + 1].x - data[i].x;
	}

	// obtain ai's using ai = f(xi)
	for (int i = 0; i < N; i++)
	{
		interp[i].a = data[i].y;
	}

	// solve [A]{X} = {C} to obtain ci's, where [A] is a tri_diagonal matrix

	// allocate memories for three vectors (a, b, c) in the tri_diagonal matrix [A]
	tri_diagonal_t* tri_diagonal = (tri_diagonal_t*)safe_malloc(N * sizeof(tri_diagonal_t));
	memset(tri_diagonal, 0, N * sizeof(tri_diagonal_t));

	// set values of three vectors
	// set vector a
	// a1 and aN is equal to 1
	tri_diagonal[0].a = 1;
	tri_diagonal[N - 1].a = 1;

	// ai = 2(hi-1 + hi) for i = 2,3,...,N-1
	for (int i = 1; i < n; i++)
	{
		tri_diagonal[i].a = 2 * (h[i - 1] + h[i]);
	}

	// set vector b, bi = 0, h1, h2...hn-1, 0
	tri_diagonal[0].b = 0;
	tri_diagonal[N - 1].b = 0;
	for (int i = 1; i < n; i++)
	{
		tri_diagonal[i].b = h[i];
	}

	// set vector c, ci = 0, h0, h1, ... hn-2, 0
	tri_diagonal[0].c = 0;
	tri_diagonal[N - 1].c = 0;
	for (int i = 1; i < n; i++)
	{
		tri_diagonal[i].c = h[i - 1];
	}

	// allocate memories for {C}
	double* vectorC = (double*)safe_malloc(N * sizeof(double));
	memset(vectorC, 0, N * sizeof(double));
	// set values of {C}, the first and the last element of {C} is 0
	vectorC[0] = 0;
	vectorC[n] = 0;

	for (int i = 1; i < n; i++)
	{
		// Ci = 3/hi*(ai+1 - ai) + 3/hi-1*(ai-1 - ai)
		vectorC[i] = 3 / h[i] * (interp[i + 1].a - interp[i].a) + 3 / h[i - 1] * (interp[i - 1].a - interp[i].a);
	}

	// use Thomas algorithm to solve the tri-diagonal system
	// rewrite the matrix
	for (int i = 1; i < N; i++)
	{
		// ai* = ai - cibi-1/ai-1*
		tri_diagonal[i].a = tri_diagonal[i].a - tri_diagonal[i].c * tri_diagonal[i - 1].b / tri_diagonal[i - 1].a;
		// Ci* = Ci - ciCi-1*/ai-1*
		vectorC[i] = vectorC[i] - tri_diagonal[i].c * vectorC[i - 1] / tri_diagonal[i - 1].a;
	}
	// find the solution of {X}, obtain the ci's of the interpolated function
	// ci = Ci*/ai* for i = N
	interp[n].c = vectorC[n] / tri_diagonal[n].a;
	for (int i = n - 1; i >= 0; i--)
	{
		// ci = (Ci*-bici+1)/ai*
		interp[i].c = (vectorC[i] - tri_diagonal[i].b * interp[i + 1].c) / tri_diagonal[i].a;
	}
	// obtain the bi's, bi = 1/hi(ai+1 - ai) - hi/3(2ci + ci+1)
	interp[n].b = 0;
	for (int i = 0; i < n; i++)
	{
		interp[i].b = 1 / h[i] * (interp[i + 1].a - interp[i].a) - h[i] / 3 * (2 * interp[i].c + interp[i + 1].c);
	}

	// obtain the di's, di = (ci+1-ci)/3hi
	interp[n].d = 0;
	for (int i = 0; i < n; i++)
	{
		interp[i].d = (interp[i + 1].c - interp[i].c) / (3 * h[i]);
	}

	fp = safe_fopen("out_interp.csv", "a");
	fprintf(fp, "xo,f(xo)\n");

	for (int i = 0; i < n; i++)
	{
		// if xo is between two data points, use the interpolated function of that interval to calculate the interpolated value
		if ((xo <= data[i].x && xo >= data[i + 1].x) || (xo >= data[i].x && xo <= data[i + 1].x))
		{
			fprintf(fp, "%.6f,%.6f\n", xo, interp[i].a + interp[i].b * (xo - data[i].x) + interp[i].c * pow((xo - data[i].x), 2) + interp[i].d * pow((xo - data[i].x), 3));
		}
	}
	fclose(fp);

	// free all the allocated memories
	free(buf);
	free(data);
	free(interp);
	free(h);
	free(tri_diagonal);
	free(vectorC);
}


void waveeqn(const char* q6_file)
{
	int Nx, out_iter;
	double c, CFL;
	FILE* fp = safe_fopen(q6_file, "r");
	char* buf = (char*)safe_malloc(MAX_BUF_LEN * sizeof(char));
	fgets(buf, MAX_BUF_LEN, fp); // skip the first line of the file
	// read values of c, Nx, CFL and out_iter
	while (fscanf(fp, "%lf,%d,%lf,%d\n", &c, &Nx, &CFL, &out_iter) == COLUMNS_OF_INFILE_Q6);
	fclose(fp);

	// calculate delta x and delta t
	double deltax = 1.0 / Nx;
	double deltat = CFL * deltax / c;

	// allocate memories for fni and fn+0.5i, no need to allocate for fn+1i because it is just the next element of fni
	// allocate fni for 1st and 2nd order method
	// rows of the 2D array correspond to different x
	double** fn_1U = (double**)safe_malloc((Nx + 1) * sizeof(double*));
	double** fn_2C = (double**)safe_malloc((Nx + 1) * sizeof(double*));
	for (int i = 0; i < (Nx + 1); i++)
	{
		// columns of the 2D array correspond to different time level
		// initialize the arrays
		fn_1U[i] = (double*)safe_malloc((out_iter + 1) * sizeof(double));
		memset(fn_1U[i], 0, (out_iter + 1) * sizeof(double));
		fn_2C[i] = (double*)safe_malloc((out_iter + 1) * sizeof(double));
		memset(fn_2C[i], 0, (out_iter + 1) * sizeof(double));

		// set initial condition for fni
		if (i * deltax >= 0 && i * deltax < 0.125)
		{
			fn_1U[i][0] = 0;
			fn_2C[i][0] = 0;
		}
		if (i * deltax >= 0.125 && i * deltax <= 0.375)
		{
			fn_1U[i][0] = 0.5 * (1 - cos(8.0 * pi * (i * deltax - 0.125)));
			fn_2C[i][0] = 0.5 * (1 - cos(8.0 * pi * (i * deltax - 0.125)));
		}
		if (i * deltax > 0.375 && i * deltax <= 1)
		{
			fn_1U[i][0] = 0;
			fn_2C[i][0] = 0;
		}
	}

	// allocate memories for fn+0.5i for both method
	double** fnplushalf_1U = (double**)safe_malloc((Nx + 1) * sizeof(double*));
	double** fnplushalf_2C = (double**)safe_malloc((Nx + 1) * sizeof(double*));
	for (int i = 0; i < (Nx + 1); i++)
	{
		// initialize the fn+0.5i arrays for both method
		fnplushalf_1U[i] = (double*)safe_malloc((out_iter + 1) * sizeof(double));
		memset(fnplushalf_1U[i], 0, (out_iter + 1) * sizeof(double));
		fnplushalf_2C[i] = (double*)safe_malloc((out_iter + 1) * sizeof(double));
		memset(fnplushalf_2C[i], 0, (out_iter + 1) * sizeof(double));
	}

	// loop over the number of timesteps, each with delta t
	for (int j = 0; j < out_iter; j++)
	{
		// loop over number of x, each with delta x, calculate fn+0.5i for both method, i.e. the function values at intermediate level n+0.5
		for (int i = 0; i < Nx + 1; i++)
		{
			fnplushalf_1U[i][j] = fn_1U[i][j] + deltat * RHS_1U(fn_1U, i, j, c, deltax);
			fnplushalf_2C[i][j] = fn_2C[i][j] + deltat * RHS_2C(fn_2C, i, j, c, deltax, Nx);
		}
		// loop over number of x, each with delta x, calculate fn+1i for both method, which is the next element of fni
		for (int i = 0; i < Nx + 1; i++)
		{
			fn_1U[i][j + 1] = fn_1U[i][j] + deltat / 2 * (RHS_1U(fn_1U, i, j, c, deltax) + RHS_1U(fnplushalf_1U, i, j, c, deltax));
			fn_2C[i][j + 1] = fn_2C[i][j] + deltat / 2 * (RHS_2C(fn_2C, i, j, c, deltax, Nx) + RHS_2C(fnplushalf_2C, i, j, c, deltax, Nx));
		}
	}

	fp = safe_fopen("out_waveeqn_1U.csv", "a");
	fprintf(fp, "x,f(x)\n");
	for (int i = 0; i < Nx + 1; i++)
	{
		// output the solution after the kth time loop executes for 1st order method, which is the last element of each rows
		fprintf(fp, "%.6f,%.6f\n", i * deltax, fn_1U[i][out_iter]);
	}
	fclose(fp);

	fp = safe_fopen("out_waveeqn_2C.csv", "a");
	fprintf(fp, "x,f(x)\n");
	for (int i = 0; i < Nx + 1; i++)
	{
		// output the solution after the kth time loop executes for 2nd order method, which is the last element of each rows
		fprintf(fp, "%.6f,%.6f\n", i * deltax, fn_2C[i][out_iter]);
	}

	fclose(fp);

	// free all the allocated memories
	free(buf);
	for (int i = 0; i < (Nx + 1); i++)
	{
		free(fn_1U[i]);
		free(fnplushalf_1U[i]);
	}
	free(fn_1U);
	free(fnplushalf_1U);

	for (int i = 0; i < (Nx + 1); i++)
	{
		free(fn_2C[i]);
		free(fnplushalf_2C[i]);
	}
	free(fn_2C);
	free(fnplushalf_2C);
}
