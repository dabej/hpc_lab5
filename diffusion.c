PEN_SOURCE 700

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

int main (int argc, char *argv[]) {
	// Argpars
	int opt, n;
	float d;
	while ((opt = getopt(argc, argv, "n:d:")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				break;
			case 'd':
				d = atof(optarg);
				break;
		}
	}
	
	// Parse init file with input values
	FILE *fp = fopen("init", "r");
	const int width, height;
	fscanf(fp, "%d %d", &width, &height);

	float *a = (float*) malloc(sizeof(float)*width*height);
	for (size_t i = 0; i < width*height; i++)
		a[i] = 0.;

	int row, col;
	float temp;
	while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
		a[row * width + col] = temp;
	}
	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


	int scatter_root = 1;
	int reduce_root = 0;

	const int sz = 10000000;
	// With this expression we round the quotient up.
	const int sz_loc = (sz - 1) / nmb_mpi_proc + 1;

	// Allocate local arrays and prepare data on the scatter root.
	double *a;
	double *a_loc = (double*) malloc(sz_loc*sizeof(double));
	if ( mpi_rank == scatter_root ) {
		a = malloc(sz*sizeof(double));
		for ( int ix = 0; ix < sz; ++ix )
			a[ix] = ix;
	}

	// The function scatterv requires arrays that indicate the positions and
	// length in the scattered vector associated with each of process.
	//
	// For the purpose of MPI they only need to be initialized at scatter_root.
	// To illustrate another use of scatter, we communicate them from there.
	int pos, poss[nmb_mpi_proc];
	int len, lens[nmb_mpi_proc];
	if ( mpi_rank == scatter_root )
		for ( int jx = 0, pos = 0; jx < nmb_mpi_proc; ++jx, pos += sz_loc) {
			poss[jx] = pos;
			lens[jx] = pos + sz_loc <= sz ? sz_loc : sz - pos;
			printf("len pos %d %d\n", lens[jx], poss[jx]);
		}

	// The sendcount argument determines the number of elements sent to EACH
	// process.
	MPI_Scatter(poss, 1, MPI_INT, &pos, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
	MPI_Scatter(lens, 1, MPI_INT, &len, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);

	MPI_Scatterv(a, lens, poss, MPI_DOUBLE, a_loc, sz_loc, MPI_DOUBLE,
			scatter_root, MPI_COMM_WORLD);

	if ( mpi_rank == scatter_root )
		free(a);


	// Compute the maximum and its location in the local part of a.
	double max = -1.;
	int loc;
	for ( int ix = 0; ix < len; ++ix )
		if ( a_loc[ix] > max ) {
			max = a_loc[ix];
			loc = ix + pos;
		}

	free(a_loc);


	// Now reduce the maxima and their locations globally.
	struct {
		double v;
		int    l;
	} maxloc, maxloc_glob;

	maxloc.v = max;
	maxloc.l = loc;
	MPI_Reduce(&maxloc, &maxloc_glob, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
			reduce_root, MPI_COMM_WORLD);

	if ( mpi_rank == reduce_root )
		printf( "Global maximum %f at location %d\n", maxloc_glob.v, maxloc_glob.l );


	MPI_Finalize();

	float *c = malloc(width*height* sizeof(float));
	// average temp	
	float sum = 0.;
	for (size_t jx=0; jx<height; ++jx)
		for (size_t ix=0; ix<width; ++ix)
			sum +=  c[jx*width+ ix];
	float avg_temp = sum/(height*width);
	printf("average: %E\n", avg_temp);


	// the absolute difference of each temperature and the average
	sum = 0.;
	for (size_t jx=0; jx<height; ++jx)
		for (size_t ix=0; ix<width; ++ix)
			sum += fabs(c[jx*width+ ix] - avg_temp);
	avg_temp = sum/(width*height);
	printf("average absolute difference: %E\n", avg_temp);

	return 0;
}

