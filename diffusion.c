#define _XOPEN_SOURCE 700

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

	int
main(
		int argc,
		char * argv[]
	)
{


	// Parse init file with input values
	FILE *fp = fopen("init", "r");
	const int width, height;
	fscanf(fp, "%d %d", &width, &height);

	float *init = (float*) malloc(sizeof(float)*width*height);
	for (size_t i = 0; i < width*height; i++)
		init[i] = 0.;

	int row, col;
	float temp;
	while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
		init[row * width + col] = temp;
	}

	fclose(fp);

	printf("done parsing file\n");

	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


	const int bcast_root = 0;
	int msg; int len = 1;
	int rows = 100, cols = 100;
	const int sz = rows * cols;
	// With this expression we round the quotient up.
	const int sz_row = rows / nmb_mpi_proc;


	float a = 1.33;
/*
	float *a;
	a = malloc(sz*sizeof(float));
	
	if ( mpi_rank == bcast_root ) {
		for ( int ix = 0; ix < sz; ++ix )
			a[ix] = init[ix];
	}
*/

	int from, froms[nmb_mpi_proc];
	int to, tos[nmb_mpi_proc];
	if ( mpi_rank == bcast_root ) {
		for ( int jx = 0, from = 0; jx < nmb_mpi_proc; ++jx, from += sz_row) {
			froms[jx] = from;
			tos[jx] = from + sz_row <= sz ? sz_row * (jx + 1) : rows - from;
			printf("from %d to %d\n", froms[jx], tos[jx]);
		}
	}

	if ( mpi_rank == bcast_root )
		msg = 1;


	MPI_Scatter(froms, 1, MPI_INT, &from, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Scatter(tos, 1, MPI_INT, &to, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);

	
	MPI_Bcast(&a, 1, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);

	if ( mpi_rank != bcast_root ){
		printf( " received at %d: %d\n", mpi_rank, msg );
		printf("from %d to %d received at: %d\n",from,to,mpi_rank);
}
	/*
	 * MPI_Bcast is a only one example of many advanced group-communication
	 * functions provided by MPI:
	 * MPI_Reduce, MPI_Allreduce
	 * MPI_Gather, MPI_Allgather, MPI_Gatherv, MPI_Allgatherv
	 * MPI_Scatter, MPI_Scatterv
	 * MPI_Reduce_scatter, MPI_Reduce_scatterv
	 * MPI_Alltoall, MPI_Alltoallv
	 */

	MPI_Finalize();


	return 0;
}
