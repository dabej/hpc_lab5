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


	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


	const int bcast_root = 0;
	float *a;

	if ( mpi_rank == bcast_root ) {

		// Parse init file with input values
		FILE *fp = fopen("init", "r");
		const int width, height;
		fscanf(fp, "%d %d", &width, &height);

		a = (float*) malloc(sizeof(float)*width*height);
		
		for (size_t i = 0; i < width*height; i++)
			a[i] = 0.;

		int row, col;
		float temp;
		while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
			a[row * width + col] = temp;
		}

		fclose(fp);
		printf("done parsing file\n");
	}

	int rows = 100, cols = 100;
	const int sz = rows * cols;
	// With this expression we round the quotient up.
	const int sz_row = rows / nmb_mpi_proc;

	if ( mpi_rank != bcast_root )
		a = malloc(sz*sizeof(float));

	int from, froms[nmb_mpi_proc];
	int to, tos[nmb_mpi_proc];
	if ( mpi_rank == bcast_root ) {
		for ( int jx = 0, from = 0; jx < nmb_mpi_proc; ++jx, from += sz_row) {
			froms[jx] = from;
			tos[jx] = from + sz_row <= sz ? sz_row * (jx + 1) : rows - from;
		}
	}

// send indices to each node
	MPI_Scatter(froms, 1, MPI_INT, &from, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Scatter(tos, 1, MPI_INT, &to, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);

	// bcast data array	
	MPI_Bcast(a, sz, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);

	if ( mpi_rank != bcast_root ){
		printf("from %d to %d received at: %d\n",from,to,mpi_rank);
		printf("%f received at %d\n",a[0],mpi_rank);
	}


// Computation here

		




	MPI_Finalize();


	return 0;
}
