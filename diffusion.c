#define _XOPEN_SOURCE 700

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

int main (int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	const int width, height;
	const int bcast_root = 0;
	float *a;

	if ( mpi_rank == bcast_root ) {

		// Parse init file with input values
		FILE *fp = fopen("init", "r");
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

	const int sz = width * height;
	// With this expression we round the quotient up.
	const int rows_loc = (height - 1) / nmb_mpi_proc + 1;

	if ( mpi_rank != bcast_root )
		a = malloc(sz*sizeof(float));

	int from, froms[nmb_mpi_proc];
	int to, tos[nmb_mpi_proc];
	if ( mpi_rank == bcast_root ) {
		for ( int jx = 0, from = 0; jx < nmb_mpi_proc; ++jx, from += width*rows_loc) {
			froms[jx] = from;
			tos[jx] = from + width*rows_loc <= sz ? from + width*rows_loc : sz;
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
	for (size_t i = from; i < to; i++) {
		float value = a[i];
		float up = 0., down = 0., left = 0., right = 0.;

		if (i%width != 0)
			left = a[i-1];
			
		if ((i+1)%width != 0)
			right = a[i+1];
		
		if (i >= width)
			up = a[i - width];
		
		if (i < width*(height-1))
			down = a[i + width];

		value += d * ((up+down+left+right)/4 - value);
		c_loc[i-from] = value;
		printf("%f\n", value);
	}

	MPI_Finalize();

	return 0;
}
