#define _XOPEN_SOURCE 700

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

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
	
	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	int scatter_root = 1;
	int reduce_root = 0;
	
	const int width, height;

	// Allocate local arrays and prepare data on the scatter root.
	float *a;
	float *c;
	if ( mpi_rank == scatter_root ) {
		a = malloc(sizeof(float)*width*height);
		c = malloc(sizeof(float)*width*height);
		FILE *fp = fopen("init", "r");
		fscanf(fp, "%d %d", &width, &height);

		for (size_t i = 0; i < width*height; i++)
			a[i] = 0.;

		int row, col;
		float temp;
		while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
			a[row * width + col] = temp;
		}

		fclose(fp);
	}
	
	MPI_Bcast(a, width*height, MPI_FLOAT, scatter_root, MPI_COMM_WORLD);
	MPI_Bcast(c, width*height, MPI_FLOAT, scatter_root, MPI_COMM_WORLD);

	const int height_loc = (height -1) / nmb_mpi_proc + 1;

	int pos, poss[nmb_mpi_proc];
	int len, lens[nmb_mpi_proc];
	if ( mpi_rank == scatter_root )
		for ( int jx = 0, pos = 0; jx < nmb_mpi_proc; ++jx, pos += width*height_loc) {
			poss[jx] = pos;
			lens[jx] = pos + width*height_loc <= width*height ? 
						width*height_loc : width*height - pos;
			printf("len pos %d %d\n", lens[jx], poss[jx]);
		}

	// The sendcount argument determines the number of elements sent to EACH
	// process.
	MPI_Scatter(poss, 1, MPI_INT, &pos, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);
	MPI_Scatter(lens, 1, MPI_INT, &len, 1, MPI_INT, scatter_root, MPI_COMM_WORLD);

	if ( mpi_rank == scatter_root )
		free(a);

	float *c_loc = malloc(sizeof(float)*width*height_loc);
	for (size_t i = pos; i < len; i++) {
		float value = a[i];
		float up, down, left, right;

		if (i%width == 0)
			left = 0.;
		else 
			left = a[i-1];
			
		if ((i+1)%width == 0)
			right = 0.;
		else
			right = a[i+1];
		
		if (i < width)
			up = 0.;
		else
			up = a[i - width];
		
		if (i >= width*(height-1))
			down = 0.;
		else
			down = a[i + width];

		value += d * ((up+down+left+right)/4 - value);
		c_loc[i-pos] = value;
		printf("%f\n", value);
	}

	free(c_loc);
/*
	// average temp	
	float sum = 0.;
	for (size_t i = 0; i < width*height; i++)
		sum +=  c[i];
	float avg_temp = sum/(height*width);
	printf("average: %E\n", avg_temp);


	// the absolute difference of each temperature and the average
	sum = 0.;
	for (size_t i = 0; i < width*height; i++)
		sum += fabs(c[i] - avg_temp);
	avg_temp = sum/(width*height);
	printf("average absolute difference: %E\n", avg_temp);

	free(a);
	free(c);
*/
	return 0;
}

