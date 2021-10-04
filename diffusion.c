#define _XOPEN_SOURCE 700

#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

int main (int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	void computation(int width, int height, int from, int to, float d, int mpi_rank, float **input, float **output_loc);
	void parseFile(float **input,int *width, int *height);

	int width = -1, height = -1;
	float d;
	int n;
	const int bcast_root = 0;
	float *input;


	// Parse arguments and init-file
	if ( mpi_rank == bcast_root ) {

		// Argpars
		int opt;
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

		//parseFile(&a,&width,&height);

		// Parse init file with input values
		FILE *fp = fopen("init", "r");
		fscanf(fp, "%d %d", &width, &height);

		input =  malloc(sizeof(float)*width*height);
		for(int i = 0; i< width*height; i++)
			input[i] = 0.;

		int row, col;
		float temp;
		while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
			input[row * width + col] = temp;
		}

		fclose(fp);
		printf("done parsing file\n");
	}

	MPI_Barrier(MPI_COMM_WORLD); // Maybe needed for larger init files..


	// Broadcast params to each node	
	MPI_Bcast(&width, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Bcast(&height, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Bcast(&d, 1, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);


	const int sz = width * height;
	const int rows_loc = (height - 1) / nmb_mpi_proc + 1;

	// Allocate input for all other nodes
	if ( mpi_rank != bcast_root )
		input = malloc(sz*sizeof(float));


	// Compute the work for each node
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


	// initialize array that will hold the result of the computations in each node
	float *output_loc = malloc((to - from) * sizeof(float));

	float *output = NULL;
	if (mpi_rank == bcast_root){
		output = malloc(sz*sizeof(float));
		for (int i = 0; i<sz; i++)
			output[i] = 0.;
	}


	// n iterations of diffusion calculations
	for (size_t ix = 0; ix < n ; ix++) {

		// Broadcast input array	
		MPI_Bcast(input, sz, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);

		// Compute results
		computation( width, height, from, to, d, mpi_rank, &input, &output_loc);

		// Gather results into output-array	
		MPI_Gather(output_loc,to-from,MPI_FLOAT,output,to-from,MPI_FLOAT,bcast_root, MPI_COMM_WORLD);

		// set output array as input for the next iter
		if (mpi_rank == bcast_root) {
			float *temp = input;
			input = output;
			output = temp;
		}

		MPI_Barrier(MPI_COMM_WORLD); // Maybe needed for larger init files..
	}

	free(output);
	free(output_loc);

	if (mpi_rank == bcast_root) {
		// average temp	
		float sum = 0.;
		for (size_t jx=0; jx<height; ++jx){
			for (size_t ix=0; ix<width; ++ix){
				//printf("%.0f ",input[jx*width+ ix]);
				sum +=  input[jx*width+ ix];
			}
			//printf("\n");
		}
		float avg_temp = sum/(height*width);
		printf("average: %E\n", avg_temp);
	}

	free(input);

	MPI_Finalize();

	return 0;
}


void computation(int width, int height, int from, int to, float d, int mpi_rank, float **input, float **output_loc){
	for (size_t i = from; i < to; i++) {
		float value = (*input)[i];
		float up = 0., down = 0., left = 0., right = 0.;

		if (i%width != 0)
			left = (*input)[i-1];

		if ((i+1)%width != 0)
			right = (*input)[i+1];

		if (i >= width)
			up = (*input)[i - width];

		if (i < width*(height-1))
			down = (*input)[i + width];

		value += d * ((up+down+left+right)/4 - value);
		//	printf("u,d,l and r are %f,%f,%f,%f. val=%f, at n= %d\n",up,down,left,right,value,mpi_rank);
		(*output_loc)[i-from] = value;
		//printf("res is %f\n",(*output_loc)[i-from]);
	}
}
/*
   for (size_t i = from; i < to; i++) {
   float value = input[i];
   float up = 0., down = 0., left = 0., right = 0.;

   if (i%width != 0)
   left = input[i-1];

   if ((i+1)%width != 0)
   right = input[i+1];

   if (i >= width)
   up = input[i - width];

   if (i < width*(height-1))
   down = input[i + width];

   value += d * ((up+down+left+right)/4 - value);
//	printf("u,d,l and r are %f,%f,%f,%f. val=%f, at n= %d\n",up,down,left,right,value,mpi_rank);
output_loc[i-from] = value;
}
*/


void parseFile(float **input,int *width, int *height){

	FILE *fp = fopen("init", "r");
	fscanf(fp, "%d %d", width, height);

	*input =  malloc(sizeof(float)**width**height);
	for(int i = 0; i< *width**height; i++)
		(*input)[i] = 0.;

	int row, col;
	float temp;
	while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
		(*input)[row * *width + col] = temp;
	}

	fclose(fp);
	printf("done parsing file\n");
}

