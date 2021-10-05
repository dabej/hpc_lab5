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

	struct timespec bench_start;
	struct timespec bench_stop;
	double bench_diff;

	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

<<<<<<< HEAD
	void printarray(float *arr,int w,int h);

	// init vars
=======
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c
	int width, height, n;
	float d;
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

		// Parse init file with input values
		FILE *fp = fopen("init", "r");
		fscanf(fp, "%d %d", &width, &height);
		width += 2;
		height += 2;

<<<<<<< HEAD
		// fill input array with zeroes
		input =  malloc(sizeof(float)*(width+2)*(height+2));
		for(int i = 0; i < (width+2)*(height+2); i++)
=======
		input =  malloc(sizeof(float)*width*height);
		for(int i = 0; i < width*height; i++)
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c
			input[i] = 0.;

		int row, col;
		float temp;
		while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
<<<<<<< HEAD
			input[width+3 + (width+2)*row + col] = temp;
=======
			input[(row+1) * width + (col+1)] = temp;
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c
		}

		fclose(fp);
	}

	MPI_Barrier(MPI_COMM_WORLD); // Maybe needed for larger init files..

	// Broadcast params to each node	
	MPI_Bcast(&width, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Bcast(&height, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Bcast(&d, 1, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);

	const int sz = width * height;
	const int sz_padded = (width+2)*(height+2);

	// Allocate input for all other nodes
	if ( mpi_rank != bcast_root )
<<<<<<< HEAD
		input = malloc(sz_padded*sizeof(float));
=======
		input = malloc(sz*sizeof(float));

	int rows_loc = (height - 3) / nmb_mpi_proc + 1;
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c

	// Compute the work for each node
	int from, froms[nmb_mpi_proc];
	int to, tos[nmb_mpi_proc];
	if ( mpi_rank == bcast_root ) {
<<<<<<< HEAD
		int rows_loc = (height - 1) / nmb_mpi_proc + 1;
		for (int jx = 0, from = 0; jx < nmb_mpi_proc; ++jx, from += rows_loc) {
			froms[jx] = from;
			tos[jx] = from + rows_loc <= height ? from + rows_loc : height;
			//printf("node %d has indices %d to %d\n",jx,from,tos[jx]);	
=======
		for ( int jx = 0, from = width+1; jx < nmb_mpi_proc; ++jx, from += width*rows_loc) {
			froms[jx] = from;
			tos[jx] = from + width*rows_loc <= sz-width-1 ? from + width*rows_loc : sz-width-1;
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c
		}
	}

	// send indices to each node
	MPI_Scatter(froms, 1, MPI_INT, &from, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);
	MPI_Scatter(tos, 1, MPI_INT, &to, 1, MPI_INT, bcast_root, MPI_COMM_WORLD);

	// initialize array that will hold the result of the computations in each node
<<<<<<< HEAD
	float *output_loc = malloc((to - from) * width * sizeof(float));
	for (int i = 0; i<(to-from)*width; i++)
		output_loc[i] = 0.;
=======
	float *output_loc = malloc((to - from - rows_loc*2 + 2) * sizeof(float));
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c

	float *output = NULL;
	if (mpi_rank == bcast_root){
		output = malloc((width-2)*(height-2)*sizeof(float));
		for (int i = 0; i < (width-2)*(height-2); i++)
			output[i] = 0.;
	}

<<<<<<< HEAD
=======
	if ( mpi_rank == bcast_root )
		timespec_get(&bench_start,TIME_UTC);
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c

	// n iterations of diffusion calculations
	for (size_t ix = 0; ix < n ; ix++) {
		// Broadcast input array	
		MPI_Bcast(input, sz_padded, MPI_FLOAT, bcast_root, MPI_COMM_WORLD);

	if ( mpi_rank == bcast_root )
		timespec_get(&bench_start,TIME_UTC);	
		// Compute results
<<<<<<< HEAD
		for (size_t row = from; row < to; row++) {
			for (size_t col = 0; col < width; col++) {
				float value = input[width+3 + row*(width+2) + col];
				float up, down, left, right;

				left = input[width+3 + row*(width+2) + col-1];
				right = input[width+3 + row*(width+2) + col+1];
				up = input[width+3 + (row-1)*(width+2) + col];
				down = input[width+3 + (row+1)*(width+2) + col];

				value += d * ((up+down+left+right)/4 - value);
				//printf("u,d,l and r are %f,%f,%f,%f. val=%f, at n= %d\n",up,down,left,right,value,mpi_rank);
				output_loc[(row-from)*width + col] = value;
			}
		}
		
	if ( mpi_rank == bcast_root )
		timespec_get(&bench_stop,TIME_UTC);	
		// Gather results into output-array	
		MPI_Gather(output_loc,width*(to-from),MPI_FLOAT,output,width*(to-from),MPI_FLOAT,bcast_root, MPI_COMM_WORLD);

		// set output array as input for the next iter
		if (mpi_rank == bcast_root) {
			for(int i=0;i<height;i++)
				for(int j=0;j<width;j++)
					input[(i*(width+2)+(width+3)+j)]=output[i*width+j];
=======
		for (size_t j = (from-1)/width; j < (to+2)/width; j++)
			for (size_t i = 1; i < width-1; i++) {
				float value = input[i + j*width];
				float up, down, left, right;
				//printf("from: %d, to: %d, i: %d, %f\n");
				left = input[i + j*width - 1];
				right = input[i + j*width + 1];
				up = input[i + (j-1)*width];
				down = input[i + (j+1)*width];

				value += d * ((up+down+left+right)/4 - value);
				//	printf("u,d,l and r are %f,%f,%f,%f. val=%f, at n= %d\n",
				//	up,down,left,right,value,mpi_rank);
				output_loc[(width-2)*(j-1)+i-1] = value;
			}

		// Gather results into output-array	
		MPI_Gather(output_loc, (to - from - rows_loc*2 + 2), MPI_FLOAT, output, 
			(to - from - rows_loc*2 + 2), MPI_FLOAT, bcast_root, MPI_COMM_WORLD);

		// set output array as input for the next iter
		if (mpi_rank == bcast_root) {
			for (size_t j = 1; j < height-1; j++)
				for (size_t i = 1; i < width-1; i++)
					input[i + j*width] = output[(i-1) + (j-1)*(width-2)];
>>>>>>> cc5f18ba713efa5f8142873a8b973cf5c6dc746c
		}
	}


	free(output_loc);

	if (mpi_rank == bcast_root) {
		// average temp	
		float sum = 0.;
		for (size_t jx=0; jx<height; ++jx){
			for (size_t ix=0; ix<width; ++ix){
				//				printf("%.0f ",output[jx*width+ ix]);
				sum +=  output[jx*width+ ix];
			}
			//			printf("\n");
		}
		float avg_temp = sum/(height*width);
		printf("average: %E\n", avg_temp);


		// the absolute difference of each temperature and the average
		sum = 0.;
		for (size_t jx=0; jx<height; ++jx)
			for (size_t ix=0; ix<width; ++ix)
				sum += fabs(output[jx*width+ ix] - avg_temp);
		avg_temp = sum/(width*height);
		printf("average absolute difference: %E\n", avg_temp);
	}

	free(input);
	free(output);

	if ( mpi_rank == bcast_root ){
		bench_diff = difftime(bench_stop.tv_sec, bench_start.tv_sec) * 1000.
			+ (bench_stop.tv_nsec - bench_start.tv_nsec) / 1000000.;
		printf("benchmark time: %.2fms\n",bench_diff/n);
	}

	MPI_Finalize();

	return 0;
}
void pad(float *out,float *in,int w, int h)
{
	int i,j;
	for(i=0;i<h+2;i++)
		for(j=0;j<w+2;j++)
			*(in+(i*(w+2)+(w+3)+j))=*(out+i*w+j);
}


void printarray(float *arr,int w,int h){
	for (int i = 0; i < h; i++){
		for (int j=0;j<w;j++){
			printf("%.1f ",arr[i*w + j]);
		}
		printf("\n");
	}
}
