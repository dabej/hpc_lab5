#define _XOPEN_SOURCE 700

#include <math.h>
// #include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

float**
malloc_shifted_matrix(
		int nrs,
		int ncs
		)
{
	float * mentries = malloc(sizeof(float) * (nrs+2) * (ncs+2));
	memset(mentries, sizeof(float)*(nrs+2)*(ncs+2), 0);
	float ** m = malloc(sizeof(float*) * (nrs+2));
	for ( int ix = 0, jx = 1; ix < nrs+2; ++ix, jx += ncs+2 )
		m[ix] = mentries + jx;
	return m+1;
}

void
free_shifted_matrix(
	float ** m
		)
{
	free(m[-1]-1);
	free(m);
}

int main(int argc, char * argv[]) {
//	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
//	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
//	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	int width, height, n;
	float d;
//	float *input;

	float ** buff0;
	float ** buff1;

	// Parse arguments and init-file
	// if (mpi_rank == 0) {

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

		FILE *fp = fopen("init", "r");
		fscanf(fp, "%d %d", &width, &height);

	        buff0 = malloc_shifted_matrix(height, width);
	        buff1 = malloc_shifted_matrix(height, width);

		width += 2;
		height += 2;
		// input =  malloc(sizeof(float)*width*height);

		int row, col;
		float temp;
		while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
			buff0[row][col] = temp;
			// input[(row+1) * width + (col+1)] = temp;
		}

		fclose(fp);
	// }

//	MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&d, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate input for all other nodes
//	if (mpi_rank != 0) {
//		input = malloc(width*height*sizeof(float));
//		memset(input, sizeof(float)*width*height, 0);
//	}

//	float *send_buffer = malloc(width*sizeof(float));
//	float *recv_buffer = malloc(width*sizeof(float));
//	int rows = (height-3) / nmb_mpi_proc + 1;
//	float *working_buffer = malloc(width*height*sizeof(float));

//	int above = mpi_rank - 1;
//	int below = mpi_rank + 1;
//	MPI_Status status;
//	int start_row = mpi_rank*rows + 1;
//	int end_row = (mpi_rank+1) * rows + 1;
//	if (end_row >= height)
//		end_row = height - 1;
	size_t start_row = 1, end_row = height-1;
	//float *up_row, *down_row;
	float c_0 = 1.f-d;
	float c_1 = 0.25f*d;
	for (size_t iter = 0; iter < n; iter++) {
		for ( int ix = 0; ix < height - 2; ++ix ) {
		// for ( size_t rx = start_row*width; rx < end_row*width; rx += width ) {
		// for (size_t row = start_row; row < end_row; row++) {
			//up_row = &input[(row-1)*width];
			//down_row = &input[(row+1)*width];
			for ( int jx = 0; jx < width - 2; ++jx ) {
			// for (size_t col = 1; col < width-1; col++) {
				// size_t c_x = row*width + col;
				float value = buff0[ix][jx];
				float up    = buff0[ix-1][jx];
				float down  = buff0[ix+1][jx];
				float left  = buff0[ix][jx-1];
				float right = buff0[ix][jx+1];
				buff1[ix][jx] = c_0*value + c_1*(up + down + left + right);
			}
		}

		float **temp = buff0;
		buff0 = buff1;
		buff1 = temp;

//		if (start_row != 1) {
//			for (size_t col = 1; col < width-1; col++)
//				send_buffer[col] = working_buffer[start_row*width + col];
//			MPI_Sendrecv(send_buffer, width, MPI_FLOAT, above, iter,
//						 recv_buffer, width, MPI_FLOAT, above, iter,
//						 MPI_COMM_WORLD, &status);	
//			for (size_t col = 1; col < width-1; col++)
//				input[(start_row-1)*width + col] = recv_buffer[col];
//		}
//
//		if (end_row != height-1) {
//			for (size_t col = 1; col < width-1; col++)
//				send_buffer[col] = working_buffer[(end_row-1)*width + col];
//			MPI_Sendrecv(send_buffer, width, MPI_FLOAT, below, iter,
//						 recv_buffer, width, MPI_FLOAT, below, iter,
//						 MPI_COMM_WORLD, &status);	
//			for (size_t col = 1; col < width-1; col++)
//				input[end_row*width + col] = recv_buffer[col];
//		}
	}

//	if (mpi_rank != 0)
//		MPI_Send(input, width*height, MPI_FLOAT, 0, mpi_rank, MPI_COMM_WORLD);
//	else {
//		for (size_t proc = 1; proc < nmb_mpi_proc; proc++) {
//			MPI_Recv(working_buffer, width*height, MPI_FLOAT, proc, proc,
//					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//			int start_row = proc*rows + 1;
//			int end_row = (proc+1) * rows + 1;
//			if (end_row >= height)
//				end_row = height - 1;
//			for (size_t row = start_row; row < end_row; row++)
//				for (size_t col = 1; col < width-1; col++)
//					input[row*width + col] = working_buffer[row*width + col];
//		}
//
//		float sum = 0.;
//		for (size_t row = 1; row < height-1; row++) {
//			for (size_t col = 1; col < width-1; col++) {
//				sum += input[row*width + col];
//				//printf("%.1f ", input[row*width + col]);
//			}
//			//printf("\n");
//		}
//		float avg_temp = sum/((width-2)*(height-2));
//		printf("average absolute difference: %E\n", avg_temp);
//
//		sum = 0.;
//		for (size_t row = 1; row < height-1; row++)
//			for (size_t col = 1; col < width-1; col++)
//				sum += fabs(input[row*width + col] - avg_temp);
//		avg_temp = sum/((width-2)*(height-2));
//		printf("average absolute difference: %E\n", avg_temp);
//	}


//	free(send_buffer);
//	free(recv_buffer);
//	free(working_buffer);
//	free(input);
	free_shifted_matrix(buff0);
	free_shifted_matrix(buff1);
//	MPI_Finalize();

	return 0;
}

