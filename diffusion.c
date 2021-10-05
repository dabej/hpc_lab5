#define _XOPEN_SOURCE 700

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char * argv[]) {
	MPI_Init(&argc, &argv);

	int nmb_mpi_proc, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nmb_mpi_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	int width, height, n;
	float d;
	float *input;

	// Parse arguments and init-file
	if (mpi_rank == 0) {

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

		FILE *fp = fopen("init_100_100", "r");
		fscanf(fp, "%d %d", &width, &height);
		width += 2;
		height += 2;

		input =  malloc(sizeof(float)*width*height);
		for(int i = 0; i < width*height; i++)
			input[i] = 0.;

		int row, col;
		float temp;
		while (fscanf(fp, "%d %d %f", &col, &row, &temp) == 3) {
			input[(row+1) * width + (col+1)] = temp;
		}

		fclose(fp);
	}

	MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate input for all other nodes
	if (mpi_rank != 0)
		input = malloc(width*height*sizeof(float));

	MPI_Bcast(input, width*height, MPI_FLOAT, 0, MPI_COMM_WORLD);

	float *send_buffer = malloc(width*sizeof(float));
	float *recv_buffer = malloc(width*sizeof(float));
	int rows = (height-3) / nmb_mpi_proc + 1;
	float *working_buffer = malloc(width*height*sizeof(float));

	int above = mpi_rank - 1;
	int below = mpi_rank + 1;
	MPI_Status status;
	float up, down, left, right, value;
	int start_row = mpi_rank*rows + 1;
	int end_row = (mpi_rank+1) * rows + 1;
	if (end_row >= height)
		end_row = height - 1;
	//float *up_row, *down_row;

	for (size_t iter = 0; iter < n; iter++) {
		for (size_t row = start_row; row < end_row; row++) {
			//up_row = &input[(row-1)*width];
			//down_row = &input[(row+1)*width];
			for (size_t col = 1; col < width-1; col++) {
				value = input[row*width + col];
				up = input[(row-1)*width + col];
				down = input[(row+1)*width + col];
				left = input[row*width + col-1];
				right = input[row*width + col+1];
				value += d * ((up + down + left + right)/4 - value);
				working_buffer[row*width + col] = value;
			}
		}

		for (size_t row = start_row; row < end_row; row++)
			for (size_t col = 1; col < width-1; col++)
				input[row*width + col] = working_buffer[row*width + col];

		if (start_row != 1) {
			for (size_t col = 1; col < width-1; col++)
				send_buffer[col] = working_buffer[start_row*width + col];
			MPI_Sendrecv(send_buffer, width, MPI_FLOAT, above, iter,
						 recv_buffer, width, MPI_FLOAT, above, iter,
						 MPI_COMM_WORLD, &status);	
			for (size_t col = 1; col < width-1; col++)
				input[(start_row-1)*width + col] = recv_buffer[col];
		}

		if (end_row != height-1) {
			for (size_t col = 1; col < width-1; col++)
				send_buffer[col] = working_buffer[(end_row-1)*width + col];
			MPI_Sendrecv(send_buffer, width, MPI_FLOAT, below, iter,
						 recv_buffer, width, MPI_FLOAT, below, iter,
						 MPI_COMM_WORLD, &status);	
			for (size_t col = 1; col < width-1; col++)
				input[end_row*width + col] = recv_buffer[col];
		}
	}

	if (mpi_rank != 0)
		MPI_Send(input, width*height, MPI_FLOAT, 0, mpi_rank, MPI_COMM_WORLD);
	else {
		for (size_t proc = 1; proc < nmb_mpi_proc; proc++) {
			MPI_Recv(working_buffer, width*height, MPI_FLOAT, proc, proc,
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int start_row = proc*rows + 1;
			int end_row = (proc+1) * rows + 1;
			if (end_row >= height)
				end_row = height - 1;
			for (size_t row = start_row; row < end_row; row++)
				for (size_t col = 1; col < width-1; col++)
					input[row*width + col] = working_buffer[row*width + col];
		}

		float sum = 0.;
		for (size_t row = 1; row < height-1; row++) {
			for (size_t col = 1; col < width-1; col++) {
				sum += input[row*width + col];
				//printf("%.1f ", input[row*width + col]);
			}
			//printf("\n");
		}
		float avg_temp = sum/((width-2)*(height-2));
		printf("average absolute difference: %E\n", avg_temp);

		sum = 0.;
		for (size_t row = 1; row < height-1; row++)
			for (size_t col = 1; col < width-1; col++)
				sum += fabs(input[row*width + col] - avg_temp);
		avg_temp = sum/((width-2)*(height-2));
		printf("average absolute difference: %E\n", avg_temp);
	}

	free(send_buffer);
	free(recv_buffer);
	free(working_buffer);
	free(input);
	MPI_Finalize();

	return 0;
}

