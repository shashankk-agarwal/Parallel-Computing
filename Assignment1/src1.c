#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){
	int myrank, size, n ;
	double stime, etime, time ;
	double maxtime;	
	n = sqrt(atoi(argv[1])) ;

	MPI_Init(&argc, &argv) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
	MPI_Comm_size(MPI_COMM_WORLD, &size) ;

	int timesteps = atoi(argv[2]) ;
	double (*mat)[n][n] = malloc((timesteps+1) * sizeof *mat);
	double (*new_mat)[4][n]=malloc((timesteps) * sizeof *new_mat);

	for(long int j=0 ; j<n ; j++){
		for(long int k=0 ; k<n ; k++){
			mat[0][j][k] = 4.00 ;
		}
	}
/*
	for(int j=0 ; j<4 ; j++){
                for(int k=0 ; k<n ; k++){
                        new_mat[0][j][k] = -1.00 ;
                }
        }
*/
	int size_n = sqrt(size) ;

	int row_p = myrank/size_n ;
	int col_p = myrank%size_n ;

	int ts = 0 ;
	
	MPI_Request (*request)[4][n] = malloc((timesteps) * sizeof *request);
	MPI_Request (*request2)[4][n] = malloc((timesteps) * sizeof *request2);	
	MPI_Status (*status)[4][n] = malloc((timesteps) * sizeof *status);

	/*
	int arr[2][2][n] ;
	int n_arr[2][2][n] ;
	MPI_Status stat[2][2][n] ;
	MPI_Request req[2][2][n] ;
	MPI_Request req1[2][2][n] ;

	if(myrank == 0){
		for(int i=0 ; i<n ; i++){
			arr[0][0][i] = 11 + i ;
			MPI_Isend(&arr[0][0][i], 1, MPI_INT, myrank+1, 99, MPI_COMM_WORLD, &req[0][0][i]) ;
		}
	}
	if(myrank == 1){
		for(int i=0 ; i<n ; i++){
			MPI_Irecv(&n_arr[0][0][i], 1, MPI_INT, myrank-1, 99, MPI_COMM_WORLD, &req1[0][0][i]) ;
			MPI_Wait(&req1[0][0][i], &stat[0][0][i]) ;
			printf("%d ", n_arr[0][0][i]) ;
		}
	}
	*/
	
	stime = MPI_Wtime() ;
	while(ts < timesteps){
		if(col_p==0){
			for(int i=0 ; i<n ; i++){
				MPI_Isend(&mat[ts][i][n-1], 1, MPI_DOUBLE, myrank+1, i, MPI_COMM_WORLD, &request[ts][2][i]) ;
				MPI_Irecv(&new_mat[ts][2][i], 1, MPI_DOUBLE, myrank+1, i, MPI_COMM_WORLD, &request2[ts][2][i]) ;
			}
			if(row_p>0){
				for(int i=0 ; i<n ; i++){
                         	      	MPI_Isend(&mat[ts][0][i], 1, MPI_DOUBLE, myrank-size_n, i, MPI_COMM_WORLD, &request[ts][1][i]) ;
					MPI_Irecv(&new_mat[ts][1][i], 1, MPI_DOUBLE, myrank-size_n, i, MPI_COMM_WORLD, &request2[ts][1][i]) ;
				}
			}
			if(row_p<size_n-1){
				for(int i=0 ; i<n ; i++){
                                        MPI_Isend(&mat[ts][n-1][i], 1, MPI_DOUBLE, myrank+size_n, i, MPI_COMM_WORLD, &request[ts][3][i]) ;
				        MPI_Irecv(&new_mat[ts][3][i], 1, MPI_DOUBLE, myrank+size_n, i, MPI_COMM_WORLD, &request2[ts][3][i]) ;
				}
			}
		}

		else if(col_p==size_n-1){
			for(int i=0 ; i<n ; i++){
                                MPI_Isend(&mat[ts][i][0], 1, MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD, &request[ts][0][i]) ;
				MPI_Irecv(&new_mat[ts][0][i], 1, MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD, &request2[ts][0][i]) ;
                        }
                        if(row_p>0){
                                for(int i=0 ; i<n ; i++){
                                        MPI_Isend(&mat[ts][0][i], 1, MPI_DOUBLE, myrank-size_n, i, MPI_COMM_WORLD, &request[ts][1][i]) ;
					MPI_Irecv(&new_mat[ts][1][i], 1, MPI_DOUBLE, myrank-size_n, i, MPI_COMM_WORLD, &request2[ts][1][i]) ;
                                }
                        }
                        if(row_p<size_n-1){
                                for(int i=0 ; i<n ; i++){
                                        MPI_Isend(&mat[ts][n-1][i], 1, MPI_DOUBLE, myrank+size_n, i, MPI_COMM_WORLD, &request[ts][3][i]) ;
					MPI_Irecv(&new_mat[ts][3][i], 1, MPI_DOUBLE, myrank+size_n, i, MPI_COMM_WORLD, &request2[ts][3][i]) ;
                                }
                        }
		}
	
		else{
			for(int i=0 ; i<n ; i++){
                                MPI_Isend(&mat[ts][i][n-1], 1, MPI_DOUBLE, myrank+1, i, MPI_COMM_WORLD, &request[ts][2][i]) ;
				MPI_Irecv(&new_mat[ts][2][i], 1, MPI_DOUBLE, myrank+1, i, MPI_COMM_WORLD, &request2[ts][2][i]) ;
                        }
			for(int i=0 ; i<n ; i++){
                                MPI_Isend(&mat[ts][i][0], 1, MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD, &request[ts][0][i]) ;
				MPI_Irecv(&new_mat[ts][0][i], 1, MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD, &request2[ts][0][i]) ;
                        }
                        if(row_p>0){
                                for(int i=0 ; i<n ; i++){
                                        MPI_Isend(&mat[ts][0][i], 1, MPI_DOUBLE, myrank-size_n, i, MPI_COMM_WORLD, &request[ts][1][i]) ;
					MPI_Irecv(&new_mat[ts][1][i], 1, MPI_DOUBLE, myrank-size_n, i, MPI_COMM_WORLD, &request2[ts][1][i]) ;
                                }
                        }
                        if(row_p<size_n-1){
                                for(int i=0 ; i<n ; i++){
                                        MPI_Isend(&mat[ts][n-1][i], 1, MPI_DOUBLE, myrank+size_n, i, MPI_COMM_WORLD, &request[ts][3][i]) ;
					MPI_Irecv(&new_mat[ts][3][i], 1, MPI_DOUBLE, myrank+size_n, i, MPI_COMM_WORLD, &request2[ts][3][i]) ;
                                }
                        }
		}
		
		if(col_p==0){
			for(int i=0 ; i<n ; i++){
				MPI_Wait(&request2[ts][2][i], &status[ts][2][i]) ;
			}
			if(row_p>0){
				for(int i=0 ; i<n ; i++){
					MPI_Wait(&request2[ts][1][i], &status[ts][1][i]) ;
				}
			}
			if(row_p<size_n-1){
				for(int i=0 ; i<n ; i++){
					 MPI_Wait(&request2[ts][3][i], &status[ts][3][i]) ;
				}
			}
		}

		else if(col_p==size_n-1){
			for(int i=0 ; i<n ; i++){
				MPI_Wait(&request2[ts][0][i], &status[ts][0][i]) ;
                        }
                        if(row_p>0){
                                for(int i=0 ; i<n ; i++){
					MPI_Wait(&request2[ts][1][i], &status[ts][1][i]) ;
                                }
                        }
                        if(row_p<size_n-1){
                                for(int i=0 ; i<n ; i++){
					MPI_Wait(&request2[ts][3][i], &status[ts][3][i]) ;
                                }
                        }
		}
	
		else{
			for(int i=0 ; i<n ; i++){
				MPI_Wait(&request2[ts][2][i], &status[ts][2][i]) ;
                        }
			for(int i=0 ; i<n ; i++){
				MPI_Wait(&request2[ts][0][i], &status[ts][0][i]) ;
                        }
                        if(row_p>0){
                                for(int i=0 ; i<n ; i++){
					MPI_Wait(&request2[ts][1][i], &status[ts][1][i]) ;
                                }
                        }
                        if(row_p<size_n-1){
                                for(int i=0 ; i<n ; i++){
					MPI_Wait(&request2[ts][3][i], &status[ts][3][i]) ;
                                }
                        }
		}
		ts++ ;

		for(int j=1 ; j<n-1 ; j++){
			if(row_p == 0){
				mat[ts][0][j] = (mat[ts-1][1][j] + mat[ts-1][0][j-1] + mat[ts-1][0][j+1])/4 ;
			}
			else{
				mat[ts][0][j] = (mat[ts-1][1][j] + mat[ts-1][0][j-1] + mat[ts-1][0][j+1] + new_mat[ts-1][1][j])/4 ;
			}

			if(row_p == size_n-1){
                                mat[ts][n-1][j] = (mat[ts-1][n-2][j] + mat[ts-1][n-1][j-1] + mat[ts-1][n-1][j+1])/4 ;
                        }
                        else{
                                mat[ts][n-1][j] = (mat[ts-1][n-2][j] + mat[ts-1][n-1][j-1] + mat[ts-1][n-1][j+1] + new_mat[ts-1][3][j])/4 ;
                        }
		}

		for(int i=1 ; i<n-1 ; i++){
                        if(col_p == 0){
                                mat[ts][i][0] = (mat[ts-1][i-1][0] + mat[ts-1][i+1][0] + mat[ts-1][i][1])/4 ;
                        }
                        else{
				mat[ts][i][0] = (mat[ts-1][i-1][0] + mat[ts-1][i+1][0] + mat[ts-1][i][1] + new_mat[ts-1][0][i])/4 ;
                        }

                        if(col_p == size_n-1){
                                mat[ts][i][n-1] = (mat[ts-1][i-1][n-1] + mat[ts-1][i+1][n-1] + mat[ts-1][i][n-2])/4 ;
                        }
                        else{
                                mat[ts][i][n-1] = (mat[ts-1][i-1][n-1] + mat[ts-1][i+1][n-1] + mat[ts-1][i][n-2] + new_mat[ts-1][2][i])/4 ;
                        }
                }

		for(int i=1 ; i<n-1 ; i++){
			for(int j=1 ; j<n-1 ; j++){
				mat[ts][i][j] = (mat[ts-1][i-1][j] + mat[ts-1][i+1][j] + mat[ts-1][i][j-1] + mat[ts-1][i][j+1])/4 ;
			}
		}

		if(col_p == 0){
			if(row_p == 0){
				mat[ts][0][0] = (mat[ts-1][0][1] + mat[ts-1][1][0])/4 ;
			}
			else{
				mat[ts][0][0] = (mat[ts-1][0][1] + mat[ts-1][1][0] + new_mat[ts-1][1][0])/4 ;
			}
		}
		else{
			if(row_p == 0){
                                mat[ts][0][0] = (mat[ts-1][0][1] + mat[ts-1][1][0] + new_mat[ts-1][0][0])/4 ;
                        }
                        else{
                                mat[ts][0][0] = (mat[ts-1][0][1] + mat[ts-1][1][0] + new_mat[ts-1][1][0] + new_mat[ts-1][0][0])/4 ;
                        }
		}

		if(col_p == 0){
                        if(row_p == size_n - 1){
                                mat[ts][n-1][0] = (mat[ts-1][n-1][1] + mat[ts-1][n-2][0])/4 ;
                        }
                        else{
                                mat[ts][n-1][0] = (mat[ts-1][n-1][1] + mat[ts-1][n-2][0] + new_mat[ts-1][3][0])/4 ;
                        }
                }
                else{
                        if(row_p == size_n-1){
                                mat[ts][n-1][0] = (mat[ts-1][n-1][1] + mat[ts-1][n-2][0] + new_mat[ts-1][0][n-1])/4 ;
                        }
                        else{
                                mat[ts][n-1][0] = (mat[ts-1][n-1][1] + mat[ts-1][n-2][0] + new_mat[ts-1][0][n-1] + new_mat[ts-1][3][0])/4 ;
                        }
                }

		if(col_p == size_n-1){
                        if(row_p == 0){
                                mat[ts][0][n-1] = (mat[ts-1][0][n-2] + mat[ts-1][1][n-1])/4 ;
                        }
                        else{
                                mat[ts][0][n-1] = (mat[ts-1][0][n-2] + mat[ts-1][1][n-1] + new_mat[ts-1][1][n-1])/4 ;
                        }
                }
                else{
                        if(row_p == 0){
                                mat[ts][0][n-1] = (mat[ts-1][0][n-2] + mat[ts-1][1][n-1] + new_mat[ts-1][2][0])/4 ;
                        }
                        else{
                                mat[ts][0][n-1] = (mat[ts-1][0][n-2] + mat[ts-1][1][n-1] + new_mat[ts-1][2][0] + new_mat[ts-1][1][n-1])/4 ;
                        }
                }

		if(col_p == size_n-1){
                        if(row_p == size_n-1){
                                mat[ts][n-1][n-1] = (mat[ts-1][n-1][n-2] + mat[ts-1][n-2][n-1])/4 ;
                        }
                        else{
                                mat[ts][n-1][n-1] = (mat[ts-1][n-1][n-2] + mat[ts-1][n-2][n-1] + new_mat[ts-1][3][n-1])/4 ;
                        }
                }
                else{
                        if(row_p == size_n-1){
                                mat[ts][n-1][n-1] = (mat[ts-1][n-1][n-2] + mat[ts-1][n-2][n-1] + new_mat[ts-1][2][n-1])/4 ;
                        }
                        else{
                                mat[ts][n-1][n-1] = (mat[ts-1][n-1][n-2] + mat[ts-1][n-2][n-1] + new_mat[ts-1][2][n-1] + new_mat[ts-1][3][n-1])/4 ;
                        }
                }

	}
	
	etime = MPI_Wtime() ;
	time = etime - stime ;
	
	
	
	 // obtain max time
          MPI_Reduce (&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
          if (!myrank) printf ("%lf\n", maxtime);
	
	free(mat);
	free(new_mat);
	free(request);
	free(request2);
	free(status);
	/*
	for(int i=0 ; i<n ; i++){
		for(int j=0 ; j<n ; j++){
			printf("%lf ", mat[timesteps][i][j]) ;
		}
		printf("\n") ;
	}
	*/
	MPI_Finalize() ;
}
