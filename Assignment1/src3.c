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
	double (*mat)[n][n]=malloc((timesteps+1) * sizeof *mat) ;
	double (*new_mat)[4][n]=malloc((timesteps) * sizeof *new_mat);

	for(int j=0 ; j<n ; j++){
		for(int k=0 ; k<n ; k++){
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

	long int size_n = sqrt(size) ;

	long int row_p = myrank/size_n ;
	long int col_p = myrank%size_n ;

	int ts = 0 ;
	
	MPI_Request request[timesteps][4] ;
	MPI_Request request2[timesteps][4] ;	
	MPI_Status status[timesteps][4] ;

	MPI_Datatype ROW ;
	MPI_Datatype COLUMN ;

	MPI_Type_contiguous(n, MPI_DOUBLE, &ROW) ;
	MPI_Type_commit(&ROW) ;

	MPI_Type_vector(n, 1, n, MPI_DOUBLE, &COLUMN) ;
	MPI_Type_commit(&COLUMN) ;

	double (*buf)[n][4] =malloc((timesteps) * sizeof *buf);

	stime = MPI_Wtime() ;
	while(ts < timesteps){
		if(col_p==0){
			MPI_Isend(&mat[ts][0][n-1], 1, COLUMN, myrank+1, 99, MPI_COMM_WORLD, &request[ts][2]) ;
		        MPI_Irecv(&new_mat[ts][2][0], n, MPI_DOUBLE, myrank+1, 99, MPI_COMM_WORLD, &request2[ts][2]) ;	
			
			if(row_p>0){
                         	MPI_Isend(&mat[ts][0][0], 1, ROW, myrank-size_n, 99, MPI_COMM_WORLD, &request[ts][1]) ;
				MPI_Irecv(&new_mat[ts][1][0], n, MPI_DOUBLE, myrank-size_n, 99, MPI_COMM_WORLD, &request2[ts][1]) ;
			}
			
			if(row_p<size_n-1){
                               MPI_Isend(&mat[ts][n-1][0], 1, ROW, myrank+size_n, 99, MPI_COMM_WORLD, &request[ts][3]) ;
				MPI_Irecv(&new_mat[ts][3][0], n, MPI_DOUBLE, myrank+size_n, 99, MPI_COMM_WORLD, &request2[ts][3]) ;
			}
		}

		else if(col_p==size_n-1){
			MPI_Isend(&mat[ts][0][0], 1, COLUMN, myrank-1, 99, MPI_COMM_WORLD, &request[ts][0]) ;
			MPI_Irecv(&new_mat[ts][0][0], n, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &request2[ts][0]) ;
			
                        if(row_p>0){
                               MPI_Isend(&mat[ts][0][0], 1, ROW, myrank-size_n, 99, MPI_COMM_WORLD, &request[ts][1]) ;
				MPI_Irecv(&new_mat[ts][1][0], n, MPI_DOUBLE, myrank-size_n, 99, MPI_COMM_WORLD, &request2[ts][1]) ;
                        }
                        
                        if(row_p<size_n-1){
                               MPI_Isend(&mat[ts][n-1][0], 1, ROW, myrank+size_n, 99, MPI_COMM_WORLD, &request[ts][3]) ;
				MPI_Irecv(&new_mat[ts][3][0], n, MPI_DOUBLE, myrank+size_n, 99, MPI_COMM_WORLD, &request2[ts][3]) ;
                        }
		}
	
		else{
			MPI_Isend(&mat[ts][0][n-1], 1, COLUMN, myrank+1, 99, MPI_COMM_WORLD, &request[ts][2]) ;
			MPI_Irecv(&new_mat[ts][2][0], n, MPI_DOUBLE, myrank+1, 99, MPI_COMM_WORLD, &request2[ts][2]) ;

			MPI_Isend(&mat[ts][0][0], 1, COLUMN, myrank-1, 99, MPI_COMM_WORLD, &request[ts][0]) ;
			MPI_Irecv(&new_mat[ts][0][0], n, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &request2[ts][0]) ;

                        if(row_p>0){
                               MPI_Isend(&mat[ts][0][0], 1, ROW, myrank-size_n, 99, MPI_COMM_WORLD, &request[ts][1]) ;
				MPI_Irecv(&new_mat[ts][1][0], n, MPI_DOUBLE, myrank-size_n, 99, MPI_COMM_WORLD, &request2[ts][1]) ;
                        }
                        
                        if(row_p<size_n-1){
                               MPI_Isend(&mat[ts][n-1][0], 1, ROW, myrank+size_n, 99, MPI_COMM_WORLD, &request[ts][3]) ;
				MPI_Irecv(&new_mat[ts][3][0], n, MPI_DOUBLE, myrank+size_n, 99, MPI_COMM_WORLD, &request2[ts][3]) ;
                        }
		}
		
		if(col_p==0){
			MPI_Wait(&request2[ts][2], &status[ts][2]) ;
			
			if(row_p>0){
				MPI_Wait(&request2[ts][1], &status[ts][1]) ;
			}
			if(row_p<size_n-1){
				MPI_Wait(&request2[ts][3], &status[ts][3]) ;
			}
		}

		else if(col_p==size_n-1){
			MPI_Wait(&request2[ts][0], &status[ts][0]) ;

                       if(row_p>0){
                               MPI_Wait(&request2[ts][1], &status[ts][1]) ;
                       }
                       
                       if(row_p<size_n-1){
                               MPI_Wait(&request2[ts][3], &status[ts][3]) ;
                       }
		}
	
		else{
			MPI_Wait(&request2[ts][2], &status[ts][2]) ;
	
			MPI_Wait(&request2[ts][0], &status[ts][0]) ;
			
                        if(row_p>0){
                               MPI_Wait(&request2[ts][1], &status[ts][1]) ;
                        }
                        if(row_p<size_n-1){
                               MPI_Wait(&request2[ts][3], &status[ts][3]) ;
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
	free(buf);


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

