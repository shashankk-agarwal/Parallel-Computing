#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){
	int  myrank, size, n ;
	double stime, etime, time ;
	double maxtime;
	n = sqrt(atoi(argv[1])) ;

	MPI_Init(&argc, &argv) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
	MPI_Comm_size(MPI_COMM_WORLD, &size) ;

	int timesteps = atoi(argv[2]) ;
	double (*mat)[n][n] =malloc((timesteps+1) * sizeof *mat);
	double (*new_mat)[4][n]=malloc((timesteps) * sizeof *new_mat);

	for(int j=0 ; j<n ; j++){
		for(int k=0 ; k<n ; k++){
			mat[0][j][k] = 4.00 ;
		}
	}

	long int  size_n = sqrt(size) ;

	long int  row_p = myrank/size_n ;
	long int  col_p = myrank%size_n ;

	int ts = 0 ;
	
	MPI_Request request[timesteps][4] ;
	MPI_Request request2[timesteps][4] ;	
	MPI_Status status[timesteps][4] ;

	double (*buffer)[4][n] = malloc((timesteps) * sizeof *buffer);
	double (*output)[4][n] = malloc((timesteps) * sizeof *output);

	int position ;
	stime = MPI_Wtime() ;

	while(ts < timesteps){
		if(col_p==0){
		
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Pack(&mat[ts][i][n-1], 1, MPI_DOUBLE, buffer[ts][2], 8*n, &position, MPI_COMM_WORLD);
			}
			MPI_Isend(buffer[ts][2], position, MPI_PACKED, myrank+1, 99, MPI_COMM_WORLD, &request[ts][2]) ;
			MPI_Irecv(output[ts][2], 8*n, MPI_PACKED, myrank+1, 99, MPI_COMM_WORLD, &request2[ts][2]) ;
			
			if(row_p>0){
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Pack(&mat[ts][0][i], 1, MPI_DOUBLE, buffer[ts][1], 8*n, &position, MPI_COMM_WORLD);
				}
				MPI_Isend(buffer[ts][1], position, MPI_PACKED, myrank-size_n, 99, MPI_COMM_WORLD, &request[ts][1]) ;
				MPI_Irecv(output[ts][1], 8*n, MPI_PACKED, myrank-size_n, 99, MPI_COMM_WORLD, &request2[ts][1]) ;
			}
			
			if(row_p<size_n-1){
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Pack(&mat[ts][n-1][i], 1, MPI_DOUBLE, buffer[ts][3], 8*n, &position, MPI_COMM_WORLD);
				}
				MPI_Isend(buffer[ts][3], position, MPI_PACKED, myrank+size_n, 99, MPI_COMM_WORLD, &request[ts][3]) ;
				MPI_Irecv(output[ts][3], 8*n, MPI_PACKED, myrank+size_n, 99, MPI_COMM_WORLD, &request2[ts][3]) ;
			}
		}

		else if(col_p==size_n-1){
			
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Pack(&mat[ts][i][0], 1, MPI_DOUBLE, buffer[ts][0], 8*n, &position, MPI_COMM_WORLD);
			}
			MPI_Isend(buffer[ts][0], position, MPI_PACKED, myrank-1, 99, MPI_COMM_WORLD, &request[ts][0]) ;                      
			MPI_Irecv(output[ts][0], 8*n, MPI_PACKED, myrank-1, 99, MPI_COMM_WORLD, &request2[ts][0]) ;
			
                        if(row_p>0){
                        	position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Pack(&mat[ts][0][i], 1, MPI_DOUBLE, buffer[ts][1], 8*n, &position, MPI_COMM_WORLD);
				}
				MPI_Isend(buffer[ts][1], position, MPI_PACKED, myrank-size_n, 99, MPI_COMM_WORLD, &request[ts][1]) ;
                               MPI_Irecv(output[ts][1], 8*n, MPI_PACKED, myrank-size_n, 99, MPI_COMM_WORLD, &request2[ts][1]) ;
                        }
                        
                        if(row_p<size_n-1){
                        	position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Pack(&mat[ts][n-1][i], 1, MPI_DOUBLE, buffer[ts][3], 8*n, &position, MPI_COMM_WORLD);
				}
				MPI_Isend(buffer[ts][3], position, MPI_PACKED, myrank+size_n, 99, MPI_COMM_WORLD, &request[ts][3]) ;
                               MPI_Irecv(output[ts][3], 8*n, MPI_PACKED, myrank+size_n, 99, MPI_COMM_WORLD, &request2[ts][3]) ;
                        }
		}
	
		else{
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Pack(&mat[ts][i][n-1], 1, MPI_DOUBLE, buffer[ts][2], 8*n, &position, MPI_COMM_WORLD);
			}
			MPI_Isend(buffer[ts][2], position, MPI_PACKED, myrank+1, 99, MPI_COMM_WORLD, &request[ts][2]) ;
			MPI_Irecv(output[ts][2], 8*n, MPI_PACKED, myrank+1, 99, MPI_COMM_WORLD, &request2[ts][2]) ;
		
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Pack(&mat[ts][i][0], 1, MPI_DOUBLE, buffer[ts][0], 8*n, &position, MPI_COMM_WORLD);
			}
			MPI_Isend(buffer[ts][0], position, MPI_PACKED, myrank-1, 99, MPI_COMM_WORLD, &request[ts][0]) ;                      
			MPI_Irecv(output[ts][0], 8*n, MPI_PACKED, myrank-1, 99, MPI_COMM_WORLD, &request2[ts][0]) ;			
                        
                        if(row_p>0){
                        	position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Pack(&mat[ts][0][i], 1, MPI_DOUBLE, buffer[ts][1], 8*n, &position, MPI_COMM_WORLD);
				}
				MPI_Isend(buffer[ts][1], position, MPI_PACKED, myrank-size_n, 99, MPI_COMM_WORLD, &request[ts][1]) ;
                               MPI_Irecv(output[ts][1], 8*n, MPI_PACKED, myrank-size_n, 99, MPI_COMM_WORLD, &request2[ts][1]) ;
                        }
                        
                        if(row_p<size_n-1){
                        	position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Pack(&mat[ts][n-1][i], 1, MPI_DOUBLE, buffer[ts][3], 8*n, &position, MPI_COMM_WORLD);
				}
				MPI_Isend(buffer[ts][3], position, MPI_PACKED, myrank+size_n, 99, MPI_COMM_WORLD, &request[ts][3]) ;
                               MPI_Irecv(output[ts][3], 8*n, MPI_PACKED, myrank+size_n, 99, MPI_COMM_WORLD, &request2[ts][3]) ;
                        }
		}
		
		if(col_p==0){
			
			MPI_Wait(&request2[ts][2], &status[ts][2]) ;
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Unpack(output[ts][2], 8*n, &position, &new_mat[ts][2][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
			}
			
			if(row_p > 0){	
				MPI_Wait(&request2[ts][1], &status[ts][1]) ;
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Unpack(output[ts][1], 8*n, &position, &new_mat[ts][1][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
				}
			}
			
			if(row_p < size_n - 1){
				MPI_Wait(&request2[ts][3], &status[ts][3]) ;
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Unpack(output[ts][3], 8*n, &position, &new_mat[ts][3][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
				}
			}
		}

		else if(col_p==size_n-1){
			
			MPI_Wait(&request2[ts][0], &status[ts][0]) ;
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Unpack(output[ts][0], 8*n, &position, &new_mat[ts][0][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
			}
			
                        if(row_p>0){
                               MPI_Wait(&request2[ts][1], &status[ts][1]) ;
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Unpack(output[ts][1], 8*n, &position, &new_mat[ts][1][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
				}
                        }
                        
                        if(row_p < size_n - 1){
				MPI_Wait(&request2[ts][3], &status[ts][3]) ;
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Unpack(output[ts][3], 8*n, &position, &new_mat[ts][3][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
				}
			}
		}
	
		else{
			MPI_Wait(&request2[ts][2], &status[ts][2]) ;
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Unpack(output[ts][2], 8*n, &position, &new_mat[ts][2][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
			}
			
			MPI_Wait(&request2[ts][0], &status[ts][0]) ;
			position = 0 ;
			for(int i=0 ; i<n ; i++){
				MPI_Unpack(output[ts][0], 8*n, &position, &new_mat[ts][0][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
			}
                        
                        if(row_p>0){
                               MPI_Wait(&request2[ts][1], &status[ts][1]) ;
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Unpack(output[ts][1], 8*n, &position, &new_mat[ts][1][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
				}
                        }
                        
                        if(row_p < size_n - 1){
				MPI_Wait(&request2[ts][3], &status[ts][3]) ;
				position = 0 ;
				for(int i=0 ; i<n ; i++){
					MPI_Unpack(output[ts][3], 8*n, &position, &new_mat[ts][3][i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
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
	free(buffer);
	free(output);
	
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


