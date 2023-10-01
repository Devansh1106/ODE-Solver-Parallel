#include <stdio.h>
#include "Master.hpp"
#include <stdlib.h>
#include "mpi.h"
#include "dmumps_c.h"

void Solver(int argc, char** argv)
{
	MUMPS_INT n;
	int lines_irn =0;
	int lines_jcn =0;
	int lines_rhs =0;
	MUMPS_INT8 lines_a =0;
	int* irn;
	int* jcn;
	double* rhs;
	double* a;
	double time_taken,time_taken_2; 
	DMUMPS_STRUC_C id;
	int myid, ierr;
	int error =0;
  	FILE *file_1 = fopen("irn.txt", "r");
  	FILE *file_2 = fopen("jcn.txt", "r");
  	FILE *file_3 = fopen("rhs.txt", "r");
  	FILE *file_4 = fopen("a.txt", "r");
  	if(file_1!=NULL && file_2!=NULL && file_3!=NULL && file_4!=NULL)
  	{
  		char ch;
	  	while((ch = fgetc(file_1)) != EOF)
	  	{
	  		if(ch == '\n')
	  		{
	  			lines_irn++;
	  		}
	  	}
	  	ch =fgetc(file_2);
	  	while((ch = fgetc(file_2)) != EOF)
	  	{
	  		if(ch == '\n')
	  		{
	  			lines_jcn++;
	  		}
	  	}
	  	ch =fgetc(file_3);
	  	while((ch = fgetc(file_3)) != EOF)
	  	{
	  		if(ch == '\n')
	  		{
	  			lines_rhs++;
	  		}
	  	}
	  	n = lines_rhs;
	  	ch =fgetc(file_4);
	  	while((ch = fgetc(file_4)) != EOF)
	  	{
	  		if(ch == '\n')
	  		{
	  			lines_a++;
	  		}
	  	}
	  	rewind(file_1);
	  	rewind(file_2);
	  	rewind(file_3);
	  	rewind(file_4);
		irn = (int*)malloc(lines_irn * sizeof(int));
		jcn = (int*)malloc(lines_jcn * sizeof(int));
		rhs = (double*)malloc(lines_rhs * sizeof(double));
		a = (double*)malloc(lines_a * sizeof(double));

	  	ierr = MPI_Init(&argc, &argv);		//MPI initialization
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		
		if (irn == NULL && jcn == NULL && rhs == NULL && a == NULL)
		{
			printf("Memory can't be allocated.");
		}
		else 
		{
		  	int i=0;
		  	double num;
		  	while(fscanf(file_1, "%lf", &num)!=EOF)
		  	{
		  		irn[i] = num;
		  		i++;
		  	}
		  	i=0;
		  	while(fscanf(file_2, "%lf", &num)!=EOF)
		  	{
		  		jcn[i] = num;
		  		i++;
		  	}
		  	i=0;
		  	while(fscanf(file_3, "%lf", &num)!=EOF)
		  	{
		  		rhs[i] = num;
		  		i++;
		  	}
		  	i=0;
		  	while(fscanf(file_4, "%lf", &num)!=EOF)
		  	{
		  		a[i] = num;
		  		i++;
		  	}
		  	// for(int i =0; i < 50; i++)
		  	// {
		  	// 	printf("%d\n", rhs[i]);
		  	// }
		}
  	}
  	else
  	{
  		printf("Can't open a file!");
  	}

  	fclose(file_1);
  	fclose(file_2);
  	fclose(file_3);
  	fclose(file_4);

	time_taken = MPI_Wtime();		//record the start time

	id.par =1; id.sym=0;
	id.job = -1;
	dmumps_c(&id);

	if(myid == 0)
	{
		id.n = n; id.nnz = lines_a; id.irn = irn; id.jcn = jcn;
		id.a = a; id.rhs = rhs;
	}
	#define ICNTL(I) icntl[(I)-1] 
	id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0; //Supressing error msgs

	//Call the MUMPS package (analyze , factorization and solve)
	id.job = 6;
	dmumps_c(&id);

	if (id.infog[0] < 0)
	{
		printf("(PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
		myid, id.infog[0], id.infog[1]);
		error = 1;
	}
	//Terminate instance
	id.job = -2;
	dmumps_c(&id);
	if (myid == 0)
	{
		if (!error)
	    {
	    	FILE *file_sol = fopen("model_prob_sol2.txt", "w");
	    	if(file_sol != NULL)
	    	{
	    		for(int i=0; i < lines_rhs; i++)
	    		{
	    			float d= (float)(i)/(float)(lines_rhs);
	    			fprintf(file_sol, "%f %lf\n",d, rhs[(int)(i)]);
	    		}
	    		fclose(file_sol);
	    		// printf("%d\n",lines_rhs);
	    		printf("Solution is in model_prob_sol2.txt!");
	    	}
	    	else
	    	{
	    		printf("Error opening solution file!");
	    	}
	    }
	    else
	    {
	    	printf("An error occured, please check error code returnd by MUMPS.\n");
	    }
	}
	ierr = MPI_Finalize();

	time_taken_2 = MPI_Wtime(); //end time
	printf("\nTime Taken = %f",time_taken_2 - time_taken);
}