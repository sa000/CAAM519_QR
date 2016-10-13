

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>


#include <stdio.h>
static char help[]= "QR stuf \n\n";
int main(int argc, char **args) {

 
    Mat A,Q,B,E, QT,R, R_calc, A_calc;
    Vec a,x,y, u,e, a_const, e1, e2;
    PetscScalar dot;
    PetscScalar neg=-1;
    PetscInt n=3, b=0, i=1, i0=0, i1=1, i2=2;

    
    PetscReal *norms;
    PetscMalloc1(n,&norms);
    PetscScalar uvals;

    const PetscScalar *vals;
    const PetscInt    *cols;
	const PetscScalar *evals;
    PetscInt ncols, rstart, rend;
    PetscInt col[]={0,1,2};
    PetscScalar val0[]= {1,1,0};
    PetscScalar val1[]= {1,0,1};
    PetscScalar val2[]= {0,1,1};
    PetscInitialize(&argc, &args, (char*)0, help);

    MatCreate(MPI_COMM_SELF, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(A);
    MatSetUp(A);

    MatCreate(MPI_COMM_SELF, &Q);
    MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(Q);
    MatSetUp(Q);

    MatCreate(MPI_COMM_SELF, &QT);
    MatSetSizes(QT, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(QT);
    MatSetUp(QT);

    MatCreate(MPI_COMM_SELF, &B);
    MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(B);
    MatSetUp(B);


    MatCreate(MPI_COMM_SELF, &R);
    MatSetSizes(R, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(R);
    MatSetUp(R);

    MatCreate(MPI_COMM_SELF, &R_calc);
    MatSetSizes(R_calc, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(R_calc);
    MatSetUp(R_calc);

    	//Create sample vector
	 VecCreate(PETSC_COMM_WORLD,&x);
	 VecSetSizes(x,n,n);
	 VecSetFromOptions(x);

	 VecCreate(PETSC_COMM_WORLD,&a_const);
	 VecSetSizes(a_const,n,n);
	 VecSetFromOptions(a_const);

	 VecCreate(PETSC_COMM_WORLD,&u);
	 VecSetSizes(u,n,n);
	 VecSetFromOptions(u);

	 VecCreate(PETSC_COMM_WORLD,&e);
	 VecSetSizes(e,n,n);
	 VecSetFromOptions(e);


	 VecCreate(PETSC_COMM_WORLD,&e1);
	 VecSetSizes(e1,n,n);
	 VecSetFromOptions(e1);

	 VecCreate(PETSC_COMM_WORLD,&e2);
	 VecSetSizes(e2,n,n);
	 VecSetFromOptions(e2);

	 VecCreate(PETSC_COMM_WORLD,&a);
	 VecSetSizes(a,n,n);
	 VecSetFromOptions(a);

	 MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);

    MatCreate(MPI_COMM_SELF, &E);
    MatSetSizes(E, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(E);
    MatSetUp(E);

    MatCreate(MPI_COMM_SELF, &A_calc);
    MatSetSizes(A_calc, PETSC_DECIDE, PETSC_DECIDE, n,n);
    MatSetFromOptions(A_calc);
    MatSetUp(A_calc);
    //Build matrix, #of rows, their indicies, # of columns, their indices, and values
   	MatSetValues(A,1, &i0, n, &col, val0, INSERT_VALUES);
   	MatSetValues(A,1, &i1, n, &col, val1, INSERT_VALUES);
   	MatSetValues(A,1, &i2, n, &col, val2, INSERT_VALUES);

 	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
 	printf("Matrix A:\n");
 	MatView(A, PETSC_VIEWER_DEFAULT);

 	//Norm of each column
 	MatGetColumnNorms(A,NORM_2,norms);

 	MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
 	MatTranspose(A, MAT_INITIAL_MATRIX , &B);
 	MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
 	//printf("here is B, the tranpose");
 	//MatView(B, PETSC_VIEWER_DEFAULT);

	

	

	//grab the row of the matrix
	PetscReal      nrm;
	PetscReal      tmpnrm;
	PetscScalar yz;
	//Algorithm setup
	for(int i=0; i<n; i++){
		MatGetRow(B, i, NULL, NULL, &vals);
		VecSetValues(a,n,col,&vals[0], INSERT_VALUES);
		//Calculate u
		//Go thru subloop 

		if (i>0){
			VecCopy(a,a_const);
			for (int j=0; j<i; j++){
				//Grab the previous col
				if (i==2){
					if(j==0){
						VecGetValues(e1,3,col,evals);
					}
					else{
						VecGetValues(e2,3,col,evals);

					}
				}
				else{ 
				MatGetRow(QT, j, NULL, NULL, &evals);
				}
				VecSetValues(e,n,col,&evals[0], INSERT_VALUES);
				//PetscScalarView(n,&evals[0],PETSC_VIEWER_STDOUT_WORLD);
			//	VecView(e,PETSC_VIEWER_DEFAULT 	);
				//We need to grab the original a
				VecDot(a_const,e,&dot);
				VecScale(e,neg);
				VecAXPY(a,dot,e);

			}
		}
		else{	
			//Use this as u for the first case	
		}
		VecCopy(a,u);
		VecNorm(u,NORM_2,&nrm);
		tmpnrm=1.0/nrm;
		VecScale(u,tmpnrm);//this is now e, this is u/norm(u)
		//VecView(u,PETSC_VIEWER_DEFAULT 	);
		VecGetValues(u,3,col, &uvals);//transfers values to uvals to insert into q
		MatSetValues(Q,1, &i, n, &col, &uvals, INSERT_VALUES);
		if(i==0){
			VecCopy(u,e1);
			//VecView(e1,PETSC_VIEWER_DEFAULT 	);
		}
		if(i==1){
			VecCopy(u,e2);
			//VecView(e1,PETSC_VIEWER_DEFAULT 	);
			//VecView(e2,PETSC_VIEWER_DEFAULT 	);
		}
		MatSetValues(QT,n, &col, 1, &i, &uvals, INSERT_VALUES);

		//MatTranspose(Q, MAT_INITIAL_MATRIX , &QT);
			MatAssemblyBegin(QT, MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd(QT, MAT_FINAL_ASSEMBLY);	
		MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);


	}
	//Assemble R
	MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY);
	MatZeroEntries(R);
	for (int i=0; i<n; i++){
		MatGetRow(Q, i, NULL, NULL, &vals); //get the ith column of the Original Q matrix. 
		VecSetValues(e,n,col,&vals[0], INSERT_VALUES);
		//VecView(e,PETSC_VIEWER_DEFAULT 	);
		for (int j=0; j<n; j++){
			MatGetRow(B, j, NULL, NULL, &vals); //get the jth row. 
			VecSetValues(a,n,col,&vals[0], INSERT_VALUES);
			//VecView(a,PETSC_VIEWER_DEFAULT 	);
			if (j>=i){

				VecDot(a,e,&dot);
				MatSetValue(R, i, j, dot, INSERT_VALUES);
			}
		}
	}
 	printf("\n");
	MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY);
	printf("QR Factorization yields the following \n");
	printf("Here is R \n\n");
	MatView(R, PETSC_VIEWER_DEFAULT);
	printf("\n");	
 		//Switched Q & QT
 		printf("Here is Q \n");
		MatView(QT, PETSC_VIEWER_DEFAULT);
		printf("\n");
 	//PetscScalarView(3,&vals[5], PETSC_VIEWER_DEFAULT);

 	//MatView(Q, PETSC_VIEWER_DEFAULT);

	MatAssemblyBegin(R_calc, MAT_FINAL_ASSEMBLY);
	MatMatMult(A, QT, MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R_calc);
	MatTranspose(R_calc,MAT_REUSE_MATRIX, &R_calc);
	MatAssemblyEnd(R_calc, MAT_FINAL_ASSEMBLY);	
	printf("We need to verify our factorization is correct \n");
 	printf("Here is R Calculated via Q^T*A \n");
  	MatView(R_calc, PETSC_VIEWER_DEFAULT);
  	printf("\n");

 	MatAssemblyBegin(A_calc, MAT_FINAL_ASSEMBLY);
	MatMatMult(QT, R, MAT_INITIAL_MATRIX,PETSC_DEFAULT,&A_calc);
    MatAssemblyEnd(R_calc, MAT_FINAL_ASSEMBLY);
	
 	printf("Here is A Calculated via Q*R \n");
  	MatView(A_calc, PETSC_VIEWER_DEFAULT);
  	printf("\n");
  	//More verification
  	MatAXPY(A_calc, neg, A, SAME_NONZERO_PATTERN);
  	printf("\n");
  	printf("here is the normed difference ");
  	//MatView(A, PETSC_VIEWER_DEFAULT);
  	MatView(A_calc, PETSC_VIEWER_DEFAULT);



     PetscFinalize();
	}
