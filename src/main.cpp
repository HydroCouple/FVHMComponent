#include "stdafx.h"
#include "sparsemat.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <HYPRE.h>
#include <HYPRE_utilities.h>
#include <HYPRE_parcsr_ls.h>
#include <vector>
#include <HYPRE_krylov.h>
#include <HYPRE_seq_mv.h>
#include <HYPRE_parcsr_mv.h>
#include "_hypre_parcsr_mv.h"


double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}


void solutionExample()
{
  //  HYPRE_IJMatrix A;

  //  HYPRE_IJVector b;
  //  HYPRE_IJVector x;

  //  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, 1000 - 1, 0, 1000 -1, &A);
  //  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  //  HYPRE_IJMatrixInitialize(A);
  //  HYPRE_IJMatrixSetOMPFlag(A,1);

  //  HYPRE_IJVectorCreate(MPI_COMM_WORLD,0,1000 - 1,&b);
  //  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  //  HYPRE_IJVectorInitialize(b);

  //  HYPRE_IJVectorCreate(MPI_COMM_WORLD,0,1000 - 1,&x);
  //  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  //  HYPRE_IJVectorInitialize(x);

  //  //The arrays ncols and rows are of dimension nrows
  //  //jjand contain the number of columns in each row and the row indices, respective
  ////  int ncols[1000];
  ////  int rows[1000];
  ////  double values[1000];
  //  double bValues[1000];
  //  double xValues[1000];

  ////  std::vector<int> columnIndexes(1000000)

  //  SparseMatrix matrix(0,999,1000, 1000);

  //  for(int r = 500; r > -1 ; r--)
  //  {
  //    matrix.appendValue(r,r,1);//fRand(0.9,1.1));
  //    bValues[r] = fRand(-1000,1000);
  //    xValues[r] = 0;
  //  }

  //  for(int r = 501; r < 1000 ; r++)
  //  {
  //    matrix.appendValue(r-1,r,1.0);
  //    matrix.appendValue(r,r,1.0);
  //    bValues[r] = fRand(0,1000);
  //    xValues[r] = 0;
  //  }

  //  matrix.print();

  //  double* values = nullptr;
  //  int* colIndexes = nullptr;

  //  matrix.getDataByRow(values);
  //  matrix.getColumnIndexes(colIndexes);

  //  HYPRE_IJMatrixSetValues(A,matrix.numRows(), matrix.colsPerRow(), matrix.rows() ,colIndexes,values);
  //  HYPRE_IJMatrixAssemble(A);

  //  delete[] values;
  //  delete[] colIndexes;

  //  HYPRE_IJVectorSetValues(b,matrix.numRows(),matrix.rows(),bValues);
  //  HYPRE_IJVectorAssemble(b);

  //  HYPRE_IJVectorSetValues(x,matrix.numRows(),matrix.rows(),xValues);
  //  HYPRE_IJVectorAssemble(x);

  //  HYPRE_IJMatrixPrint(A,"AMatrix.txt");
  //  HYPRE_IJVectorPrint(b,"BVector.txt");


  //  HYPRE_ParCSRMatrix parcsr_A;
  //  HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

  //  HYPRE_ParVector par_b;
  //  HYPRE_IJVectorGetObject(b,(void**) &par_b);

  //  HYPRE_ParVector par_x;
  //  HYPRE_IJVectorGetObject(x,(void**) &par_x);


  //  /* Choose a solver and solve the system */

  //  int solver_id = 1;
  //  HYPRE_Solver solver, precond;

  //  /* AMG */
  //  if (solver_id == 0)
  //  {
  //    int num_iterations;
  //    double final_res_norm;

  //    /* Create solver */
  //    HYPRE_BoomerAMGCreate(&solver);

  //    /* Set some parameters (See Reference Manual for more parameters) */
  //    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
  //    HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
  //    HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
  //    HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
  //    HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
  //    HYPRE_BoomerAMGSetTol(solver, 1e-7);      /* conv. tolerance */
  //    /* Now setup and solve! */
  //    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
  //    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

  //    /* Run info - needed logging turned on */
  //    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
  //    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  //    printf("\n");
  //    printf("Iterations = %d\n", num_iterations);
  //    printf("Final Relative Residual Norm = %e\n", final_res_norm);
  //    printf("\n");


  //    /* Destroy solver */
  //    HYPRE_BoomerAMGDestroy(solver);
  //  }
  //  /* PCG */
  //  else if (solver_id == 50)
  //  {
  //    int num_iterations;
  //    double final_res_norm;

  //    /* Create solver */
  //    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

  //    /* Set some parameters (See Reference Manual for more parameters) */
  //    HYPRE_PCGSetMaxIter(solver, 10000); /* max iterations */
  //    HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
  //    HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  //    HYPRE_PCGSetPrintLevel(solver, 3); /* prints out the iteration info */
  //    HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

  //    /* Now setup and solve! */
  //    HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
  //    HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

  //    /* Run info - needed logging turned on */
  //    HYPRE_PCGGetNumIterations(solver, &num_iterations);
  //    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  //    printf("\n");
  //    printf("Iterations = %d\n", num_iterations);
  //    printf("Final Relative Residual Norm = %e\n", final_res_norm);
  //    printf("\n");

  //    /* Destroy solver */
  //    HYPRE_ParCSRPCGDestroy(solver);
  //  }
  //  /* PCG with AMG preconditioner */
  //  else if (solver_id == 1)
  //  {
  //    int num_iterations;
  //    double final_res_norm;

  //    /* Create solver */
  //    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

  //    /* Set some parameters (See Reference Manual for more parameters) */
  //    HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
  //    HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
  //    HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  //    HYPRE_PCGSetPrintLevel(solver, 1); /* print solve info */
  //    HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

  //    /* Now set up the AMG preconditioner and specify any parameters */
  //    HYPRE_BoomerAMGCreate(&precond);
  //    HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
  //    HYPRE_BoomerAMGSetCoarsenType(precond, 6);
  //    HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
  //    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
  //    HYPRE_BoomerAMGSetTol(precond, 1e-15); /* conv. tolerance zero */
  //    HYPRE_BoomerAMGSetMaxIter(precond, 100); /* do only one iteration! */

  //    /* Set the PCG preconditioner */
  //    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
  //                        (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

  //    /* Now setup and solve! */
  //    HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
  //    HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

  //    /* Run info - needed logging turned on */
  //    HYPRE_PCGGetNumIterations(solver, &num_iterations);
  //    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  //    HYPRE_ParVector par_residualVector;
  //    HYPRE_PCGGetResidual(solver,(void **)&par_residualVector);

  //    printf("\n");
  //    printf("Iterations = %d\n", num_iterations);
  //    printf("Final Relative Residual Norm = %e\n", final_res_norm);
  //    printf("Residuals\n");

  //    HYPRE_IJVectorGetValues(x,1000,matrix.rows(),xValues);
  //    HYPRE_IJVectorPrint(x,"XVector.txt");
  //    HYPRE_ParVectorPrint(par_residualVector,"ResVector.txt");

  //    double* residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

  //    for(int i = 0; i < 1000; i++)
  //    {
  //      printf("%d: Value:%f, Residual:%e\n", i,  xValues[i], residualValues[i]);
  //    }

  //    printf("\n");

  //    /* Destroy solver and preconditioner */
  //    HYPRE_ParCSRPCGDestroy(solver);
  //    HYPRE_BoomerAMGDestroy(precond);

  //  }
  //  /* PCG with Parasails Preconditioner */
  //  else if (solver_id == 8)
  //  {
  //    int    num_iterations;
  //    double final_res_norm;

  //    int      sai_max_levels = 1;
  //    double   sai_threshold = 0.1;
  //    double   sai_filter = 0.05;
  //    int      sai_sym = 1;

  //    /* Create solver */
  //    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

  //    /* Set some parameters (See Reference Manual for more parameters) */
  //    HYPRE_PCGSetMaxIter(solver, 100000); /* max iterations */
  //    HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
  //    HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  //    HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
  //    HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

  //    /* Now set up the ParaSails preconditioner and specify any parameters */
  //    HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);

  //    /* Set some parameters (See Reference Manual for more parameters) */
  //    HYPRE_ParaSailsSetParams(precond, sai_threshold, sai_max_levels);
  //    HYPRE_ParaSailsSetFilter(precond, sai_filter);
  //    HYPRE_ParaSailsSetSym(precond, sai_sym);
  //    HYPRE_ParaSailsSetLogging(precond, 3);

  //    /* Set the PCG preconditioner */
  //    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
  //                        (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);

  //    /* Now setup and solve! */
  //    HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
  //    HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);


  //    /* Run info - needed logging turned on */
  //    HYPRE_PCGGetNumIterations(solver, &num_iterations);
  //    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  //    printf("\n");
  //    printf("Iterations = %d\n", num_iterations);
  //    printf("Final Relative Residual Norm = %e\n", final_res_norm);
  //    printf("\n");


  //    /* Destory solver and preconditioner */
  //    HYPRE_ParCSRPCGDestroy(solver);
  //    HYPRE_ParaSailsDestroy(precond);
  //  }
  //  /* Flexible GMRES with  AMG Preconditioner */
  //  else if (solver_id == 61)
  //  {
  //    int    num_iterations;
  //    double final_res_norm;
  //    int    restart = 30;
  ////    int    modify = 1;


  //    /* Create solver */
  //    HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

  //    /* Set some parameters (See Reference Manual for more parameters) */
  //    HYPRE_FlexGMRESSetKDim(solver, restart);
  //    HYPRE_FlexGMRESSetMaxIter(solver, 1000); /* max iterations */
  //    HYPRE_FlexGMRESSetTol(solver, 1e-7); /* conv. tolerance */
  //    HYPRE_FlexGMRESSetPrintLevel(solver, 2); /* print solve info */
  //    HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */


  //    /* Now set up the AMG preconditioner and specify any parameters */
  //    HYPRE_BoomerAMGCreate(&precond);
  //    HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
  //    HYPRE_BoomerAMGSetCoarsenType(precond, 6);
  //    HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
  //    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
  //    HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
  //    HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

  //    /* Set the FlexGMRES preconditioner */
  //    HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
  //                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);


  //    //    if (modify)
  //    //      /* this is an optional call  - if you don't call it, hypre_FlexGMRESModifyPCDefault
  //    //        is used - which does nothing.  Otherwise, you can define your own, similar to
  //    //        the one used here */
  //    //      HYPRE_FlexGMRESSetModifyPC( solver,  (HYPRE_PtrToModifyPCFcn) hypre_FlexGMRESModifyPCAMGExample);


  //    /* Now setup and solve! */
  //    HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
  //    HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

  //    /* Run info - needed logging turned on */
  //    HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
  //    HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);


  //    printf("\n");
  //    printf("Iterations = %d\n", num_iterations);
  //    printf("Final Relative Residual Norm = %e\n", final_res_norm);
  //    printf("\n");


  //    /* Destory solver and preconditioner */
  //    HYPRE_ParCSRFlexGMRESDestroy(solver);
  //    HYPRE_BoomerAMGDestroy(precond);


  //  }

  //  HYPRE_IJVectorPrint(x,"XVector.txt");

  //  HYPRE_IJMatrixDestroy(A);
  //  HYPRE_IJVectorDestroy(b);
  //  HYPRE_IJVectorDestroy(x);
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  solutionExample();

  MPI_Finalize();
}
