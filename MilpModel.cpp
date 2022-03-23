// Replicate simple MILP problem solved in Quinton2020 paper - Section 4: Flexible Cyclic jobshop scheduling problem
//
#define _CRTDBG_MAP_ALLOC
#include <iostream>

#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <crtdbg.h>


int main()
{
    // number of tasks - plus 2 tasks (start and end node) 
     const int noTasks = 13; // SMAPI Data


    // number of available machines - added pseudo machines to process start and end nodes
     const int noMachines = 11; // SMAPI Data

    // number of batches considered
    int noBatch;

    // parameters that denote the best cycle time found, the best number of batches and the current cycle time
    double bestCycleT = 1.0e30;
    double cycleT;
    double bestNb = 2;

    // define the optimization sense and objective offset
    const int sense_minimization = 1;
    const int sense_maximization = -1;
    int sense = sense_maximization;
    const double offset = 0;

    // declare set E(i,j) - this should be automated using the start times of each task
    // values of E(i,j) are equal to parameter p(i,j) - distance between starting times of two tasks
    

    // startT denotes the starting time of each task in the recipe
    double startT[noTasks] = { 0, 0, 9, 15, 15.5, 16, 19.5, 20.5, 24, 25, 43, 45, 45 }; // SMAPI Data

    double E[noTasks][noTasks];
    //initialize
    for (int i = 0; i < noTasks; i++)
    {
        for (int j = 0; j < noTasks; j++)
        {
            E[i][j] = -1;
        }
    }

    E[0][1] = 0; // first task always starts with start node
    E[noTasks-1][0] = 0; // end node always connected to start node
    E[noTasks-2][noTasks-1] = 0; // last task always starts with end node

    for (int i = 1; i < noTasks-2; i++)
    {
        for (int j = i+1; j < noTasks-1; j++)
        {
            E[i][j] = startT[j] - startT[i];
        }
    }

    // SE(i)=1 when a task is a start or an end node
    //int* SE = new int [noTasks];
    int SE[noTasks];
    
    for (int i = 0; i < noTasks; i++)
    {
            SE[i] = 0;
    }
    SE[0] = 1;
    SE[noTasks-1] = 1;

    // declare set D(i,j) - is 1 when tasks i and j can be processed in the same machine
  
    int D[noTasks][noTasks];

    // declare set H(i,j) - is the height between two nodes. 
    // regarding this problem all values are zero except the one between the end and the start node, which denotes the work in progress
    
    int H[noTasks][noTasks];
    H[noTasks-1][0] = 4;

   
    // Data SMAPI
    double R[noTasks][noMachines] = { {-1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   0,  -1},
                                      {16,  16,  16,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1},
                                      {-1,  -1,  -1,  -1,	3,   3,  -1,  -1,	3,  -1,  -1},
                                      {-1,  -1,  -1,  -1,  -1,  -1,	  2,   2,  -1,  -1,  -1},
                                      {-1,  -1,  -1,  -1,	3,   3,  -1,  -1,	3,  -1,  -1},
                                      {-1,  -1,  -1,  -1, 4.5, 4.5,  -1,  -1, 4.5,  -1,  -1},
                                      {-1,  -1,  -1,  -1,	2,   2,  -1,  -1,	2,  -1,  -1},
                                      {-1,  -1,  -1,  -1, 4.5, 4.5,  -1,  -1, 4.5,  -1,  -1},
                                      {-1,  -1,  -1,  -1,	2,   2,  -1,  -1,	2,  -1,  -1},
                                      {20,  20,  20,  20,  -1,  -1,  -1,  -1,  -1,  -1,  -1},
                                      {-1,  -1,  -1,  -1,  -1,  -1,	  4,   4,  -1,  -1,  -1},
                                      {-1,  -1,  -1,  -1,  -1,  -1,	  7,   7,  -1,  -1,  -1},
                                      {-1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,	  0}
    };

    


    // declare set RR(i,j,r). It takes the value 1 when a machine r can process both tasks i and j

    double RR[noTasks][noTasks][noMachines];

    // define nodes with disjunctive connection - if D(i,j) = 1 then they i,j share same equipment
    // value of D denotes the position of the disjunctive variable K(i,j,l,ll) in the Amatrix  
    int k = 0;
    for (int i = 0; i < noTasks; i++) 
    {
        for (int j = 0; j < noTasks; j++) 
        {
            D[i][j] = 0;
            for (int r = 0; r < noMachines; r++) 
            {
                if ((R[j][r] >= 0) && (R[i][r] >= 0)) 
                {
                    RR[i][j][r] = 1;
                    D[i][j] = 1;
                }
                else
                {
                    RR[i][j][r] = 0;
                }
            }
        }
    }

    // Outer loop over the different number of batches to find minimum cycle time
    for (noBatch = 1; noBatch < 7; noBatch++) 
    {
        // initialize number of disjunctives variables K(i,j,l,ll)
        int nDisj = 0;

        // initialize Ri - denotes the number of m(i,r,l) and y(i,r,l) variables
        int nRi = 0;

        // initialize number of columns (variables) 
        int numCol = 0;

        // initialize number of rows (constraints)
        int numRow = 0;

        // initialize number of nonzero elements in the A matrix 
        int numNz = 0;

        // calculate number of columns (variables)
        // 
        // calculate number of u(i,l) variables
        int nU = noTasks * noBatch; 

        // initialize 4-dimensional array used to denote position of variables K(i,j,l,ll) in the A matrix
        int**** posK = new int***[noTasks];
        for (int i = 0; i < noTasks; i++)
        {
            posK[i] = new int**[noTasks];
            for (int j = 0; j < noTasks; j++)
            {
                posK[i][j] = new int*[noBatch];
                for (int k = 0; k < noBatch; k++)
                {
                    posK[i][j][k] = new int[noBatch];
                }
            }
        }

        // calculate number of K(i,j,l,ll) variables and populate posK(i,j,l,ll)
        k = 0;
        for (int i = 0; i < noTasks; i++)
        {
            for (int j = 0; j < noTasks; j++)
            {
                for (int l = 0; l < noBatch; l++)
                {
                    for (int ll = 0; ll < noBatch; ll++)
                    {
                        if ((D[i][j] == 1) && ((i != j) || (i == j && l != ll)))
                        {
                            k++;
                            nDisj++;
                            posK[i][j][l][ll] = k;
                        }
                    }
                }
            }
        }
 
        // initialize 3-dimensional array used to denote position of variables m(i,r,l) and y(i,r,l) in the A matrix
        k = 0;
        int*** posR = new int** [noTasks];
        for (int i = 0; i < noTasks; i++)
        {
            posR[i] = new int* [noMachines];
            for (int j = 0; j < noMachines; j++)
            {
                posR[i][j] = new int [noBatch];

            }
        }

        // calculate number of m(i,r,l) and y(i,r,l) variables and populate posR(i,r,l)
        for (int i = 0; i < noTasks; i++) 
        {
            for (int r = 0; r < noMachines; r++) 
            {
                if (R[i][r] >= 0) 
                {
                    for (int l = 0; l < noBatch; l++) 
                    {
                        k++;
                        nRi++;
                        posR[i][r][l] = k;
                    }
                }
            }
        }

        // number of columns = number of u vars (ui,l) plus number of disjunctive vars (K(i,j,l,ll)) plus number of variables y(i,r,l) and m(i,r,l) plus one (inverse of cycle time)
        numCol = nU + nDisj + 2 * nRi + 1;

        // Calculate number of rows

        //constraint 1: for all i,j in E(i,j) and i or j or both not in SE
        //constraint 2: for all i,j in E(i,j) and both i and j in SE
        //constraint 3: for all i,j,r,l,ll that are in RR(i,j,r) and i<>j OR i==j and l<>ll
        //constraint 4: for all i,l
        //constraint 5: for all i,r in R(i,r) and for all l
        //constraint 6: for all i,l and i not in SE
        //constraint 7: for all i,j,l,ll with D(i,j)=1 and (i<j OR i=j and l<ll)
        //constraint 8: for all l (but not for the last batch!)
            
        for (int i = 0; i < noTasks; i++) 
        {
            for (int j = 0; j < noTasks; j++) 
            {
                for (int r = 0; r < noMachines; r++) 
                {
                    for (int l = 0; l < noBatch; l++) 
                    {
                        if (r==0 && E[i][j] >= 0 && ((SE[i] + SE[j]) != 2 || (SE[i] == 1 && SE[j] == 1)))
                        {
                            numRow++; //constraints 1 and 2 - NOT for r, so r==0 
                        }
                        if (r==0 && j == 0)
                        {
                            numRow++; // constraint 4 - for all i,l, so r=j=-
                            if (SE[i] == 0) 
                            {
                                numRow++; // objective - for all i,l, so r=j=-
                                numRow++; // constraint 6 - for all i,l, so r=j=-
                            }
                            else if (i == 0 && l < noBatch - 1)
                            {
                                numRow++; // constraint 8 - for all l so i=j=r=0
                            }
                        }
                        if (R[i][r] >= 0 && j == 0) 
                        {
                            numRow++; //constraint 5 - for i and r NOT for j, that is why we use j=1
                        }
                        for (int ll = 0; ll < noBatch; ll++) 
                        {
                            if (RR[i][j][r] == 1 && ((i != j) || (i == j && l != ll))) 
                            {
                                numRow++; // constraint 3
                            }
                            if (D[i][j] == 1 && (i < j || (i == j && l < ll)) && r == 0)
                            {
                                numRow++; //constraint 7 - for i,j,l,ll, NOT for r, therefore additional condition r=0
                            }
                        }
                    }
                }
            }
        }

        // initialize A matrix
        double** A = new double* [numRow];
        for (int i = 0; i < numRow; i++)
        {
            A[i] = new double [numCol];
            for (int j = 0; j < numCol; j++)
            {
                A[i][j] = 0;
            }
        }

        double* row_lower = new double [numRow];
        double* row_upper = new double[numRow];

        // Generate A matrix
        k = 0;
        

        // Objective: T - SUM(r,m(i,r,l)/pt(i,r)) <= 0
        //Constraint 8: u('0',l+1) - u('0',l) = 1/noBatch;
        for (int i = 0; i < noTasks; i++)
        {    
            for (int j = 0; j < noTasks; j++)
            {
                for (int l = 0; l < noBatch; l++)
                {

                    if (E[i][j] >= 0 && (SE[i] + SE[j]) != 2) // constraint 1
                    {
                        A[k][i * noBatch + l] = -1; 
                        A[k][j * noBatch + l] = 1;
                        for (int r = 0; r < noMachines; r++) 
                        {
                            if (R[i][r] >= 0) 
                            {
                               A[k][nU + nDisj + nRi + posR[i][r][l]] = -E[i][j];
                            }
                        }
                        row_upper[k] = 0;
                        row_lower[k] = 0;
                        k++;
                    }
                    else if (E[i][j] >= 0 && (SE[i] == 1 && SE[j] == 1)) // constraint 2
                    {
                        A[k][i * noBatch + l] = -1; // u(i,l)-> -1
                        A[k][j * noBatch + l] = 1; // u(j,l) -> 1
                        for (int r = 0; r < noMachines; r++) {
                            if (R[i][r] >= 0) {
                                A[k][nU + nDisj + nRi + posR[i][r][l]] = -E[i][j];
                            }
                        }
                        row_upper[k] = 1.0e30;
                        row_lower[k] = -H[i][j];
                        k++;
                    }

                    if (i == 0 && j==0 && l < noBatch - 1) // constraint 8
                    {
                        A[k][l + 1] = 1;
                        A[k][l] = -1;
                        row_upper[k] = 1.0 / noBatch;
                        row_lower[k] = 1.0 / noBatch;
                        k++;
                    }
                    if (j==0 && SE[i] == 0) // objective
                    {
                        A[k][nU] = 1;  // T -> 1
                        for (int r = 0; r < noMachines; r++)
                        {
                            if (R[i][r] >= 0)
                            {
                                A[k][nU + nDisj + posR[i][r][l]] = -1 / R[i][r];
                            }
                        }
                        row_upper[k] = 0;
                        row_lower[k] = -1.0e30;
                        k++;
                    }
                }
            }
        }

        // Constraint 3: BM*(2-m(i,r,l)-m(j,r,ll) + u(j,ll) + K(i,j,l,ll) - U(i,l) - y(i,r,l)*pt(i) >= 0
        for (int i = 0; i < noTasks; i++) 
        {
            for (int j = 0; j < noTasks; j++) 
            {
                for (int r = 0; r < noMachines; r++) 
                {
                    for (int l = 0; l < noBatch; l++) 
                    {
                        for (int ll = 0; ll < noBatch; ll++) 
                        {
                            if (RR[i][j][r] == 1 && (i != j || (i == j && l != ll))) 
                            {
                                A[k][i * noBatch + l] = -1; // u(i,l)-> -1
                                A[k][j * noBatch + ll] = 1; // u(j,ll) -> 1
                                A[k][nU + nDisj + nRi + posR[i][r][l]] = -R[i][r]; // y(i,r,l) -> -R(i,r)
                                A[k][nU + nDisj + posR[i][r][l]] = -100; // m(i,r,l) -> -BM
                                A[k][nU + nDisj + posR[j][r][ll]] = -100; // m(j,r,ll) -> -BM
                                A[k][nU + posK[i][j][l][ll]] = 1; // K(i,j,l,ll) -> 1
                                row_upper[k] = 1.0e30;
                                row_lower[k] = -200;
                                k++;
                            }
                        }
                    }
                }
            }
        }


        //Constraint 4: SUM(r,y(i,r,l)) - T = 0
        for (int i = 0; i < noTasks; i++) 
        {
            for (int l = 0; l < noBatch; l++) 
            {
                for (int r = 0; r < noMachines; r++) 
                {
                    if (R[i][r] >= 0) 
                    {
                        A[k][nU + nDisj + nRi + posR[i][r][l]] = 1;
                    }
                }
                A[k][nU] = -1; //T-> -1
                row_upper[k] = 0;
                row_lower[k] = 0;
                k++;
            }
        }



        //Constraint 5: y(i,r,l) - BM*m(i,r,l) <=0
        for (int i = 0; i < noTasks; i++) 
        {
            for (int r = 0; r < noMachines; r++) 
            {
                if (R[i][r] >= 0) 
                {
                    for (int l = 0; l < noBatch; l++) 
                    {
                        A[k][nU + nDisj + nRi + posR[i][r][l]] = 1; // y(i,r,l) -> 1
                        A[k][nU + nDisj + posR[i][r][l]] = -100; // m(i,r,l) -> -BM
                        row_upper[k] = 0;
                        row_lower[k] = -1.0e30;
                        k++;
                    }
                }
            }
        }


        //Constraint 6: SUM(r,m(i,r,l)) = 1
        for (int i = 0; i < noTasks; i++) 
        {
            if (SE[i] == 0) 
            {
                for (int l = 0; l < noBatch; l++) 
                {
                    for (int r = 0; r < noMachines; r++) 
                    {
                        if (R[i][r] >= 0) 
                        {
                            A[k][nU + nDisj + posR[i][r][l]] = 1; // m(i,r,l) -> 1
                        }
                    }
                    row_upper[k] = 1;
                    row_lower[k] = 1;
                    k++;
                }
            }
        }

        //Constraint 7: K(i,j,l,ll) + K(j,i,ll,l) =1
        for (int i = 0; i < noTasks; i++) 
        {
            for (int j = i; j < noTasks; j++) 
            {
                if (D[i][j] == 1) 
                {
                    for (int l = 0; l < noBatch; l++) 
                    {
                        for (int ll = 0; ll < noBatch; ll++) 
                        {
                            if (i < j || (i == j && l < ll)) 
                            {
                                A[k][nU + posK[i][j][l][ll]] = 1; // K(i,j,l,ll) -> 1
                                A[k][nU + posK[j][i][ll][l]] = 1; // K(j,i,ll,l) -> 1

                                row_upper[k] = 1;
                                row_lower[k] = 1;
                                k++;
                            }
                        }
                    }
                }
            }
        }

        // calculate number of nonzero elements in the A matrix
        for (int i = 0; i < numRow; i++) 
        {
            for (int j = 0; j < numCol; j++) 
            {
               // if (A[i * numCol + j] != 0)
                if (A[i][j] != 0)
                {
                    numNz++;
                }
            }
        }

        // Define the column costs - only cycle time in the objective function 
       // min f = a // So all coefficients are zero except of col_cost[0] = 1

        double* col_cost = new double[numCol];
        for (int i = 0; i < numCol; i++) 
        {
            if (i == nU) 
            {
                col_cost[i] = 1;
            }
            else 
            {
                col_cost[i] = 0;
            }
        }

        // Define the variable lower bounds
        double* col_lower = new double[numCol];

        for (int i = 0; i < numCol; i++) 
        {
            if (i > nU && i <= nU + nDisj) 
            {
                col_lower[i] = -1.030;
            }
            else 
            {
                col_lower[i] = 0;
            }

        }

        // Define the variable upper bounds - no upper bounds for hte variables of this optimization problem 
        // t<=min_t, u(i)<=+inf, K(i,j)<=1
        double* col_upper = new double[numCol];
        for (int i = 0; i < numCol; i++) 
        { //make this just a for statement
            if (i <= (nU + nDisj + nRi) && i > (nU + nDisj)) 
            {
                col_upper[i] = 1;
            }
            else 
            {
                col_upper[i] = 1.0e30;
            }
        }

        // Define the constraint matrix column-wise
        // * The indices of the nonnzeros in the vectors of A are stored in a_index
        //
        // * The values of the nonnzeros in the vectors of A are stored in a_value
        //
        // * The position in a_index/a_value of the index/value of the first
        // nonzero in each vector is stored in a_start
        int a_format = 1;
        int* a_start = new int[numRow];
        int* a_index = new int[numNz];
        double* a_value = new double[numNz];


        k = 0;
        int l = 0;
        for (int j = 0; j < numCol; j++) {
            bool n = true;
            for (int i = 0; i < numRow; i++) {
                //if (A[i * numCol + j] != 0 and n)
                if (A[i][j] != 0 and n)
                {
                    a_start[l] = k;
                    n = false;
                    l++;
                }
                //if (A[i * numCol + j] != 0)
                if (A[i][j] != 0)
                {
                 //   a_value[k] = A[i * numCol + j];
                    a_value[k] = A[i][j];
                    a_index[k] = i;
                    k++;
                }
            }
        }

        double objective_value;

        double* col_value = new double[numCol];
        double* col_dual = new double[numCol];
        double* row_value = new double[numRow];
        double* row_dual = new double[numRow];
        int* col_basis_status = new int[numCol];
        int* row_basis_status = new int[numRow];

        int model_status;
        int run_status;

        int* integrality = new int[numCol];
        for (int i = 0; i < numCol; i++) {
            integrality[i] = 0;
        }
        for (int i = nU + 1; i <= nU + nDisj + nRi; i++) {
            integrality[i] = 1;
        }


        run_status = Highs_mipCall(numCol, numRow, numNz, a_format,
            sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
            a_start, a_index, a_value,
            integrality,
            col_value, row_value,
            &model_status);


        // The run must be successful, and the model status optimal
        assert(run_status == 0);
        assert(model_status == 7);

        printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

        objective_value = offset;
        // Report the column primal and dual values, and basis status
        for (int i = 0; i < numCol; i++) {
            printf("Col%d = %lf; dual = %lf; status = %d\n", i, col_value[i], col_dual[i], col_basis_status[i]);
            objective_value += col_value[i] * col_cost[i];
        }


        // Report the row primal and dual values, and basis status
        for (int i = 0; i < numRow; i++) {
            printf("Row%d = %lf; dual = %lf; status = %d\n", i, row_value[i], row_dual[i], row_basis_status[i]);
        }
        printf("Optimal objective value = %g\n", objective_value);
        cycleT = 1 / (objective_value * noBatch);
        if (cycleT > bestCycleT) {
            break;
        }
        else {
            bestCycleT = cycleT;
            bestNb = noBatch;
        }

        // delete dynamic arrays
        // delete posK
        for (int i = 0; i < noTasks; i++)
        {
            for (int j = 0; j < noTasks; j++)
            {
                for (int k = 0; k < noBatch; k++)
                {
                    delete[] posK[i][j][k];
                }
                delete[] posK[i][j];
            }
            delete[] posK[i];
        }
        delete[] posK;

        // delete posR
        for (int i = 0; i < noTasks; i++)
        {
            for (int j = 0; j < noMachines; j++)
            {
                delete[] posR[i][j];
            }
            delete[] posR[i];
        }
        delete[] posR;

        // delete A
        for (int i = 0; i < numRow; i++)
        {
            delete[] A[i];
        }
        delete[] A;

        // delete rest
        delete[] row_lower;
        delete[] row_upper;
        delete[] col_cost;
        delete[] col_lower;
        delete[] col_upper;
        delete[] a_start;
        delete[] a_index;
        delete[] a_value;
        delete[] col_value;
        delete[] col_dual;
        delete[] row_value;
        delete[] row_dual;
        delete[] col_basis_status;
        delete[] row_basis_status;
        delete[] integrality;  
    }

    printf("Optimal cycle time = %g\n", bestCycleT);
    printf("Optimal number of batches = %g\n", bestNb);

    _CrtDumpMemoryLeaks();

}