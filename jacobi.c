#include "spkmeans.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>



double eps = 0.00001;
int max_iter=100;

void max_off_diag(int n, double **A, int pivot[2]);/*A is n*n */
int sign(double a);
void rotation(int n,double s,double c, int pivot [2], double **A, double **new_A); /* A and new_A is n*n */
void ID(int n, double **res); /*res is n*n*/
void dot(int n, double **a , double **b);/*n*n matrix */
double sum_dot(int n, double *row, double *col); /* arrays of size n*/
void calc_col(int j, int n, double **matrix, double *col ); /*matrix size n*n, col is array size n*/
int check_convergeJ(int n, double **A , double **new_A); /* A and new_A is n*n */
void add_eigenV(int n, double **A , double **V, double **res ); /* A and V is n*n res is n+1*n */
void create_matrix(int n, int m, double ***matrix_ptr);



void jacobi_c(int n, double **A , double **res){
    int iter,i,j;
    double **V = NULL;
    double **P = NULL;
    double **new_A = NULL;
    create_matrix(n,n,&V); /*added this line*/
    create_matrix(n,n,&P); /*added this line*/
    create_matrix(n,n,&new_A); /*added this line*/
    
    ID(n,V);
    for(iter=0;iter<max_iter;iter++){
        double theta, t, c, s;
        int pivot[2] ;
        
        max_off_diag(n,A, pivot);
        if (fabs(A[pivot[0]][pivot[1]])<=eps){
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    new_A[i][j] = A[i][j];
            }
        }
            break;
        }
        
        theta = (A[pivot[1]][pivot[1]]-A[pivot[0]][pivot[0]])/(2*A[pivot[0]][pivot[1]]);
        t = ((sign(theta))/(fabs(theta)+pow(pow(theta,2)+1,0.5)));
        c = 1/ pow(pow(t,2)+1,0.5);
        s = t*c;

        rotation(n , s, c ,pivot, A, new_A);

        ID(n,P);
        P[pivot[0]][pivot[0]] = c;
        P[pivot[1]][pivot[1]] = c;
        P[pivot[1]][pivot[0]] = -s;
        P[pivot[0]][pivot[1]] = s;
        
        dot(n,V,P);  

        if (check_convergeJ(n,A,new_A)){
            break;}

        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                A[i][j] = new_A[i][j];
            }
        }

              
    }
    add_eigenV(n,new_A,V, res);

    free_matrix(new_A,n);
    free_matrix(V,n);
    free_matrix(P,n);
    
    
}


void max_off_diag(int n, double **A, int pivot[2]){
    int i,j;
    double max;
    max = -1;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                continue;
            }
            if(fabs(A[i][j]) > max){
                max = fabs(A[i][j]);
                pivot[0] = i;
                pivot[1] = j;
            }
        }
    }
    
}

int sign(double a){
    if(a<0){
        return -1;
    }
    else{
        return 1;
    }
}

void rotation(int n,double s, double c, int pivot [2], double **A , double **new_A){
    int r,i,j,k,l;
    double **res = NULL;
    i = pivot[0];
    j = pivot[1];
    create_matrix(n,n,&res); /*added this line*/
    

    for (k=0;k<n;k++){
        for (l=0;l<n;l++){
            res[k][l] = A[k][l];
        }
    }
 
    for(r=0; r<n; r++){                        /*#1 & #2*/
        if (r != i && r != j){
            res[r][i] = c*A[r][i] - s*A[r][j];
            res[i][r] = res[r][i];

            res[r][j] = c*A[r][j] + s*A[r][i];
            res[j][r] = res[r][j];
        }
    }
    

    res[i][i] = pow(c,2)*A[i][i] + pow(s,2)*A[j][j]-2*s*c*A[i][j];   /*#3*/
    res[j][j] = pow(s,2)*A[i][i] + pow(c,2)*A[j][j] + 2*c*s*A[i][j]; /*#4*/ 

    
    
    res[i][j] = 0;                                                   /*#5*/
    res[j][i] = 0;
    for(l=0; l<n; l++){
        for(k=0; k<n; k++){
            new_A[l][k] = res[l][k];
        }
    }
    free_matrix(res,n);
}

void ID(int n, double **res){
    int i,j;
    for (i=0; i<n;i++){
        for(j=0;j<n;j++){
            if (i == j){
                res[i][j] = 1;
            }
            else{res[i][j] = 0;}
        }
    }
}

void dot(int n, double **a , double **b){
    double **c = NULL;
    int i,j;
    double *col = NULL;
    col = (double*)calloc(n, sizeof(double));
    create_matrix(n,n,&c); /*added this line*/

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            
            calc_col(j, n, b, col);
            
            c[i][j] = sum_dot(n,a[i],col);
            
        }
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            a[i][j] = c[i][j];
        }
    }
    free_matrix(c,n);
    free(col);
}

double sum_dot(int n, double *row, double *col){
    double sum;
    int i;
    sum =0;
    for(i=0;i<n;i++){
        sum = sum + row[i]*col[i];
    }
    return sum;
    
}

void calc_col(int j, int n, double **matrix, double *col){
    int i;
    for (i=0; i<n; i++){
        col[i] = matrix[i][j];
    }
}

int check_convergeJ(int n, double **A , double **new_A){
    int i,j;
    double d1,d2;
    d1 = 0;
    d2 = 0;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            if(i==j){continue;}
            else{
                d1 = d1 + pow(A[i][j],2);
            }
        }
    }
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            if(i==j){continue;}
            else{
                d2 = d2 + pow(new_A[i][j],2);
            }
        }
    }
    if (fabs(d1-d2)<=eps){return 1;}
    return 0;
}

void add_eigenV(int n, double **A , double **V, double **res){ /*res is n+1*n  */
    
    int i,j;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            res[i][j] = V[i][j];
        }
    }
    for (i=0;i<n;i++){
        res[n][i] = A[i][i];}
}


void create_matrix(int n, int m, double ***matrix_ptr){
    int i;
    double **matrix = (double**)calloc(n, sizeof(double*));
    for(i=0;i<n;i++){
        matrix[i] = (double*)calloc(m, sizeof(double));
    }
    *matrix_ptr = matrix;
}

