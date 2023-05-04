# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
#include "spkmeans.h"


double sum_row (int n, double *arr);

void ddg_c(int n, double **W, double **D){ /*W and D is n*n matrix */
    int i, j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i!=j){
                D[i][j]=0;
            }
            else{
                D[i][j]= sum_row(n, W[i]);
            }
        }
    }
}

double sum_row (int n, double *arr ){
    int i;
    double ret;
    ret=0;
    for (i=0; i<n ;i++){
        ret = ret + arr[i];
    }
    return ret;

}
