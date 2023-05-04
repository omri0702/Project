# ifndef SPKMEANS_H_
# define SPKMEANS_H_

struct cord
{
    double value;
    struct cord *next;
};
struct vector
{
    struct vector *next;
    struct cord *cords;
};




void wam_c(int n, int m, double **arr , double **output_array);

void ddg_c(int n, double **W, double **D);

void gl_c(int n, double **D, double **W, double **output_array);

void jacobi_c(int n, double **A, double **output_array);

void free_matrix(double **matrix, int n);

# endif
