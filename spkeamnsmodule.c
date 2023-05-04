#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>



PyObject* toPyJ(int n, double **input);
PyObject* toPy(int n, double **input);
void create_vec(int n, double *curr, PyObject* row);
void create_matrix_m(int n, int m, double ***matrix_ptr);

struct cord* createVector(struct cord *x);
struct vector* createEmptyCentroids(int k, int size_of_vector);
double euclideanDistance(struct cord *x1, struct cord *x2);
int getCordLength(struct vector *dp_head);
void divide_centroid(struct vector *new_centroid, int *num_of_dp_in_cluster, int k);
void check_converge(struct vector *new_head_vec, struct vector *head_vec1, int *convergence_array, int k, double e);
void addCords(struct cord *curr_cord1, struct cord *min_cord);
void create_Array(int *array, int len);
struct cord* createEmptyCord(int size_of_vector);
struct vector* copyOf(struct vector *old_centroid);
void freeSpaceCord(struct cord *Cord);
void freeSpaceVec(struct vector *vec);
struct cord* create_cord_from_arr(int size_of_vec, double *cord);
struct vector* arr_to_ll(int num_points, int size_of_vec, double **point_arr,int bool);
PyObject* kmeans_c(int k, int size_of_vec , int num_of_dp,int iter, double e, double **dp_arr, double **centroids_arr);
PyObject* kmeans(struct vector *dp_head, struct vector *old_centroids_head_vec,int k, int iter, int num_of_dp,double e, int size_of_vec);
PyObject* ll_to_arr(struct vector *centroids, int k, int size_of_vec);
PyObject* create_vec_from_arr(struct cord *cord,  int size_of_vec);

static PyObject* wrap_spk(PyObject *self, PyObject *args){
    int i,j;
    int num_of_dp, size_of_vec;
    PyObject *dp_arr_py, *res;
    PyObject *item1, *item2;
    double val;
    res = NULL;

    if(!PyArg_ParseTuple(args, "iiO", &num_of_dp, &size_of_vec,&dp_arr_py)) {
        return NULL; 
    }

    double *cords;
    double **dp_arr;
    double **wam_array;
    double **ddg_array;
    double **gl_array;
    double **output_array;
    
    cords = (double*)calloc(size_of_vec, sizeof(double));
    create_matrix_m(num_of_dp,size_of_vec,&dp_arr);
    create_matrix_m(num_of_dp,num_of_dp,&wam_array);
    create_matrix_m(num_of_dp,num_of_dp,&ddg_array);
    create_matrix_m(num_of_dp,num_of_dp,&gl_array);
    create_matrix_m(num_of_dp+1, num_of_dp, &output_array);

 
    for(i=0; i<num_of_dp; i++){
        item1 = PyList_GetItem(dp_arr_py,i);
        for(j=0; j<size_of_vec; j++){
            item2 = PyList_GetItem(item1, j);
            val = PyFloat_AsDouble(item2);
            cords[j] = val;
        }
        for(j=0; j<size_of_vec; j++){
            dp_arr[i][j] = cords[j];
        }

    }
    wam_c(num_of_dp,size_of_vec,dp_arr,wam_array);
    ddg_c(num_of_dp, wam_array, ddg_array);
    gl_c(num_of_dp,ddg_array,wam_array,gl_array);
    jacobi_c(num_of_dp, gl_array, output_array);
    
    res =  toPyJ(num_of_dp, output_array);
    return res;
     
}  

static PyObject* wrap_wam(PyObject *self, PyObject *args)
{
    

    int i,j;
    int num_of_dp, size_of_vec;
    PyObject *dp_arr_py, *res;
    PyObject *item1, *item2;
    double val;
    res= NULL;
    item1 = NULL;
    item2 = NULL;
    if(!PyArg_ParseTuple(args, "iiO", &num_of_dp, &size_of_vec,&dp_arr_py)) {
        return NULL; 
    }
    
    

    double *cords;
    double **dp_arr;
    double **output_array;
    
    cords = (double*)calloc(size_of_vec, sizeof(double));
    create_matrix_m(num_of_dp,size_of_vec,&dp_arr);
    create_matrix_m(num_of_dp,num_of_dp,&output_array);

    
    for(i=0; i<num_of_dp; i++){
        item1 = PyList_GetItem(dp_arr_py,i);
        
        for(j=0; j<size_of_vec; j++){
            
            item2 = PyList_GetItem(item1, j);
            
            val = PyFloat_AsDouble(item2);
            
            cords[j] = val;
        }
        for(j=0; j<size_of_vec; j++){
            dp_arr[i][j] = cords[j];
        }

    }
    
    wam_c(num_of_dp,size_of_vec,dp_arr,output_array);
    
    res = toPy(num_of_dp, output_array);
    return res;
}

static PyObject* wrap_ddg(PyObject *self, PyObject *args)
{
    int i,j;
    int num_of_dp, size_of_vec;
    PyObject *dp_arr_py;
    PyObject *item1, *item2;
    double val;

    if(!PyArg_ParseTuple(args, "iiO", &num_of_dp, &size_of_vec,&dp_arr_py)) {
        return NULL;
    }

    double*cords;
    double **dp_arr;
    double **wam_array;
    double **output_array;

    cords = (double*)calloc(size_of_vec, sizeof(double));
    create_matrix_m(num_of_dp,size_of_vec,&dp_arr);
    create_matrix_m(num_of_dp,num_of_dp,&wam_array);
    create_matrix_m(num_of_dp,num_of_dp,&output_array);

    for(i=0; i<num_of_dp; i++){
        item1 = PyList_GetItem(dp_arr_py,i);
        for(j=0; j<size_of_vec; j++){
            item2 = PyList_GetItem(item1, j);
            val = PyFloat_AsDouble(item2);
            cords[j] = val;
        }
        for(j=0; j<size_of_vec; j++){
            dp_arr[i][j] = cords[j];
        }

    }
    wam_c(num_of_dp,size_of_vec,dp_arr,wam_array);
    ddg_c(num_of_dp, wam_array, output_array);
    return toPy(num_of_dp, output_array);
}

static PyObject* wrap_gl(PyObject *self, PyObject *args)
{
    int i,j;
    int num_of_dp, size_of_vec;
    PyObject *dp_arr_py;
    PyObject *item1, *item2;
    double val;

    if(!PyArg_ParseTuple(args, "iiO", &num_of_dp, &size_of_vec,&dp_arr_py)) {
        return NULL; 
    }

    double *cords;
    double **dp_arr;
    double **wam_array;
    double **ddg_array;
    double **output_array;

    cords = (double*)calloc(size_of_vec, sizeof(double));
    create_matrix_m(num_of_dp,size_of_vec,&dp_arr);
    create_matrix_m(num_of_dp,num_of_dp,&wam_array);
    create_matrix_m(num_of_dp,num_of_dp,&ddg_array);
    create_matrix_m(num_of_dp,num_of_dp,&output_array);

    for(i=0; i<num_of_dp; i++){
        item1 = PyList_GetItem(dp_arr_py,i);
        for(j=0; j<size_of_vec; j++){
            item2 = PyList_GetItem(item1, j);
            val = PyFloat_AsDouble(item2);
            cords[j] = val;
        }
        for(j=0; j<size_of_vec; j++){
            dp_arr[i][j] = cords[j];
        }

    }
    wam_c(num_of_dp,size_of_vec,dp_arr,wam_array);
    ddg_c(num_of_dp, wam_array, ddg_array);
    gl_c(num_of_dp,ddg_array,wam_array,output_array);
    return toPy(num_of_dp, output_array);
}

static PyObject* wrap_jacobi(PyObject *self, PyObject *args){
    
    int n,i,j;
    PyObject *dp_arr_py, *res;
    PyObject *item1, *item2;
    double val;
    res = NULL;
    if(!PyArg_ParseTuple(args, "iO", &n ,&dp_arr_py)) {
        return NULL; 
    }
    double *row;
    double **matrix;
    double **output_array;

    row =  (double*)calloc(n, sizeof(double));
    create_matrix_m(n, n, &matrix);
    create_matrix_m(n+1, n, &output_array);



    for(i=0; i<n; i++){
        item1 = PyList_GetItem(dp_arr_py,i);
        for(j=0; j<n; j++){
            item2 = PyList_GetItem(item1, j);
            val = PyFloat_AsDouble(item2);
            row[j] = val;
        }
        for(j=0; j<n; j++){
            matrix[i][j] = row[j];
        }
    }

    jacobi_c(n, matrix, output_array);
    
    res = toPyJ(n, output_array);
    return res;
}

static PyObject* fit(PyObject *self, PyObject *args){
    int k;
    int iter;
    double e;
    PyObject *dp_arr_py;
    int num_of_dp;
    int size_of_vec;
    PyObject *centroids_arr_py;
    PyObject *item1, *item2;
    double val;
    int i, j;
    

    if(!PyArg_ParseTuple(args, "iidOOii", &k, &iter, &e, &dp_arr_py, &centroids_arr_py, &size_of_vec, &num_of_dp)) {
        return NULL; 
    }

    double **dp_arr;
    double **centroids_arr;
    double *cords;
    
    create_matrix_m(num_of_dp,size_of_vec,&dp_arr);
    create_matrix_m(k,size_of_vec, &centroids_arr);
    cords = (double*)calloc(size_of_vec, sizeof(double));


    for(i=0; i<num_of_dp; i++){
        item1 = PyList_GetItem(dp_arr_py,i);
        for(j=0; j<size_of_vec; j++){
            item2 = PyList_GetItem(item1, j);
            val = PyFloat_AsDouble(item2);
            cords[j] = val;
        }
        for(j=0; j<size_of_vec; j++){
            dp_arr[i][j] = cords[j];
        }

    }

    for(i=0; i<k; i++){
        item1 = PyList_GetItem(centroids_arr_py,i);
        for(j=0; j<size_of_vec;j++){
            item2 = PyList_GetItem(item1, j);
            val = PyFloat_AsDouble(item2);
            cords[j] = val;
        }
        for(j=0; j<size_of_vec; j++) {
            centroids_arr[i][j] = cords[j];
        }
    }



    return kmeans_c( k, size_of_vec, num_of_dp, iter, e, dp_arr, centroids_arr); 

}

static PyMethodDef spkmeansMethods[] = {
{"wrap_spk",
(PyCFunction) wrap_spk,
METH_VARARGS,
PyDoc_STR("calcs the spkmeans algo from given data points")
},

{"wrap_wam",
(PyCFunction) wrap_wam,
METH_VARARGS,
PyDoc_STR("calculates the wheighted adjecency matrix for given data points using the exp squared euclidean distance")},

{"wrap_ddg",
(PyCFunction) wrap_ddg,
METH_VARARGS,
PyDoc_STR("calcs the diagoinal degree matrix of the WAM of given data points")},

{"wrap_gl",
(PyCFunction) wrap_gl,
METH_VARARGS,
PyDoc_STR("calcs the Graph Laplacian for given data points")},

{"wrap_jacobi",
(PyCFunction) wrap_jacobi,
METH_VARARGS,
PyDoc_STR("calc the eigenvalue and eigenvectors of a given matrix")},

{"fit",
(PyCFunction) fit,
METH_VARARGS,
PyDoc_STR("runs the kmean algo from given data points, and centroids")
},

{NULL, NULL , 0, NULL}
};

static struct PyModuleDef spkmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "myspkmeans",
    NULL,
    -1,
    spkmeansMethods
};

PyMODINIT_FUNC PyInit_myspkmeans(void){
    PyObject *m;
    m = PyModule_Create(&spkmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

PyObject* toPyJ(int n, double **input){
    int i;
    PyObject *row, *output;
    output = PyList_New(n+1);
    for(i=0; i<n+1; i++){
        row = PyList_New(n);
        create_vec(n, input[i], row);
        
        PyList_SetItem(output, i, row);
    
        }
    return output;
}

PyObject* toPy(int n, double **input){
    int i;
    PyObject *row, *output;
    output = PyList_New(n);
    
    for(i=0; i<n; i++){
        row = PyList_New(n);

        create_vec(n, input[i], row);
        PyList_SetItem(output, i, row);

    }
    return output;
}


void create_vec(int n, double *curr, PyObject* row){
    int i;
    PyObject* py_double;
    for(i=0; i<n; i++){
        py_double = Py_BuildValue("d", curr[i]);
        PyList_SetItem(row, i, py_double);
    }
}

void create_matrix_m(int n, int m, double ***matrix_ptr){
    int i;
    double **matrix = (double**)calloc(n, sizeof(double*));
    for(i=0;i<n;i++){
        matrix[i] = (double*)calloc(m, sizeof(double));
    }
    *matrix_ptr = matrix;
}







PyObject* kmeans_c(int k, int size_of_vec , int num_of_dp,int iter, double e, double **dp_arr, double **centroids_arr)
{
    struct vector *dp_head;
    struct vector *centroids_head_vec;

    /*create the linked lists from the array*/

    dp_head = arr_to_ll(num_of_dp, size_of_vec, dp_arr, 0);
    centroids_head_vec = arr_to_ll(k, size_of_vec, centroids_arr, 1);

    return kmeans(dp_head, centroids_head_vec, k, iter, num_of_dp, e, size_of_vec);
    }

PyObject* kmeans(struct vector *dp_head, struct vector *old_centroids_head_vec,int k, int iter, int num_of_dp,double e, int size_of_vec){

    int cord_length;

    int *convergence_array; /* the array of the convergence*/

    int *num_of_dp_in_cluster; /* the array of the division*/

    int convergenceCNT;/* = 0; // for first while loop*/
    int iterCNT;/* = 0;*/


    struct vector *head_vec1, *curr_vec1; /* points to data-point struct*/
    struct cord  *curr_cord1;/**head_cord1,*/


    struct vector *head_vec2, *curr_vec2; /* points to old-centroids struct*/
    struct cord *curr_cord2;/**head_cord2,*/


    struct vector *new_head_vec, *new_curr_vec; /* the new centorids struct*/
    struct cord *new_curr_cord, *min_cord; /* *new_head_cord uninitiallized,  min_cord - indicates the cord in new centroid*/

    PyObject* final_centroids_py;

    cord_length = getCordLength(dp_head);
    convergence_array = calloc(num_of_dp , (sizeof(int)));
    create_Array(convergence_array, num_of_dp);
    num_of_dp_in_cluster = calloc(k , (sizeof(int)));
    create_Array(num_of_dp_in_cluster, k);
    convergenceCNT = 0;
    iterCNT = 0;
    head_vec2 = old_centroids_head_vec;
    curr_vec2 = head_vec2;
    curr_cord2 = head_vec2->cords;

    new_head_vec = NULL;
    min_cord = NULL;


    while(convergenceCNT<k && iterCNT<iter){ /*while not all centroids converged and the iteration is smaller than iter*/
        int cnt;
        int l;
        create_Array(num_of_dp_in_cluster, k); /* resets the nuber of dp in cluster*/

        head_vec1 = dp_head;
        curr_vec1 = head_vec1;
        curr_cord1 = head_vec1->cords;
        if(iterCNT!=0){
            freeSpaceVec(new_head_vec);
        }

        new_head_vec = createEmptyCentroids(k, cord_length);

        while(curr_vec1->cords != NULL){ /* iterates over Data Points*/
            int i;/*=0; // to where in new centroids we insert*/
            int j;
            double distance;
            double min_distance;

            new_curr_vec = new_head_vec;
            new_curr_cord = new_head_vec->cords;
            curr_vec2 = head_vec2;
            curr_cord2 = curr_vec2->cords;

            i=0;
            j=0;
            min_distance = -1;
            /*min_cord = new_head_cord;*/

            while(curr_vec2 != NULL){ /* iterates over old centroids*/
                distance = euclideanDistance(curr_cord1, curr_cord2);
                if((min_distance ==-1) || (distance < min_distance)){
                    min_cord = new_curr_cord;
                    j = i;
                    min_distance = distance;
                }

                i++;
                if((curr_vec2->next == NULL)){
                    break;
                }
                curr_vec2 = curr_vec2->next;
                curr_cord2 = curr_vec2->cords;
                if(new_curr_vec->next != NULL){
                    new_curr_vec = new_curr_vec->next;
                    new_curr_cord = new_curr_vec->cords;}

            }
            num_of_dp_in_cluster[j] += 1;
            addCords(curr_cord1,min_cord);
            if(curr_vec1->next!= NULL){
                if(curr_vec1->next->cords == NULL){
                    break;
                }
            }
            curr_vec1 = curr_vec1->next;
            curr_cord1 = curr_vec1->cords;

        }
        divide_centroid(new_head_vec, num_of_dp_in_cluster, k);
        check_converge(new_head_vec, head_vec1, convergence_array, k, e);
        cnt = 0;

        /*sums the number of convergd centroids*/
        for (l=0; l<k; l++){
            if(convergence_array[l]==1){
                cnt+=1;
            }
        }
        convergenceCNT = cnt;
        iterCNT +=1;

        /*frees the old centroid Linked list*/
        freeSpaceVec(head_vec2);
        curr_vec2 = copyOf(new_head_vec);
        head_vec2 = curr_vec2;



    }

    /*creates the final centroids PyObjectArray from the Linked List*/
    final_centroids_py = ll_to_arr(new_head_vec, k, size_of_vec);

    /*free all allocated memory*/
    freeSpaceVec(head_vec2);
    freeSpaceVec(dp_head);
    freeSpaceVec(new_head_vec);
    free(convergence_array);
    free(num_of_dp_in_cluster);

    return final_centroids_py;

}

/*calcs the euclidean distance of two given vectors*/
double euclideanDistance(struct cord *x1, struct cord *x2){

    double d = 0;
    while (x1 != NULL){
        d = d + (x1->value - x2->value)*(x1->value - x2->value);

        if(x1->next ==NULL){
            break;
        }
        x1 = x1->next;
        x2 = x2->next;
    }
    d = pow(d,0.5);
    return d;
}

/*create a linked list of k vectors with cords = 0*/
struct vector* createEmptyCentroids(int k, int size_of_vector){
    struct vector *vec, *curr_vec;
    int i;
    vec = calloc(1,sizeof(struct vector));
    curr_vec = vec;
    vec->next = NULL;

    for(i = 0; i<k; i++){
        curr_vec->cords = createEmptyCord(size_of_vector);
        if(i != (k-1)){
            curr_vec->next = calloc(1,sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
        }
    }
    return vec;
}

/*returns a copy of a cord*/
struct cord* createVector(struct cord *x){
    struct cord *res;
    struct cord *curr1, *curr2;
    res = calloc(1,sizeof(struct cord));
    curr1 = res;
    curr2 = x;
    while(curr2->next != NULL){
        curr1->value = curr2->value;
        curr1->next = calloc(1,sizeof(struct cord));
        curr1 = curr1->next;
        curr2 = curr2->next;
    }
    curr1->value = curr2->value;
    curr1->next = NULL;
    return res;
}

/*return the number of coordiates in a vector*/
int getCordLength(struct vector *dp_head){
    int res = 0;
    struct cord *curr_cord;
    curr_cord = dp_head->cords;
    while(curr_cord != NULL){
        res++;
        curr_cord = curr_cord->next;
    }
    return res;
}

/*for every centroids in centroids linked list - divides every cooridate by the number of dp in the current cluster*/
void divide_centroid(struct vector *new_centroid, int *num_of_dp_in_cluster, int k){
    struct vector *curr_vec;
    struct cord *curr_cord;
    int i;
    curr_vec = new_centroid;
    curr_cord = new_centroid->cords;
    for(i=0;i<k; i++){
        while(curr_cord != NULL){
            curr_cord->value = (curr_cord->value)/(num_of_dp_in_cluster[i]);
            curr_cord = curr_cord->next;
        }

        if(curr_vec->next == NULL){
            break;
        }
        curr_vec = curr_vec->next;
        curr_cord = curr_vec->cords;
    }
}

/*updates convergenec array by the convergence of each centroid*/
void check_converge(struct vector *new_head_vec, struct vector *head_vec1, int *convergence_array, int k, double e){
    int i;
    double d;
    struct vector *curr_vec1, *curr_vec2;
    curr_vec1 = new_head_vec;
    curr_vec2 = head_vec1;
    

    for(i=0; i<k; i++){
        if(convergence_array[i]==1){
            curr_vec1 = curr_vec1->next;
            curr_vec2 = curr_vec2->next;
            continue;
        }
        d = euclideanDistance(curr_vec1->cords,curr_vec2->cords);

        if( d< e){
            convergence_array[i] = 1;
        }
        curr_vec1 = curr_vec1->next;
        curr_vec2 = curr_vec2->next;
    }
}

/*add the coordinates of curr_cord1 to the coorinates of min_cord in place*/
void addCords(struct cord *curr_cord1, struct cord *min_cord){
    struct cord *curr_min_cord, *curr_cord;
    curr_cord = curr_cord1;
    curr_min_cord = min_cord;
    while(curr_cord != NULL){
        curr_min_cord->value += curr_cord->value;
        if((curr_cord->next == NULL) || (curr_min_cord->next == NULL)){
            break;
        }
        curr_cord = curr_cord->next;
        curr_min_cord = curr_min_cord->next;
    }
}

/*creates an array of zeros in the length of len*/
void create_Array(int *array, int len){
    int i;
    for( i=0; i<len; i++){
        array[i] = 0;
    }
}

/*returns a cord with zeros as the value of each coordinate*/
struct cord* createEmptyCord(int size_of_vector){
    struct cord *head_cord, *curr_cord;
    int i;
    head_cord = calloc(1,sizeof(struct cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;
    for(i=1;i<size_of_vector; i++){
        curr_cord->value = 0;
        if(i!=size_of_vector){
        curr_cord->next = calloc(1,sizeof(struct cord));
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
        }
    }
    return head_cord;
}

/*return a copy of an given linked list*/
struct vector* copyOf(struct vector *old_centroid){
    struct vector *new_centroids_head, *curr_vec1, *curr_vec2;

    new_centroids_head = calloc(1,sizeof(struct vector));
    new_centroids_head->next = NULL;
    new_centroids_head->cords = NULL;
    curr_vec1 = new_centroids_head;
    curr_vec2 = old_centroid;

    while(curr_vec2 != NULL){
        curr_vec1->cords = createVector(curr_vec2->cords);
        if(curr_vec2->next != NULL){
            curr_vec1->next = calloc(1,sizeof(struct vector));
            curr_vec1 = curr_vec1->next;
            curr_vec1->next = NULL;
        }
        if(curr_vec2->next == NULL){
            break;
        }
        curr_vec2 = curr_vec2->next;
    }
    return new_centroids_head;
}

/*free the allocated memory of a given struct vector*/
void freeSpaceVec(struct vector *vec){
    struct vector *curr, *tmp;
    curr = vec;
    while (curr != NULL)
    {
        freeSpaceCord(curr->cords);
        curr = curr->next;
    }
    curr = vec;
    tmp = vec;
    while (curr != NULL){
        tmp = curr;
        curr = curr->next;
        free(tmp);
    }
}

/*free the allocated memory of a given struct cord*/
void freeSpaceCord(struct cord *cord){
    struct cord *curr, *tmp;
    curr = cord;
    while(curr != NULL){
        tmp = curr;
        curr = curr->next;
        free(tmp);
    }
}



/*returns a given struct vector as a PyObject array of arrays*/
PyObject* ll_to_arr(struct vector *centroids, int k, int size_of_vec){
    int i;
    struct vector *curr_vec;
    PyObject* res = PyList_New(k);
    PyObject* cord = PyList_New(size_of_vec);
    curr_vec = centroids;
    for(i=0; i<k; i++){
        cord = create_vec_from_arr(curr_vec->cords, size_of_vec);
        PyList_SetItem(res, i, cord);
        if(i == k-1){
            break;
        }
        curr_vec = curr_vec->next;
    }
    return res;
}

/*returns a given struct cord as a PyObject array*/
PyObject* create_vec_from_arr(struct cord *cord,  int size_of_vec){
    struct cord* curr_cord;
    int i;
    PyObject* python_double;
    PyObject* res = PyList_New(size_of_vec);
    curr_cord = cord;
    for(i=0; i<size_of_vec; i++){
        python_double = Py_BuildValue("d",curr_cord->value);
        PyList_SetItem(res, i, python_double);
        if(i == size_of_vec-1){
            break;
        }
        curr_cord = curr_cord->next;
    }
    return res;
}

/*returns a given array of arrays as a struct vetor*/
struct vector* arr_to_ll(int num_points, int size_of_vec, double **point_arr,int bool){
    struct vector *head_vec, *curr_vec;
    head_vec = calloc(1,sizeof(struct vector));
    curr_vec = head_vec;
    int i;
    for(i=0; i<num_points; i++){
        curr_vec->cords = create_cord_from_arr(size_of_vec, point_arr[i]);
        if((bool) && (i == num_points-1)){
            break;
        }
        curr_vec->next = calloc(1,sizeof(struct vector));;
        curr_vec = curr_vec->next;
    }
    return head_vec;
}

/*returns a given array as a struct cord*/
struct cord* create_cord_from_arr(int size_of_vec, double *cord){
    struct cord *head_cord, *curr_cord;
    head_cord = calloc(1,sizeof(struct cord));
    curr_cord = head_cord;
    int i;
    for(i=0; i<size_of_vec; i++){
        curr_cord->value = cord[i];
        if(i == size_of_vec-1){
            break;
        }
        curr_cord->next = calloc(1,sizeof(struct cord));
        curr_cord = curr_cord->next;
    }
    return head_cord;
}

