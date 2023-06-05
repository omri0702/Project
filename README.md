# General description:
In this project I implemented the Spectral Clustering Algorithm.
The algorithm is going through some steps- for a given 'n' data points as an input,
it forms weighted adjacency matrix W, compute the graph Laplacian L, 
and according to Jacobi algorithm, it calculates its eigenvalues and eigenvectors (let J be the matrix of eignvectors).
Then, we use "Eigengap Heuristic" in order to find the best 'k' to divide our data to k groups, and to reduce the dimension of the data.
With this k, we reduce J to be a n * k size matrix. Treating each row as a point in Rk, we finally cluster them into k clusters using K-means++, 
findind the best centroids that will represent the different groups of data.

# How to use it:
Make sure you cloned all the files from this repository.
As you can see, I created a library implemented in C and import to Pyton.
In this library there is "wam","ddg","gl","jacobi" and "spk" functions. 
You can use each of the functions that performes different steps in the algorithm, while "spk" perform its all.
You can use both the python implementation and C, but C is much faster.

In order to compile the C files- write in the terminal:
```shell
make
```

In order to create the module (C-py Api)- write in the terminal:
```shell
python3 setup.py build_ext --inplace
```

Now you can run it
Create your own datapoint file (lets call it input.txt) and write in the terminal:
```shell
python3 spkmeans.py spk input.txt
```

The returning values will be the indices of centroids and centroids themselves.

For any other part of the algorithm you can replace the "spk" with other func I implemented in the library.
