import sys
import math
import pandas as pd
import numpy as np
import myspkmeans
np.random.seed(0)


def main():
    funcs = ["spk","wam","ddg","gl","jacobi"]
    k, func, fileName = split_program_args()
    if func == 0:
        print("An Error Has Occurred")
        return
    # we need to check valid for k only before kmean++
    if func not in funcs:
        print("An Error Has Occurred")
        return
    
    df = pd.read_csv(fileName, header=None)
    n = len(df)
    #df = df.drop(df.columns[0],axis=1)
    dp = np.array(df.values.tolist())
    m = len(dp[0])
    res = None
    dp = dp.tolist()
    if func == "spk":
        res = np.array(myspkmeans.wrap_spk(n,m,dp))
        res = res.T
        sort_key = lambda row: row[-1]
        sorted_matrix = np.array(sorted(res, key=sort_key))
        sorted_matrix = sorted_matrix.T
        
            
        if k == -1:
            max = -1
            for i in range(math.floor(len(sorted_matrix[n])/2)):
                if abs(sorted_matrix[n][i]-sorted_matrix[n][i+1]) > max:
                    max = abs(sorted_matrix[n][i]-sorted_matrix[n][i+1])
                    k = i+1
        sorted_matrix = sorted_matrix[:-1]
        dp = sorted_matrix[:,:k]
        
        centroids, list_of_indx = create_centroids(dp, k)
        
        dp = dp.tolist()
        centroids = np.array(centroids).tolist()
        last_centroids = myspkmeans.fit(k, 300, 0, dp, centroids, len(dp[0]), len(dp))
        output_print_spk(last_centroids, list_of_indx,len(dp[0]))
        return
    
    elif func == "wam":
        res = myspkmeans.wrap_wam(n,m,dp)
    elif func == "ddg":
        res = myspkmeans.wrap_ddg(n,m,dp)
    elif func == "gl":
        res = myspkmeans.wrap_gl(n,m,dp)
    elif func == "jacobi":
        res = myspkmeans.wrap_jacobi(n,dp)
        
    
    output_print(res)


def create_centroids(dp,k):
    val_of_indx = [i for i in range(len(dp))]
    centroids = []
    list_of_indx = []
    centroids.append(first_choose(dp, list_of_indx))
    for i in range(k-1):
        prob_arr = calc_probability(dp, centroids)
        cen , correct_indx = choose_centroid(dp,centroids,prob_arr, val_of_indx)
        centroids.append(dp[correct_indx])
        list_of_indx.append(cen)
    return centroids, list_of_indx

def choose_centroid(dp,centroids,prob_arr,val_of_indx):
    cen = np.random.choice(val_of_indx, p = prob_arr)
    correct_indx = val_of_indx.index(cen)
    return cen,correct_indx
    
def calc_probability (dp,centroids):
    cnt = 0
    dist=[0 for i in range(len(dp))]
    for a in dp:
        min_dist = np.inf
        for b in centroids:
            curr_dist = euclideanDistance(a,b)
            if curr_dist < min_dist:
                min_dist = curr_dist
        dist[cnt] = min_dist
        cnt +=1
    total = sum(dist)
    probability_arr = [dist[i]/total for i in range(len(dist))]
    return probability_arr

def first_choose(dp, list_of_indx):

    a = np.random.choice(range(0,len(dp)))
    list_of_indx.append(a)
    return dp[a]

def euclideanDistance(x1,x2):
    d = 0
    for i in range(len(x1)):
        d += (x1[i]-x2[i])**2
    d = d**0.5
    return d    

def split_program_args():
    if (len(sys.argv)<2) or (len(sys.argv)>4):
        return 0,0,0
    argv_len = len(sys.argv)
    if argv_len == 4:
        k = int(sys.argv[1])
        func = sys.argv[2]
        fileName = sys.argv[3]
    else:
        k = -1
        func = sys.argv[1]
        fileName = sys.argv[2]

    k = int(k)
    return k, func, fileName

def output_print_spk(centroids, list_of_indx,size_of_vec):
    indx_str = ""
    for i in range(len(list_of_indx)):
        added = str(list_of_indx[i])
        indx_str += added
        if i==len(list_of_indx)-1:
            break
        indx_str += ","
    print(indx_str)

    for x in centroids:
        line = ''
        for i in range(size_of_vec):
            if i!= 0:
                line += ","
            add = x[i]
            add = "{:.4f}".format(add)
            add = str(add)
            line += add.__str__()
        print(line)

def output_print(matrix):
    if abs(len(matrix) - len(matrix[0])) == 1:
        vec = matrix[len(matrix) - 1]
        matrix = np.array(matrix)
        last_row = matrix[-1]
        matrix = np.delete(matrix, -1, axis=0)
        matrix = np.insert(matrix, 0, last_row, axis=0)

        
    for i in range(len(matrix)):
        line = ''
        for j in range(len(matrix[0])):
            if j!= 0:
                line += ","
            add = matrix[i][j]
            add = "{:.4f}".format(add)
            add = str(add)
            line += add.__str__()
        print(line)
    


if __name__ == "__main__":
    main()


