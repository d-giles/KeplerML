"""
DBSCAN Clustering
"""
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

def eps_est_recursive(data):
    
    # distance array containing all distances
    nbrs = NearestNeighbors(n_neighbors=int(np.ceil(.2*len(data))), algorithm='ball_tree').fit(data)
    distances, indices = nbrs.kneighbors(data)
    # Distance to 2*N/100th instead of 4th because: ... reasons
    neighbors = int(np.ceil(.02*len(data)))
    distArr = distances[:,neighbors]
    distArr.sort()
    pts = range(len(distArr))

    # The following looks for the first instance (past the mid point)
    # where the mean of the following [number] points
    # is at least (cutoff-1)*100% greater than the mean of the previous [number] points.
    # Alternatively, we could consider the variance of the points and draw conclusions from that
    
    if len(data) <= 200:
        number = 10
    else:
        number = 50
    cutoff = 1.05
    for i in range(int(np.ceil(len(pts)/2)),len(pts)-number):
        if np.mean(distArr[i+1:i+number])>=cutoff*np.mean(distArr[i-number:i-1]):
            dbEps = distArr[i]
            break

    # Estimating nneighbors by finding the number of pts. 
    # that fall w/in our determined eps for each point.

    count = np.zeros(len(pts))
    for i in pts:    
        for dist in distances[i]:
            if dist <= dbEps:
                count[i]+=1
    average = np.median(count)
    sigma = np.std(count)
    neighbors = average/2 # Divide by 2 for pts on the edges of clusters
    print("""
    Epsilon is in the neighborhood of %s, 
    with an average of %s neighbors within epsilon,
    %s neighbors in half circle (neighbors/2).
    """%(dbEps,average,neighbors))
    return dbEps,neighbors

def eps_est(data,n=4,verbose=True):
    
    # distance array containing all distances
    if verbose:print("Calculating nearest neighbor distances...")
    nbrs = NearestNeighbors(n_neighbors=int(max(n+1,100)), algorithm='ball_tree',n_jobs=-1).fit(data)
    distances, indices = nbrs.kneighbors(data)
    del nbrs
    distArr = distances[:,n]
    distArr.sort()
    del distances
    pts = range(len(distArr))

    # The following looks for the first instance (past the mid point)
    # where the mean of the following [number] points
    # is at least (cutoff-1)*100% greater than the mean of the previous [number] points.
    
    number = int(np.ceil(len(data)/500))
    cutoff = 1.05
    if verbose:print("Finding elbow...")
    for i in range(int(np.ceil(len(pts)/2)),len(pts)-number):
        if np.mean(distArr[i+1:i+number])>=cutoff*np.mean(distArr[i-number:i-1]):
            dbEps = distArr[i]
            pt=pts[i]
            break

    if verbose:
        print("""
        Epsilon is in the neighborhood of {:05.2f}.
        """.format(dbEps))
    return dbEps,distArr

def dbscan_w_outliers(data,min_n=4,check_tabby=False,verbose=True):
    # numpy array of dataframe for fit later
    X=np.array([np.array(data.loc[i]) for i in data.index])
    if verbose:print("Estimating Parameters...")
    if len(X)>10000:
        # Large datasets have presented issues where a single high density cluster 
        # leads to an epsilon of 0.0 for 4 neighbors.
        # We adjust for this by calculating epsilon for a sample of the data,
        # then we scale min_neighbors accordingly
        if verbose:print("Sampling data for parameter estimation...")
        X_sample = data.sample(n=10000)
    else:
        X_sample = data
    dbEps,distArr = eps_est(X_sample,n=min_n,verbose=verbose)
    if len(X)>10000:
        if verbose:print("Scaling density...")
        min_n = len(X)/10000*min_n
    
    if verbose:print("Clustering data with DBSCAN, eps={:05.2f},min_samples={}...".format(dbEps,min_n))
    #est = DBSCAN(eps=dbEps,min_samples=min_n,n_jobs=-1)
    #est.fit(X)
    #clusterLabels = est.labels_
    # Outlier score: distance to 4th neighbor?
    nbrs = NearestNeighbors(n_neighbors=min_n+1, algorithm='ball_tree',n_jobs=-1).fit(data)
    distances, indices = nbrs.kneighbors(data)
    del nbrs
    distArr = distances[:,min_n]

    d = {True:-1,False:0}
    clusterLabels = np.array([d[pt>dbEps] for pt in distArr])
    for i,label in enumerate(clusterLabels):
        # For all identified outliers (pts w/o enough neighbors):
        if label == -1:
            j=1
            # for the neighbors within epsilon
            while distances[i,j]<dbEps:
                # if a neighbor is labeled as part of the cluster,
                if clusterLabels[indices[i,j]] == 0:
                    # then this pt is an edge point
                    clusterLabels[i]=1
                    break
                j+=1
        
    if check_tabby:
        if data.index.str.contains('8462852').any():
            tabbyInd = list(data.index).index(data[data.index.str.contains('8462852')].index[0])
            if clusterLabels[tabbyInd] == -1:
                print("Tabby has been found to be an outlier in DBSCAN.")
            else:
                print("Tabby has NOT been found to be an outlier in DBSCAN")
        else:
            print("MISSING: Tabby is not in this data.")
                     
    numout = len(clusterLabels[clusterLabels==-1])
    numclusters = max(clusterLabels+1)
    if verbose:
        if numclusters==1:
            print("There was {:d} cluster and {:d} total outliers".format(numclusters,numout))
        else:
            print("There were {:d} clusters and {:d} total outliers".format(numclusters,numout))

    return clusterLabels

if __name__=="__main__":
    """
    If this is run as a script, the following will parse the arguments it is fed, 
    or prompt the user for input.
    
    python keplerml.py path/to/filelist path/to/fits_file_directory path/to/output_file
    """
    # fl - filelist, a txt file with file names, 1 per line
    if sys.argv[1]:
        f = sys.argv[1]
    else:
        while not f:
            f = raw_input("Input path: ")
     
    print("Reading %s..."%f)
    
    df = pd.read_csv(f,index_col=0)
    
    if sys.argv[2]:
        of = sys.argv[2]
    else:
        of = raw_input("Output path: ")
    if not of:
        print("No output path specified, saving to 'output.npy' in local folder.")
        of = 'output'
    
    np.save(of,dbscan_w_outliers(df[['tsne_x','tsne_y']]))
    print("Done.")