import pandas as pd
import pickle
import keplerml
import km_outliers
import db_outliers
import quarterTools as qt

def import_gen(filedir="/home/dgiles/Documents/KeplerLCs/output/",suffix="_output.p",fitsdir="/home/dgiles/Documents/KeplerLCs/fitsFiles/",out_file_ext='.coo'):
    """
    Purpose:
        Creates a function to import quarters with common suffixes in a common directory (like "_output.csv" or "_PaperSample.csv"). 
    Args:
        filedir (str) - path to common directory
        fitsdir (str) - path to fitsfile directory (containing seperate directories for each quarter).
        suffix (str) - common suffix for output files
    Returns:
        lambda QN - a function which can be called with a specific Quarter specifified as a string.
        
    Example:
        importer = importGen(filedir='./output/',suffix='_output.p')
        Q1_cluster_object = importer('Q1')
        Q1_cluster_object is the clusterOutlier object for the Quarter 1 output features.
    """
    return lambda QN: clusterOutliers(filedir+QN+suffix,fitsdir+QN+"fitsfiles",output_file=filedir+QN+out_file_ext)

class clusterOutliers(object):
    def __init__(self,feats,fitsDir,output_file='out.coo'):
        # feats may be provided as the dataframe itself or a filepath to a pickle or csv file containing a single dataframe.
        # the following interprets the form of feats given and defines self.data as the dataframe of the data.
        if type(feats) == pd.core.frame.DataFrame:
            self.data=feats
        else:
            try:
                with open(feats,'rb') as file:
                    self.data = pickle.load(file)
            except:
                try:
                    self.data = pd.read_csv(feats,index_col=0)
                except:
                    print("File format not recognized. Feature data must be stored as a dataframe in either a pickle or a csv.")
            assert type(self.data) == pd.core.frame.DataFrame, 'Feature data must be formatted as a pandas dataframe.'
        self.feats = self.data # self.feats defined as an alias for the data.
        if fitsDir[-1]=="/":
            self.fitsDir = fitsDir
        else:
            self.fitsDir = fitsDir+"/"
        self.output_file = output_file
        self.files = self.data.index # A list of all files.
        # Initializing the data and files samples with 1000 entries.
        self.sample(1000,df='self',tabby=False,replace=True,rs=42) # Initializes self.dataSample and self.filesSample
        # Storing all reductions related to this object's data in its own reductions dictionary.
        self.reductions = dict()
    
    def sample(self, numLCs=10000, df='self',tabby=True, replace=True,rs=False):
        """
        Args:
            numLCs (int) - size of sample
            df ('self' or pd.DataFrame) - if self, will sample the object, if given a dataframe, will sample the given dataframe
            tabby (boolean) - if true, will ensure Boyajian's star is part of the sample.
            replace (boolean) - if true will replace the existing self.dataSample (used primarily for visualization)
            rs (boolean or int) - if int, will provide a set random state for the sample.
        Returns:
            sample (pd.DataFrame) - a randomly sampled subset from the givne dataframe
        """
        if type(df)==str:
            df = self.data
        
        assert (numLCs < len(df.index)),"Number of samples greater than the number of files."
        
        if type(rs)==int:
            sample = df.sample(n=numLCs,random_state=rs)
        else:
            sample = df.sample(n=numLCs)
        
        if tabby:
            if not sample.index.str.contains('8462852').any():
                sample = sample.drop(sample.index[0])
                sample = sample.append(self.data[self.data.index.str.contains('8462852')])
            
        if replace:
            self.dataSample = sample
            self.filesSample = sample.index          
        return sample
    
    def km_out(self,df='self',k=1):
        if type(df)==str:
            df = self.dataSample.iloc[:,0:60]
        labels = km_outliers.kmeans_w_outliers(df,k)
        self.KM_labels = labels
        return labels
    
    def db_out(self,df='self',neighbors=4,check_tabby=False,verbose=True):
        if type(df)==str:
            df = self.dataSample.iloc[:,0:60]
        labels = db_outliers.dbscan_w_outliers(data=df,min_n=neighbors,check_tabby=check_tabby,verbose=verbose)
        self.DB_labels = labels
        return labels
    
    def pca_red(self,df='self',red_name='PCA90',var_rat=0.9,scaled=False,verbose=True):
        if type(df)!=pd.core.frame.DataFrame:
            df = self.data
        pcaRed = qt.pca_red(df,var_rat,scaled,verbose)
        self.reductions[red_name]=pcaRed
        return pcaRed
    
    def save(self,of=False):
        """
        Pickles the whole object
        Args:
            of (str) - Defaults to out.coo. File to save object to. Using .coo to demarcate Cluster-Outlier Objects. 
        """
        
        if type(of)!=str:
            of=self.output_file
        else:
            self.output_file=of
        try:
            with open(of,'wb') as file:
                pickle.dump(self,file)  
        except:
            print("Something went wrong, check output file path.")

    def plot(self,df=False,pathtofits=False,clusterLabels='db_out',reduced=False):
        if type(df)!=pd.core.frame.DataFrame:
            df = self.data
        if reduced:
            data = df.iloc[:,[0,1]]
        else:
            dataReduced = qt.pca_red(df)
            data = dataReduced.iloc[:,[0,1]]
            
        if type(pathtofits)!=str:
            pathtofits = self.fitsDir
            
        qt.interactive_plot(data,pathtofits,clusterLabels)