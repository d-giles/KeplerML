KeplerML
========
# Introduction
The work contained here has been developed for the purpose of outlier detection on light curves from the Kepler Data. 
There are methods to calculate features based on light curve variability, 
methods to identify outliers based on cluster membership, 
methods to score the relative outlying nature of each point within a dataset, 
and some user interface methods to explore the data. 
Clustering and scoring can be performed independently, though are intrinsically related to one another. 
Both require features to be calculated for the light curves as a prerequisite. 
The user interface tools have various prerequisites.

It is assumed that these methods will be applied to very large sets of data; each Kepler quarter contains roughly 160k light curves. Light curve features can be calculated for any number of light curves, and most methods for clustering and scoring will work with any number of files, but a minimum of 1k is recommended. The "Cluster Outlier Object" developed as the analysis container of this work cannot be created with less than 1k objects.

# Feature Calculation:
Recommended prerequisite: Filelist with lightcurve filenames.

Recommendation (example given for long cadence lightcurves from a Kepler quarter):
1. Collect lightcurves into a single directory.
2. Generate a filelist of the contents of that directory (the lightcurves)

    for f in \*llc.fits; do echo $f >> filelist; done
        
NOTE: For unknown reasons the code is having issues processing a whole quarter at once. Recommend splitting 
the files into at least 2 groups. An easy option is the splitting the files starting in 00 and 01
    
2. (alternate)

    for f in kplr00\*llc.fits; do echo $f >> Q??\_00filelist; done
    
    for f in kplr01\*llc.fits; do echo $f >> Q??\_01filelist; done
        
Where ?? is replaced by the quarter number.
        
Ways to use:
1. Run keplerml.py to calculate the lightcurve features, this will output a numpy array with the calculated features for each lightcurve.
    
        python keplerml.py path/to/filelist path/to/fitsfiles path/to/outputfile
2. Open Feature Calculation Example.ipynb in jupyter notebook to see examples of how to run feature calculation in a notebook.

Note: Using a 48-2.70GHz core linux computer (using 47 of the cores), processing 114,948 files took 54m:48s, which translates to 1.344 seconds to process a single file on one core. If you have less cores (most computers have 1-8 cores), multiply the number of files by the time to process a single file, and divide by the number of cores in the computer for an estimate on how long it will take to process.

# Clustering:
Prerequisite: Calculated features saved as a Pandas Dataframe or Numpy array.
Recommendation: Use the feature data in a Cluster Outlier Object.

See clustering example.

# Scoring:
Prerequisite: Calculated features saved as a Pandas Dataframe or Numpy array.
Recommendation: Use the feature data in a Cluster Outlier Object.

See scoring example.