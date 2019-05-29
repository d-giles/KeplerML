import numpy as np
from multiprocessing import Pool,cpu_count
import pandas as pd
import pyfits
import sklearn
from sklearn import preprocessing
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN, KMeans

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn.decomposition import PCA

from accelerate.cuda import cuda

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

import keplerml
import km_outliers
import db_outliers

class clusterOutliers(object):
    def __init__(self,feats,fitsDir):
        self.feats = feats
        self.data = pd.read_csv(feats,index_col=0)
        if fitsDir[-1]=="/":
            self.fitsDir = fitsDir
        else:
            self.fitsDir = fitsDir+"/"
        self.files = self.data.index
        # Initializing the data and files samples with the first 100 entries.
        self.dataSample = self.data.iloc[:100]
        self.filesSample = self.files[:100]
        self.importedForPlotting = False
        self.sampleGenerated = False
        self.sampleTSNE = False
        
    def read_kepler_curve(self,file):
        """
        Given the path of a fits file, this will extract the light curve and normalize it.
        """
        lc = pyfits.getdata(file)
        t = lc.field('TIME')
        f = lc.field('PDCSAP_FLUX')
        err = lc.field('PDCSAP_FLUX_ERR')

        err = err[np.isfinite(t)]
        f = f[np.isfinite(t)]
        t = t[np.isfinite(t)]
        err = err[np.isfinite(f)]
        t = t[np.isfinite(f)]
        f = f[np.isfinite(f)]
        err = err/np.median(f)
        nf = f / np.median(f)

        return t, nf, err
    
    def sample(self, numLCs=10000, df='self',tabby=True):
        """
        Rerunning this, or randSample will replace the previous random sample.
        """
        if type(df)==str:
            df = self.data
            
        assert (numLCs < len(df.index)),"Number of samples greater than the number of files."
        print("Creating random file list...")
        sample = df.sample(n=numLCs)
        
        if tabby:
            print("Checking for Tabby...")
            if not sample.index.str.contains('8462852').any():
                print("Adding Tabby...")
                sample = sample.drop(sample.index[0])
                sample = sample.append(self.data[self.data.index.str.contains('8462852')])
        
        self.dataSample = sample
        self.filesSample = sample.index
        self.sampleGenerated = True
        return sample
    
    def scale_data(self,data):
        scaler = preprocessing.StandardScaler().fit(data)
        scaledData = scaler.transform(data)
        return scaledData
    
    def tsne_fit(self,data,perplexity='auto',scaled=False):
        """
        Performs a t-SNE dimensionality reduction on the data sample generated.
        Uses a PCA initialization and the perplexity given, or defaults to 1/10th the amount of data.
        
        Appends the dataSample dataframe with the t-SNE X and Y coordinates
        Returns tsneX and tsneY
        """
        if type(perplexity)==str:
            perplexity=len(data)/10
        if not scaled:
            scaler = preprocessing.StandardScaler().fit(data)
            scaledData = scaler.transform(data)
        tsne = TSNE(n_components=2,perplexity=perplexity,init='pca',verbose=True)
        fit=tsne.fit_transform(scaledData)
        tsne_x = fit.T[0]
        tsne_y= fit.T[1]
        # Goal is to minimize the KL-Divergence
        if sklearn.__version__ == '0.18.1':
            print("KL-Divergence was %s"%tsne.kl_divergence_ )
        print("Done.")
        return tsne_x,tsne_y
    
    def pca_fit(self,df='self'):
        """
        Performs a 2D PCA on the given dataframe
        """
        if type(df)==str:
            df = self.dataSample.iloc[:,0:60]
        scaler = preprocessing.StandardScaler().fit(df)
        scaledData = scaler.transform(df)
        pca = PCA(n_components=2)
        pca_fit = pca.fit_transform(scaledData)
        return pca_fit
    
    def km_out(self,df='self',k=1):
        if type(df)==str:
            df = self.dataSample.iloc[:,0:60]
        labels = km_outliers.kmeans_w_outliers(df,k)
        return labels
    
    def db_out(self,df='self',neighbors=4,check_tabby=False,verbose=True):
        if type(df)==str:
            df = self.dataSample.iloc[:,0:60]
        labels = db_outliers.dbscan_w_outliers(data=df,min_n=neighbors,check_tabby=check_tabby,verbose=verbose)
        return labels
    
    def save(self,df='self',of='self'):
        """
        By default saves the modified data to the original file.
        """
        if type(df)==str:
            of=self.feats
        if type(df)==str:
            df=self.data
        df.to_csv(of)
        
    def plot(self,df='self',pathtofits='self',
                    clusterLabels='db_out',reduction_method="tsne"):
        
        if type(df) == str:
            if df == 'self':
                df = self.dataSample
        if type(pathtofits)==str:
            pathtofits = self.fitsDir
        
        root = Tk.Tk()
        root.wm_title("Scatter")
        """--- import light curve data ---"""
        files = df.index
        
        if type(clusterLabels) == str:
            labels = df[clusterLabels]
        else:
            labels = clusterLabels

        cNorm  = colors.Normalize(vmin=0, vmax=max(labels))
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='jet')
        
        data_out = df[labels==-1]
        files_out = data_out.index
        data_cluster = df[labels!=-1]
        files_cluster = data_cluster.index

        if reduction_method=='tsne':
            data = np.array(df[['tsne_x','tsne_y']])
            
            outX = data_out.tsne_x
            outY = data_out.tsne_y

            clusterX = data_cluster.tsne_x
            clusterY = data_cluster.tsne_y
            
        elif reduction_method=='pca':
            data = np.array(df[['pca_x','pca_y']])

            outX = data_out.pca_x
            outY = data_out.pca_y
            
            clusterX = data_cluster.pca_x
            clusterY = data_cluster.pca_y

        """--- Organizing data and Labels ---"""


        if df.index.str.contains('8462852').any():
            tabbyInd = list(df.index).index(df[df.index.str.contains('8462852')].index[0])            
        else:
            tabbyInd = 0
            
        fig = Figure(figsize=(20,10))

        # a tk.DrawingArea
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        # Toolbar to help navigate the data (pan, zoom, save image, etc.)
        toolbar = NavigationToolbar2TkAgg(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        gs = gridspec.GridSpec(2,100)

        with sns.axes_style("white"):
            # empty subplot for scattered data
            ax = fig.add_subplot(gs[0,:62])
            # empty subplot for lightcurves
            ax2 = fig.add_subplot(gs[1,:])
            # empty subplot for center detail
            ax3 = fig.add_subplot(gs[0,67:])


        @cuda.jit
        def distance_cuda(dx,dy,dd):
            bx = cuda.blockIdx.x # which block in the grid?
            bw = cuda.blockDim.x # what is the size of a block?
            tx = cuda.threadIdx.x # unique thread ID within a blcok
            i = tx + bx * bw
            if i>len(dd):
                return
            
            # d_ph is a distance placeholder
            d_ph = (dx[i]-dx[0])**2+(dy[i]-dy[0])**2

            dd[i]=d_ph**.5
            return
        
        def distances(pts,ex,ey):
            # Calculates distances between N points
            pts = np.array(pts)
            N=len(pts)

            # Allocate host memory arrays
            # Transpose pts array to n_dims x n_pts, each index of x contains all of a dimensions coordinates
            XT = np.transpose(pts)
            x = np.array(XT[0])
            x = np.insert(x,0,ex)
            y = np.array(XT[1])
            y = np.insert(y,0,ey)
            d = np.zeros(N)


            # Allocate and copy GPU/device memory
            d_x = cuda.to_device(x)
            d_y = cuda.to_device(y)
            d_d = cuda.to_device(d)

            threads_per_block = 128
            number_of_blocks =N/128+1 

            distance_cuda [ number_of_blocks, threads_per_block ] (d_x,d_y,d_d)

            d_d.copy_to_host(d)
            return d[1:]   
        
        def calcClosestDatapoint(X, event):
            """Calculate which data point is closest to the mouse position.

            Args:
                X (np.array) - array of points, of shape (numPoints, 2)
                event (MouseEvent) - mouse event (containing mouse position)
            Returns:
                smallestIndex (int) - the index (into the array of points X) of the element closest to the mouse position
            """
            ex,ey = ax.transData.inverted().transform((event.x,event.y))
            
            #distances = [distance (XT[:,i], event) for i in range(XT.shape[1])]
            Ds = distances(X,ex,ey)
            return np.argmin(Ds)

        def drawData(index):
            # Plots the lightcurve of the point chosen
            ax2.cla()
            f = pathtofits+df.index[index]
            t,nf,err=self.read_kepler_curve(f)

            x=t
            y=nf

            axrange=0.55*(max(y)-min(y))
            mid=(max(y)+min(y))/2
            yaxmin = mid-axrange
            yaxmax = mid+axrange
            if yaxmin < .95:
                if yaxmax > 1.05:
                    ax2.set_ylim(yaxmin,yaxmax)
                else:
                    ax2.set_ylim(yaxmin,1.05)
            elif yaxmax > 1.05:
                ax2.set_ylim(.95,yaxmax)
            else:
                ax2.set_ylim(.95,1.05)

            if files[index] in files_cluster:
                color = 'blue'
            else:
                color = 'red'
            ax2.plot(x, y, 'o',markeredgecolor='none', c=color, alpha=0.2)
            ax2.plot(x, y, '-',markeredgecolor='none', c=color, alpha=0.7)
            #ax2.set_title(files[index][:13],fontsize = 20)
            ax2.set_xlabel('Time (Days)',fontsize=22)
            ax2.set_ylabel(r'$\frac{\Delta F}{F}$',fontsize=30)

            fig.suptitle(files[index][:13],fontsize=30)

            canvas.draw()

        def annotatePt(XT, index):
            """Create popover label in 3d chart

            Args:
                X (np.array) - array of points, of shape (numPoints, 3)
                index (int) - index (into points array X) of item which should be printed
            Returns:
                None
            """
            x2, y2 = XT[index][0], XT[index][1]
            # Either update the position, or create the annotation
            if hasattr(annotatePt, 'label'):
                annotatePt.label.remove()
                annotatePt.emph.remove()
            if hasattr(annotatePt, 'emphCD'):
                annotatePt.emphCD.remove()
            x_ax = ax.get_xlim()[1]-ax.get_xlim()[0]
            y_ax = ax.get_ylim()[1]-ax.get_ylim()[0]
            # Get data point from array of points X, at position index
            annotatePt.label = ax.annotate( "",
                xy = (x2, y2), xytext = (x2+.1*x_ax, y2+.2*y_ax),
                arrowprops = dict(headlength=20,headwidth=20,width=6,shrink=.1,color='black'))
            annotatePt.emph = ax.scatter(x2,y2,marker='o',s=50,c='orange')
            if files[index] in files_cluster:
                annotatePt.emphCD = ax3.scatter(x2,y2,marker='o',s=150,c='orange')
            else:
                annotatePt.emphCD = ax.scatter(x2,y2,marker='o',s=50,c='orange')
            canvas.draw()


        def onMouseClick(event, X):
            """Event that is triggered when mouse is clicked. Shows lightcurve for data point closest to mouse."""
            #XT = np.array(X.T) # array organized by feature, each in it's own array
            closestIndex = calcClosestDatapoint(X, event)
            drawData(closestIndex)

        def onMouseRelease(event, X):
            #XT = np.array(X.T)
            closestIndex = calcClosestDatapoint(X, event)
            annotatePt(X,closestIndex)
            #for centerIndex in centerIndices:
            #    annotateCenter(XT,centerIndex)

        def connect(X):
            if hasattr(connect,'cidpress'):
                fig.canvas.mpl_disconnect(connect.cidpress)
            if hasattr(connect,'cidrelease'):
                fig.canvas.mpl_disconnect(connect.cidrelease)

            connect.cidpress = fig.canvas.mpl_connect('button_press_event', lambda event: onMouseClick(event,X))
            connect.cidrelease = fig.canvas.mpl_connect('button_release_event', lambda event: onMouseRelease(event, X))

        def redraw():       
            # Clear the existing plots
            ax.cla()
            ax2.cla()
            ax3.cla()
            # Set those labels
            ax.set_xlabel("Reduced X",fontsize=18)
            ax.set_ylabel("Reduced Y",fontsize=18)
            # Scatter the data
            ax.scatter(outX, outY,c="red",s=30,alpha=.5)

            ax.hexbin(clusterX,clusterY,mincnt=5,bins="log",cmap="viridis",gridsize=35)
            hb = ax3.hexbin(clusterX,clusterY,mincnt=5,bins="log",cmap="viridis",gridsize=35)
            cb = fig.colorbar(hb)
            """
            ax.scatter(clusterX,clusterY,s=30,c='b')
            ax3.scatter(clusterX,clusterY)
            """
            ax3.set_title("Center Density Detail",fontsize=14)
            ax3.set_xlabel("Reduced X",fontsize=18)
            ax3.set_ylabel("Reduced Y",fontsize=18)
            ax.tick_params(axis='both',labelsize=18)
            ax2.tick_params(axis='both',labelsize=18)
            ax3.tick_params(axis='both',labelsize=18)
            #for centerIndex in centerIndices:
            #    annotateCenter(currentData1,centerIndex)

            if hasattr(redraw,'cidenter'):
                    fig.canvas.mpl_disconnect(redraw.cidenter)
                    fig.canvas.mpl_disconnect(redraw.cidexit)
            connect(data)

            annotatePt(data,tabbyInd)

            drawData(tabbyInd)
            #fig.savefig('Plots/Q16_PCA_kmeans/Tabby.png')
            canvas.draw()
            canvas.show()
        print("Plotting.")
        redraw()            
            
        def _delete_window():
            print("window closed.")
            root.destroy()
            sys.exit()
    
        root.protocol("WM_DELETE_WINDOW",_delete_window)
        root.mainloop()