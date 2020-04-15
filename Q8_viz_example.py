# Data prep for visualization
from bokeh.io import curdoc

import pickle
import sys
sys.path.append('python')
import clusterOutliers
import numpy as np
import pandas as pd

# Import full quarter coo
file_path    = 'data/output/Q8_sample.coo'
score_column = 'k1_1000x10'
reduction_name = 'PCA90'
Quarter = 8

with open(file_path,'rb') as file:
    coo = pickle.load(file)

ft_data = coo.data
# choose subset of data to plot, the full quarter will typically be too much.
# Recommendation: include all outlying points (determined however) and a sampling of the normal data
scores_df           = coo.scores.copy()
scores_df['KIC']    = [int(i[4:13]) for i in scores_df.index] # KIC as integers for easier parsing
scores_df['files']  = scores_df.index # keeping the files handy
scores_df.set_index('KIC',inplace=True) # setting the index to KIC integers
top_out_KIC         = scores_df.sort_values(score_column,ascending=False).head(1000).index # 
seed = 138633
np.random.seed(seed) # setting the random seed to the integer generated earlier
# 1000 pts sampled from the least outlying points (all points excluding the top outlying points)
rand_bottom_out_KIC = scores_df.sort_values(score_column,ascending=True).head(-1000).sample(1000).index
sample              = np.append(top_out_KIC,rand_bottom_out_KIC) # numpy array of KICs for our sample
samp_df             = scores_df.loc[sample,:]

# x and y coordinates 
# Using the PCA90 reduction (only one included in our example)
reduct_df    = coo.reductions[reduction_name].loc[samp_df.files,:]
samp_df['x'] = list(reduct_df.iloc[:,0])
samp_df['y'] = list(reduct_df.iloc[:,1])

# colors_for_plot maps rgba color values to the range of values in a numerical array, whether continuous or discrete
from quarterTools import colors_for_plot
rgbas = [colors_for_plot(samp_df['DB_outliers'],'color_blind'),
         colors_for_plot(-samp_df['k1_1000x10'],'viridis'), # flipping scores to match outlier colors
         colors_for_plot(-np.log(samp_df['k1_1000x10']),'viridis')]

# converted to hex values for bokeh
from matplotlib.colors import to_hex
cmaps  = ["Outliers","Scores - Linear","Scores - Log"] # choose sensible and descriptive names for the color maps

for i,rgba in enumerate(rgbas):
    samp_df[cmaps[i]] = [to_hex(c) for c in rgba]

samp_df["colors"] = samp_df[cmaps[0]] # setting default color values for plot to outliers

from bokeh.plotting import figure
from bokeh.layouts import row,gridplot,column
from bokeh.models import ColumnDataSource, Circle
from bokeh.models.widgets import Button, Select
from bokeh.events import Tap
from bokeh.palettes import Colorblind4,Category10
from bokeh.transform import factor_cmap, linear_cmap, log_cmap

from lightkurve import search_lightcurvefile

##### Data Import #####
s1 = ColumnDataSource(samp_df)

s2=ColumnDataSource(data={
    't'  : [],
    'nf' : []
})

##### Bokeh Setup #####
# Tools to use
TOOLS     = "pan,wheel_zoom,reset,tap,box_select,poly_select"
alphas    = {1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1} # for lc plot

# Set up callbacks
def auto_update_plot(attr,old,new):
    inds = new
    if len(inds)==0:
        inds=[0]
    elif len(inds)>10:
        print("Only showing 10 of {} selected.".format(len(inds)))
        inds = inds[:10] # too many lightcurves are illegible, 10 is pushing it.

    alpha = alphas[len(inds)]
    lc_files = list(s1.data['files'][inds])
    new_ts = [[]]*len(inds)
    new_nfs = [[]]*len(inds)
    reset_data = {'t':new_ts,'nf':new_nfs,'colors':np.array(Category10[10][:len(inds)])}
    s2.data=reset_data
    plc.title.text = "KIC "
    for i,ind in enumerate(inds):
        lc = lc_files[i][4:13]

        # download Kepler lighcurve via Lightkurve
        lcf = search_lightcurvefile(lc, mission="Kepler", quarter=Quarter).download()
        # use the normalized PDCSAP flux 
        nlc = lcf.PDCSAP_FLUX.normalize()
        new_ts[i] = nlc.time
        new_nfs[i] = nlc.flux+i
        plc.title.text += lc + " "
    newdata = {'t':new_ts,'nf':new_nfs,'colors':np.array(Category10[10][:len(inds)])}
    s2.data = newdata
    lg.glyph.line_alpha=alpha

# Idea for multi-line plotting with different colors:
ts = [[]] # empty arrays for now 
nfs = [[]]
s2 = ColumnDataSource(data={
    't':ts,
    'nf':nfs,
    'colors':[Category10[10][0]]
})

plc = figure(tools="pan,wheel_zoom,reset",plot_width=1000,plot_height=200)
plc.title.text = "Select a point to plot a light curve"
lg = plc.multi_line('t','nf',color='colors',source=s2)
# create a plot for the light curve and add a line renderer

def update_colors(attrname, old, new):
    s1.data['colors'] = s1.data[select.value]

# Bokeh data
# create a column data source for the plots to share

select = Select(title="Color Scheme:", value=cmaps[0], options=cmaps)
select.on_change('value',update_colors)

# create a new plot and add a renderer
left = figure(tools=TOOLS, plot_width=1000, plot_height=600, title=None)
scatter = left.circle(x='x', y='y', fill_color='colors', line_color=None, size=4, source=s1)
s1.selected.on_change('indices',auto_update_plot)
# Planning to incorporate a detailed view of the cluster center on a right plot in the future

# Set up layouts and add to document
inputs = row(select)
#p = gridplot([[left,right]]) #future planning
layout = column(inputs, left, plc)
curdoc().add_root(layout)
