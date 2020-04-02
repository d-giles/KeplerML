import pickle
import sys
sys.path.append('python')
import clusterOutliers
import numpy as np
from quarterTools import colors_for_plot
import pandas as pd
from lightkurve import search_lightcurvefile
from matplotlib.colors import to_hex

from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.layouts import row,gridplot,column
from bokeh.models import ColumnDataSource, Circle
from bokeh.models.widgets import Button, Select
from bokeh.events import Tap
from bokeh.palettes import Colorblind4,Category10
from bokeh.transform import factor_cmap, linear_cmap, log_cmap

##### Data Import #####

with open('data/output/Q8_sample.coo','rb') as file:
    Q8_coo = pickle.load(file)

data               = Q8_coo.data
Quarter            = 8
pca                = Q8_coo.reductions['PCA90']
x                  = pca.iloc[:,0]
y                  = pca.iloc[:,1]
files              = data.index
labels             = np.array(['KIC '+str(int(lc[4:13])) for lc in files])
dblabels           = Q8_coo.scores['DB_outliers']
dbcolors           = colors_for_plot(dblabels,'color_blind')
dblabels_as_string = np.array(dblabels,dtype=str)
scores             = Q8_coo.scores['k1_1000x10']
lin_colors         = colors_for_plot(-scores,'viridis')
log_colors         = colors_for_plot(-np.log(scores),'viridis')



##### Bokeh Setup #####
# Tools to use
TOOLS     = "pan,wheel_zoom,reset,tap,box_select,poly_select"
palette   = Colorblind4[0], Colorblind4[1], Colorblind4[3] # Color palette for outlier cmap
alphas    = {1:1,2:.9,3:.8,4:.75,5:.7,6:.65,7:.6,8:.55,9:.5,10:.5}
cmap_dict = {
    "Base"            : palette[0],
    "Outliers"        : [to_hex(c) for c in dbcolors],
    "Scores - Linear" : [to_hex(c) for c in lin_colors],
    "Scores - Log"    : [to_hex(c) for c in log_colors],
}

# Set up callbacks
def update_plot():
    try:
        inds = s1.selected.indices
        if len(inds)==0:
            inds=[0]
        if len(inds)>10:
            inds = inds[:10] # too many lightcurves would be illegible, 10 is pushing it.
    except:
        inds = [0]
    alpha = alphas[len(inds)]
    s2.data = {'t':[],'nf':[]}
    lcs = list(s1.data['desc'][inds])
    for i,ind in enumerate(inds):
        lc = lcs[i]

        # download Kepler lighcurve via Lightkurve
        lcf = search_lightcurvefile(lc, mission="Kepler", quarter=Quarter).download()
        # use the normalized PDCSAP flux 
        nlc = lcf.PDCSAP_FLUX.normalize()
        newdata = {'t':nlc.time,'nf':nlc.flux}
        s2.stream(newdata)
        lg.glyph.line_alpha=alpha
        #lg.glyph.line_color=palette[i]
        #s2.data = {'t':t,'nf':nf}

def update_colors(attrname, old, new):
    """ 
    Bokeh was being uncooperative with it's color mapper combined with the selection tool.
    I could get the colors to change when a different color selection was made,
    but it wouldn't update the color mapping to selected data. Whenever switching to a new
    color selection, selecting data would revert to the previous color mapping except where
    colors are given explicitly. So, I'm defining the colors explicitly using matplotlib
    instead of using the more convenient Bokeh mapper function which I can't find a way to save the 
    mapping as an array for the life of me. It's clunky and I hate it but it works.
    """
    if select.value == "Base":
        s1.data['colors'] = [cmap_dict[select.value]]*len(s1.data['x'])
    else:
        s1.data['colors'] = cmap_dict[select.value]

# Bokeh data
# create a column data source for the plots to share
s1 = ColumnDataSource(data={
    'x'     : x,
    'y'     : y,
    'desc'  : labels,
    'files' : files,
    'labels': dblabels_as_string,
    'scores': scores,
    'colors': [cmap_dict['Base']]*len(x) # See update_colors method for comments
}) 

s2=ColumnDataSource(data={
    't'  : [],
    'nf' : []
})

"""
# Idea for multi-line plotting with different colors:
ts = [[]]*10 # empty arrays for now 
nfs = [[]]*10
s2 = ColumnDataSource(data={
    't':ts,
    'nf':nfs,
    'colors':Category10[10]
})

plc = figure(plot_width=1000,plot_height=200,title=None)
lg = plc.multiline('t','nf',color='colors',source=s2)
def update_plot():
    try:
        inds = s1.selected.indices
        if len(inds)==0:
            inds=[0]
        if len(inds)>10:
            inds = inds[:10] # too many lightcurves would be illegible, 10 is pushing it.
    except:
        inds = [0]
    alpha = alphas[len(inds)]
    lcs = list(s1.data['desc'][inds])
    new_ts = ts
    new_nfs = nfs
    for i,ind in enumerate(inds):
        lc = lcs[i]

        # download Kepler lighcurve via Lightkurve
        lcf = search_lightcurvefile(lc, mission="Kepler", quarter=Quarter).download()
        # use the normalized PDCSAP flux 
        nlc = lcf.PDCSAP_FLUX.normalize()
        new_ts[i] = nlc.time
        new_nfs[i] = nlc.flux
    newdata = {'t':new_ts,'nf':new_nfs}
    s2.stream(newdata)
    lg.glyph.line_alpha=alpha
"""
# Set up widgets
button = Button(label='Plot Selected')
button.on_click(update_plot)

select = Select(title="Color Scheme:", value="Base", options=["Base", "Outliers", "Scores - Linear", "Scores - Log"])
select.on_change('value',update_colors)

# create a new plot and add a renderer
left = figure(tools=TOOLS, plot_width=1000, plot_height=600, title=None)
scatter = left.circle(x='x', y='y', fill_color='colors', line_color=None, size=4, source=s1)
left.on_event(Tap,update_plot)
# Planning to incorporate a detailed view of the cluster center on a right plot in the future
# create a plot for the light curve and add a line renderer
plc = figure(tools=TOOLS,plot_width=1000,plot_height=200,title=None)
lg = plc.line('t','nf',source=s2)
update_plot()

# Set up layouts and add to document
inputs = row(button,select)
#p = gridplot([[left,right]]) #future planning
layout = column(inputs, left, plc)
curdoc().add_root(layout)