import pandas as pd

from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.layouts import row,gridplot,column
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Button, Select
from bokeh.events import Tap
from bokeh.palettes import Category10

gui_df = pd.read_csv('gui_df.csv')
s1 = ColumnDataSource(gui_df)

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
        s1.data['colors'] = [Category10[10][0]]*len(s1.data['x'])
    else:
        s1.data['colors'] = s1.data[select.value]

def update_plot():
    try:
        ind = s1.selected.indices[0] # can only plot one at a time
    except:
        ind = 0

    lc_path = s1.data['files'][ind]
    lc_image = plc.image_url(url=[lc_path],x=0,y=0,w=10,h=3,anchor="bottom_left")

# Set up widgets
button = Button(label='Plot Selected')
button.on_click(update_plot)

select = Select(title="Color Scheme:", value="Base", options=["Base", "Outliers", "Scores - Linear", "Scores - Log"])
select.on_change('value',update_colors)

# create a new plot and add a renderer
TOOLS     = "pan,wheel_zoom,reset,tap,box_select,poly_select"
left = figure(tools=TOOLS, plot_width=1000, plot_height=600, title=None)

scatter = left.circle(x='x', y='y', fill_color='colors', line_color=None, size=4, source=s1)
left.on_event(Tap,update_plot)
plc = figure(plot_width=1000,plot_height=300,x_range=(0,10),y_range=(0,3))
plc.xaxis.axis_line_color = None
plc.yaxis.axis_line_color = None
plc.xaxis.major_tick_line_color = None
plc.xaxis.minor_tick_line_color = None
plc.yaxis.major_tick_line_color = None
plc.yaxis.minor_tick_line_color = None
plc.xaxis.major_label_text_color = None
plc.yaxis.major_label_text_color = None
lc_path = s1.data['files'][0]
lc_image = plc.image_url(url=[lc_path],x=0,y=0,w=10,h=3,anchor="bottom_left")

# Set up layouts and add to document
inputs = row(button,select)
#p = gridplot([[left,right]]) #future planning
layout = column(inputs, left, plc)
curdoc().add_root(layout)
    
