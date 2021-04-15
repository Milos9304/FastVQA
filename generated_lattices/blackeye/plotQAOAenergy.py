#andromeda

import matplotlib.pyplot as plt
import numpy as np
import sys
from glob import glob

from bokeh.document import Document
from bokeh.embed import file_html
from bokeh.layouts import gridplot
from bokeh.models import (Title, BasicTicker, Circle, ColumnDataSource, DataRange1d,
                          Grid, LinearAxis, PanTool, Plot, WheelZoomTool,)
from bokeh.resources import INLINE
from bokeh.sampledata.iris import flowers
from bokeh.util.browser import view

xdr = DataRange1d(bounds=None)
ydr = DataRange1d(bounds=None)

#for filename in glob("statsfiles/*"):

def make_plot(dim, inst, xax=False, yax=False):
    
    xdr = DataRange1d(bounds=None)
    ydr = DataRange1d(bounds=None)

    mbl = 40 if yax else 0
    mbb = 40 if xax else 0

    plot = Plot(
        x_range=xdr, y_range=ydr, background_fill_color="#efe8e2",
        border_fill_color='white', plot_width=400 + mbl, plot_height=400 + mbb,
        min_border_left=2+mbl, min_border_right=2, min_border_top=2, min_border_bottom=2+mbb)
    #circle = Circle(x=xname, y=yname, fill_color="color", fill_alpha=0.2, size=4, line_color="color")
    circle = Circle(x="iters", y="energy", fill_alpha=0.2, size=4 )
    
    with open("statsfiles/statsfile_q_"+str(dim)+"_1_20_b_"+str(inst)+".txt", encoding='utf8') as f:

        skip = 0
        energies = []
        for line in f:
            if skip % 10 == 0:
                energies.append(float(line))
                skip = 0
            skip += 1

        data = {"iters": np.linspace(1,len(energies),len(energies)),
                "energy": energies}

    source = ColumnDataSource(data=data)

    r = plot.add_glyph(source, circle)

    xdr.renderers.append(r)
    ydr.renderers.append(r)

    xticker = BasicTicker()
    if xax:
        xaxis = LinearAxis()
        xaxis.axis_label = str(inst)
        plot.add_layout(xaxis, 'below')
        xticker = xaxis.ticker
    plot.add_layout(Grid(dimension=0, ticker=xticker))

    yticker = BasicTicker()
    if yax:
        yaxis = LinearAxis()
        yaxis.axis_label = str(dim)
        yaxis.major_label_orientation = 'vertical'
        plot.add_layout(yaxis, 'left')
        yticker = yaxis.ticker
    plot.add_layout(Grid(dimension=1, ticker=yticker))

    plot.add_tools(PanTool(), WheelZoomTool())

    return plot

def write_html_plot(dims, insts):

    plots = []
    doc = Document()

    for dim in dims:
        row = []
        for inst in insts:
            xax = (dim == dims[-1])
            yax = (inst == insts[0])
            plot = make_plot(dim, inst, xax, yax)
            row.append(plot)
        plots.append(row)

    grid = gridplot(plots)

    doc.add_root(grid)

    doc.validate()
    filename = "plots/plot_"+str(dims[0])+"-"+str(dims[-1])+"_inst"+str(insts[0])+"-"+str(insts[-1])+".html"
    with open(filename, "wb") as f:
        f.write(file_html(doc, INLINE, "Iris Data SPLOM").encode("utf-8"))
    print("Wrote %s" % filename)

    del doc
    del grid
    #view(filename)

if __name__ == "__main__":

    dims = [5,6,7]
    insts = [0,1,2,3,4]
    write_html_plot(dims, insts)

    dims = [8,9,10]
    insts = [0,1,2,3,4]
    write_html_plot(dims, insts)

"""

    fig, ax = plt.subplots(6, 10, sharex="col", sharey="row", figsize=(8, 8))
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.1, wspace=0.1)
    
    with open(filename, encoding='utf8') as f:
        
        energies = []

        for line in f:
            energy = float(line)
            energies.append(energy)
        
        filename_parts = filename.split('_')
        dim = filename_parts[2]
        inst = filename_parts[-1][:-4]
        
        ax[int(inst), int(dim)-5].scatter(np.linspace(1,len(energies),len(energies)), energies, c='r', s=1,alpha=0.6)
        break

# remove tick labels
for axi in ax.flat:
    for axis in [axi.xaxis, axi.yaxis]:
        axis.set_major_formatter(plt.NullFormatter())

plugins.connect(fig, plugins.LinkedBrush(ax[int(inst), int(dim)-5]))
mpld3.show()"""
