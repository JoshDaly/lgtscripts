#!/usr/bin/env python
###############################################################################
#
# __make_scatter_plot.py - make a scatter plot. 
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2014"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = ""
__status__ = "Development"

###############################################################################

import argparse
import sys
import operator
import random

from multiprocessing import Pool
from subprocess import Popen, PIPE

#import os
#import errno

import numpy as np
np.seterr(all='raise')

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def mungeLengthValue(val, mode="log"):
    if mode == "log":
        return np.log(val)

def doWork( args ):
    """ Main wrapper"""
    
    """make scatter plot from a tab delimited file 
    img_id_a        contig_a        img_id_b        contig_b        hits    length
    """
    # variables
    ids_list = [] # list to store uniquq ids
    unordered_ids_list = [] # contained unordered img Ids 
    sorted_ids_list = [] # contains IMG IDs ordered according to sorted_master
    #master_ids_list = [] # 
    ids_dict = {} # dict to store {id:[[id,hits,length],...]}
    tax_dict = {} # dict to store {id:tax_string}
    master_tax = [] #list to store img_ids taxa information
    sorted_master_tax = [] # Master sorted img IDs by tax string
    img_sorted_dict = {} # contains img_id -> x coordinate
    working_ids_list = [] # sorted img_ids 
    plot_border = args.plot_border 
    ordered_tax_string_lowest = [] # store lowest taxonomical information in order
    ordered_tax_string_phylum = [] # store phylum-level information in order
    subs_ordered_tax_string_phylum = [] #subset of unique entries in ordered_tax_string_lowest
    ticks = [] # array of numbers, for producing ticks on figure
    unique_taxa = {} #dictionary to check for unique taxa for plotting rectangles
    colours = [] #list of colours produced from input colour file
    Lengths = []
    Hits = []
    
    
    #change arguments into variables
    output_file = args.output_file + ".%s" % (args.image_format)
    image_format = args.image_format
    
        
    #-----
    """
    1/ alphabetically sort img ids by tax string
    2/ create a master list of img IDs
    A00000033       637000039       k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Brucellaceae; g__Brucella; s__abortus;
    """
        
    with open(args.master_tax_file, "r") as fh:
        for l in fh:
            img_id = l.split('\t')[1].rstrip()
            tax_string = l.split('\t')[2].rstrip()
            #master_ids_list.append(img_id)
            master_tax.append(tax_string)
            tax_dict[img_id] = tax_string
        sorted_master_tax = sorted(tax_dict.iteritems(), key=operator.itemgetter(1))
    # make a dictionary containing img ID and order
    for i, v in enumerate(sorted_master_tax):
        img_sorted_dict[v[0]] = int(i)
    
    #-----
    
    #read through file line by line
    with open(args.input_file, "r") as fh:
        # capture the file header
        header = fh.readline()
                       
        #read through the file line by line
        for l in fh:
            # assign each column in line to a variable
            id_a = l.split('\t')[0]
            id_b = l.split('\t')[2]
            hits = l.split('\t')[4].rstrip()
            length = l.split('\t')[5].rstrip()
            Hits.append(int(hits))
            Lengths.append(int(length))
            
            # check to see if key is in hash
            try:
                ids_dict[id_a][id_b] = [int(hits),int(length)]
            except KeyError:
                ids_dict[id_a] = {id_b:[int(hits),int(length)]}
                #unordered_ids_list.append(img_sorted_dict[id_a])
                
            try:                
                ids_dict[id_b][id_a] = [int(hits),int(length)]  
            except KeyError: 
                ids_dict[id_b] = {id_a:[int(hits),int(length)]}
                #unordered_ids_list.append(img_sorted_dict[id_b])  
                
        # make ids_list a numpy array
        ids_list= np.array(ids_dict.keys())
        for key in ids_dict.keys():
            unordered_ids_list.append(img_sorted_dict[key])
        # ordered img_ids  
        working_ids_list = ids_list[np.argsort(unordered_ids_list)]
    
    #-----
        
    #print out the tax string
    for i in working_ids_list:
        g = tax_dict[i].split(';')[5] #genus
        f = tax_dict[i].split(';')[4] #family
        o = tax_dict[i].split(';')[3] #order
        c = tax_dict[i].split(';')[2] #class
        p = tax_dict[i].split(';')[1] #phylum
        
        ordered_tax_string_phylum.append(p)
        
        if len(g) == 4:
            if len(f) == 4:
                if len(o) == 4:
                    if len(c) == 4:
                        if len(p) == 4:
                            print "error"
                        else:
                            #print p
                            ordered_tax_string_lowest.append(p)
                    else:
                        #print c
                        ordered_tax_string_lowest.append(c)
                else:
                    #print o
                    ordered_tax_string_lowest.append(o) 
            else:
                #print f
                ordered_tax_string_lowest.append(f)   
        else:
            #print g
            ordered_tax_string_lowest.append(g)
    
    #print ordered_tax_string_lowest
    #-----
    #randomise list
    #random.shuffle(ids_list)
    
    #-----
    """ Create scatter plot data 
    """
    
    #create scatter plot variable, x and y
    xs = [] 
    ys = [] 
    val = []
    heatmap_colours = [] # list of colours to be used for the heatmap (scatterplot)
    #BuPu Hex colours
    "#edf8fb" "#bfd3e6" "#9ebcda" "#8c96c6" "#8c6bb1" "#88419d" "#6e016b"
    #OrRd Hex colours
    "#fef0d9"    "#fdd49e"    "#fdbb84"    "#fc8d59"  "#ef6548"  "#d7301f"    "#990000"
    
    # variables for heatmap
    hits_max = int(max(Hits))
    length_max = int(max(Lengths))
    
    
    # loop over dictionary, to produce x and y
    for x in range(len(working_ids_list)): 
        try:
            fred = ids_dict[working_ids_list[x]]
            for y in range(x+1, len(working_ids_list)):
                try:
                    # Upper triangle - Hits data
                    val.append(fred[working_ids_list[y]][0])
                    xs.append(x)
                    ys.append(y)
                    if fred[working_ids_list[y]][0] <= (hits_max * 0.1):
                        heatmap_colours.append("#edf8fb")
                    if fred[working_ids_list[y]][0] > (hits_max * 0.1) and fred[working_ids_list[y]][0] <= (hits_max * 0.15):
                        heatmap_colours.append("#bfd3e6")
                    if fred[working_ids_list[y]][0] > (hits_max * 0.15) and fred[working_ids_list[y]][0] <= (hits_max * 0.2):
                        heatmap_colours.append("#9ebcda")
                    if fred[working_ids_list[y]][0] > (hits_max * 0.2) and fred[working_ids_list[y]][0] <= (hits_max * 0.3):
                        heatmap_colours.append("#8c96c6")
                    if fred[working_ids_list[y]][0] > (hits_max * 0.3) and fred[working_ids_list[y]][0] <= (hits_max * 0.4):
                        heatmap_colours.append("#8c6bb1")
                    if fred[working_ids_list[y]][0] > (hits_max * 0.4) and fred[working_ids_list[y]][0] <= (hits_max * 0.5):
                        heatmap_colours.append("#88419d")
                    if fred[working_ids_list[y]][0] > (hits_max * 0.5):
                        heatmap_colours.append("#6e016b")

                    # Lower triangle - Length data
                    val.append(fred[working_ids_list[y]][1])
                    xs.append(y)
                    ys.append(x)
                    if fred[working_ids_list[y]][1] <= (length_max * 0.01):
                        heatmap_colours.append("#fef0d9")
                    if fred[working_ids_list[y]][1] > (length_max * 0.01) and fred[working_ids_list[y]][1] <= (length_max * 0.1):
                        heatmap_colours.append("#fdd49e")
                    if fred[working_ids_list[y]][1] > (length_max * 0.1) and fred[working_ids_list[y]][1] <= (length_max * 0.15):
                        heatmap_colours.append("#fdbb84")
                    if fred[working_ids_list[y]][1] > (length_max * 0.15) and fred[working_ids_list[y]][1] <= (length_max * 0.2):
                        heatmap_colours.append("#fc8d59")
                    if fred[working_ids_list[y]][1] > (length_max * 0.2) and fred[working_ids_list[y]][1] <= (length_max * 0.3):
                        heatmap_colours.append("#ef6548")
                    if fred[working_ids_list[y]][1] > (length_max * 0.3) and fred[working_ids_list[y]][1] <= (length_max * 0.5):
                        heatmap_colours.append("#d7301f")
                    if fred[working_ids_list[y]][1] > (length_max * 0.5):
                        heatmap_colours.append("#990000")
                        
                except KeyError:
                    pass
        except KeyError:
            pass
    
       
    #-----
    val = np.sqrt(val)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # paint a rectangle over plot
    currentAxis = plt.gca()
    
    #-----
    """ Taxonomy(phylum-level) grids """
    # Create dictionary with unique taxa strings, and ordered subset of unique taxa strings
    for i,v in enumerate(ordered_tax_string_phylum):
        if v in unique_taxa:
            unique_taxa[v][1] += 1
        else:
            unique_taxa[v] = [i,1]
            subs_ordered_tax_string_phylum.append(v)    
    #-----
    # produce a random list of colours of len(working_ids_list)
    with open(args.colour_file,"r") as fh:
        for l in fh:
            colours.append(str(l.rstrip()))
    #randomise list of colours
    random.shuffle(colours)
    #print colours

    # produce rectangles for each unique taxa string
    for i in subs_ordered_tax_string_phylum:
        #y-axis rectangle
        currentAxis.add_patch(Rectangle((0-plot_border, unique_taxa[i][0]), len(working_ids_list)+plot_border-1, unique_taxa[i][1]-1, edgecolor = "grey",facecolor="none",alpha=0.1))
        #currentAxis.add_patch(Rectangle((0-plot_border, unique_taxa[i][0]), len(working_ids_list)+plot_border, unique_taxa[i][1], edgecolor = colours[v],facecolor="none" ,alpha=0.45))
        #x-axis rectangle
        currentAxis.add_patch(Rectangle((unique_taxa[i][0],0-plot_border), unique_taxa[i][1]-1, len(working_ids_list)+plot_border-1,edgecolor = "grey",facecolor="none",alpha=0.1))
        #currentAxis.add_patch(Rectangle((unique_taxa[i][0],0-plot_border), unique_taxa[i][1], len(working_ids_list)+plot_border,edgecolor = colours[v],facecolor="none",alpha=0.45))    
        
    #-----
    """ Produce scatter plot """
    ax.scatter(xs,
               ys,
               s=3,
               marker='s',
               alpha=1,
               c=heatmap_colours,
               edgecolors = 'grey',
               linewidths = 0.1,
               )
    
    
    #Plot labels
    plt.xlabel('Cumulative length of contigs (bp)')
    plt.ylabel('No. of hits')
    plt.title("Gut-Oral LGT events")
    #ax.grid(b=True,which='both')
    
    #set the number of tick labels
    for i in range(len(working_ids_list)):
        ticks.append(i)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    
    #Change label names according to taxonomy
    labels = [item.get_text() for item in ax.get_xticklabels()]
    for i,v in enumerate(labels):
        labels[i] = ordered_tax_string_lowest[i]
    ax.set_xticklabels(labels,size=args.text_size)
    ax.set_yticklabels(labels,size=args.text_size)
    
    #rotate x labels 90deg
    plt.xticks(rotation=90)
    
    #adjust margin size
    #plt.subplots_adjust(bottom=0.2)
    plt.tight_layout()
    
    #plt.xticks(np.arange(min(xs), max(xs)+1,1.0))
    #plt.yticks(np.arange(min(xs), max(xs)+1,1.0))
    # set the plot axes
    plt.axis([plot_border*-1,
              len(working_ids_list)+plot_border,
              plot_border*-1,
              len(working_ids_list)+plot_border])

    if args.show_plot == "True":
        plt.show()
    else:
        plt.savefig(output_file,dpi=args.dpi,format=str(image_format))                   
       
    """
# run somethign external in threads
pool = Pool(6)
cmds = ['ls -l', 'ls -alh', 'ps -ef']
print pool.map(runCommand, cmds)
"""

    """
# parse a file
try:
with open(filename, "r") as fh:
for line in fh:
print line
except:
print "Error opening file:", filename, exc_info()[0]
raise
"""

    """
fig = plt.figure()

#-----
# make a 3d plot
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0],
points[:,1],
points[:,2],
#edgecolors='none',
#c=colors,
#s=2,
#marker='.'
)

#-----
# make a 2d plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(points[:,0],
points[:,1],
'*g')

#-----
# show figure
plt.show()
# or save figure
plt.savefig(filename,dpi=300,format='png')

#-----
# clean up!
plt.close(fig)
del fig
"""

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_file', help="Tab delimited file containing IDs, hits and length")
    parser.add_argument('-o','--output_file', default = "scatter_plot", help="Output file")
    parser.add_argument('-f','--image_format', default='png', help="Image format; default: png")
    parser.add_argument('-d','--dpi', type=int, default=400, help="Set image dpi")
    parser.add_argument('-t','--master_tax_file', help="Master taxa file containing all img IDs")
    parser.add_argument('-pt','--plot_title', help="Plot title")
    parser.add_argument('-pb','--plot_border', type=int, default=1,help="Adjust plot border to fit data")
    parser.add_argument('-ts','--text_size', type=int, default=2,help="Set text size for x and y labels")
    parser.add_argument('-sp','--show_plot', default='False',help="Show plot in terminal options: True. default: False")
    parser.add_argument('-cf','--colour_file', help="Colours file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
