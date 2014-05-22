#!/usr/bin/env python
###############################################################################
#
# __annotation_stats__.py - Provide information regarding annotation for each transfer
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

from multiprocessing import Pool
from subprocess import Popen, PIPE

import re

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


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

def doWork( args ):
    """ Main wrapper"""
    # global variables
    transfer_annotations = {}
    COG_groups= {}
    COG_categories= {}
    
    # Add unknown to Cog_categories
    COG_categories["unknown"]="unknown"
    
    #read in COG catergories
    """
    comma separated COG categories
    """
    with open(args.cog_cats,"r") as fh:
        # no header
        for l in fh: #line by line
            commas= l.split(",")
            cog_id= commas[0]
            cog_cat= commas[1]
            annotation= commas[2].rstrip()
            COG_categories[cog_id]=cog_cat
    
    
    #read in annotation file
    with open(args.annotation_file,"r") as fh:
        # no header 
        for l in fh:
            tabs = l.split("\t")
            lines = tabs[0].split("-")
            contig = lines[0].rstrip()
            id_a = lines[1].split(":")[1].rstrip()
            genome_tree_a = lines[2].split(":")[1].rstrip()
            start = int(lines[3].split(":")[1].rstrip())
            stop = int(lines[4].split(":")[1].rstrip())
            id_b = lines[5].split(":")[1].rstrip()
            genome_tree_b = lines[6].split(":")[1].rstrip()
            unique_id = lines[7].rstrip()
            COG= tabs[5]
            annotation = tabs[8]
            if not COG:
                COG = "unknown" 
            
            try:
                transfer_annotations[id_a][id_b][unique_id] =  [genome_tree_a,genome_tree_b,start,stop,annotation,COG]
            except KeyError:
                try:
                    transfer_annotations[id_a][id_b] = {unique_id:[genome_tree_a,genome_tree_b,start,stop,annotation,COG]}
                except KeyError:
                    transfer_annotations[id_a] = {id_b:{unique_id:[genome_tree_a,genome_tree_b,start,stop,annotation,COG]}} 
    
    for id_a in transfer_annotations.keys():
        for id_b in transfer_annotations[id_a]:
            for uid in transfer_annotations[id_a][id_b]: 
                COG= transfer_annotations[id_a][id_b][uid][-1] 
                if re.search("..*",COG_categories[COG]):
                    print COG
                
                #try:
                #    COG_groups[COG_categories[COG]]+= 1
                #except KeyError:
                #    COG_groups[COG_categories[COG]]= 1
                #try: 
                #    COG_groups[COG]+= 1
                #except KeyError:
                #    COG_groups[COG]=1
    #for COG in COG_groups:
    #    print COG+"\t"+str(COG_groups[COG])
                
                
                
    """
    try:
        transfer_annotations[id_a][id_b][contig] +=  [[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation,COG]]
    except KeyError:
        try:
            transfer_annotations[id_a][id_b][contig] = [[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation,COG]]
        except KeyError:
            try: 
                transfer_annotations[id_a][id_b] = {contig:[[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation,COG]]}
            except KeyError:
                transfer_annotations[id_a] = {id_b:{contig:[[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation,COG]]}}
    """ 
  

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
    parser.add_argument('-a','--annotation_file', help="Provide file containing transfer annotations")
    parser.add_argument('-c','--cog_cats', help="File containing COG categories")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
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
