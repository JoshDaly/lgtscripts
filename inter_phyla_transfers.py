#!/usr/bin/env python
###############################################################################
#
# __inter_phyla_transfers__.py - Calculate the number of inter-phyla transfers. 
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
    #Objects
    intra_phyla= 0 
    inter_phyla= 0
    hits_accumulative= 0
    phylum_dict = {}
    
    #firmicutes = {}
    #actinobacteria = {}
    #synergistes = {}
    #epsilon_proteobacteria = {} 
    #fusobacteria = {}
    #bacteroidetes = {}
    #gamma_proteobacteria = {} 

    with open(args.phyla,"r") as fh:
        for l in fh:
            tabs= l.split()
            id= tabs[0].rstrip()
            phylum= tabs[1].rstrip()
            if phylum=="firmicutes":
                #firmicutes[id] = 0
                phylum_dict[id] = 0
            elif phylum=="actinobacteria":
                #actinobacteria[id] = 1
                phylum_dict[id] = 1
            elif phylum=="synergistes":
                #synergistes[id] = 2
                phylum_dict[id] = 2
            elif phylum=="proteobacteria":
                #epsilon_proteobacteria[id] = 3
                phylum_dict[id] = 3
            #elif phylum=="epsilon_proteobacteria":
                #epsilon_proteobacteria[id] = 3
                #phylum_dict[id] = 3
            elif phylum=="fusobacteria":
                #fusobacteria[id] = 4
                phylum_dict[id] = 4
            elif phylum=="bacteroidetes":
                #bacteroidetes[id] = 5
                phylum_dict[id] = 5
            #elif phylum=="gamma_proteobacteria":
                #gamma_proteobacteria[id] = 6
                #phylum_dict[id] = 6
    
    """img_id_a        genome_tree_id_a        contig_a        contig_length_a start_a stop_a  length_a        img_id_b        genome_tree_id_b        contig_b        contig_length_b start_b stop_b  length_b"""
    # read in transfers file
    with open(args.transfers_file,"r") as fh:
        # capture header
        header = fh.readline()
        for l in fh:
            tabs= l.split("\t")
            genome_tree_a= tabs[1]
            genome_tree_b=  tabs[8]
            hits_accumulative+= 1
            try:
                if phylum_dict[genome_tree_a] == phylum_dict[genome_tree_b]:
                    intra_phyla+= 1
                else:
                    inter_phyla+= 1
            except KeyError:
                pass
                print genome_tree_a
                print genome_tree_b
    print "Total:"+str(hits_accumulative)
    print "intra_phyla:"+str(intra_phyla)
    print "inter_phlya:"+str(inter_phyla)
    
    #for key in phylum_dict.keys():
    #    print key        
            
            
            
            
    
    

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
    parser.add_argument('-t','--transfers_file', help="File containing transfers information")
    parser.add_argument('-p','--phyla', help="File containing phyla information")
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
