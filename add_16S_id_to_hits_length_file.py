#!/usr/bin/env python
###############################################################################
#
# __add_16S_id_to_hits_length_file__.py - description!
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
    # objects
    id_percentage = {}
    
    """ read in 16S ID file"""
    with open(args.ID_file, "r") as fh:
        header = fh.readline()
        for l in fh:
            tabs = l.split("\t")
            id_a = tabs[0].rstrip()
            id_b = tabs[1].rstrip()
            identity = float(tabs[2].rstrip())
            try:
                id_percentage[id_a][id_b] = identity
            except KeyError:
                id_percentage[id_a] = {id_b:identity}
            try:
                id_percentage[id_b][id_a] = identity
            except KeyError:
                id_percentage[id_b] = {id_a:identity}
    
    """ read in hits_length file"""
    with open(args.hits_length, "r") as fh:
           # capture the file header
           header = fh.readline().rstrip()
           print header+"\t"+"ID%"               
           #read through the file line by line
           for l in fh:
               # assign each column in line to a variable
               tabs = l.split("\t")
               genome_tree_a = tabs[0]
               id_a = tabs[1]
               body_site_a = tabs[2]
               genome_tree_b = tabs[3]
               id_b = tabs[4]
               body_site_b = tabs[5]
               hits = tabs[6].rstrip()
               length = tabs[7].rstrip()
               #interaction = tabs[8].rstrip()
               #phage = tabs[9].rstrip()
               #plasmid = tabs[10].rstrip()
               #transposon = tabs[11].rstrip()
               id = 0
               try:
                   id = id_percentage[id_a][id_b]
               except KeyError:
                   print "not in percentage file"
               print "\t".join([genome_tree_a,
                                id_a,
                                body_site_a,
                                genome_tree_b,
                                id_b,
                                body_site_b,
                                hits,
                                length,
                                str(id)
                                ]) 
            
            
            
    
    

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
    parser.add_argument('-i','--hits_length', help="File containing hits and length data")
    parser.add_argument('-id','--ID_file', help="File containing 16S ID information")
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
