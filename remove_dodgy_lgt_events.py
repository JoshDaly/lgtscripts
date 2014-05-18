#!/usr/bin/env python
###############################################################################
#
# __script_name__.py - description!
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
    with open(args.input_file,"r") as fh:
        header = fh.readline().rstrip()
        print header
        for l in fh:
            "img_id_a    genome_tree_id_a    contig_a    contig_length_a    start_a    stop_a    length_a    img_id_b    genome_tree_id_b    contig_b    contig_length_b    start_b    stop_b    length_b"
            img_id_a = l.split("\t")[0].rstrip()
            genome_tree_id_a = l.split("\t")[1].rstrip()
            contig_a = l.split("\t")[2].rstrip()
            contig_length_a = l.split("\t")[3].rstrip()
            start_a = l.split("\t")[4].rstrip()
            stop_a = l.split("\t")[5].rstrip()
            length_a = l.split("\t")[6].rstrip()
            img_id_b = l.split("\t")[7].rstrip()
            genome_tree_id_b =l.split("\t")[8].rstrip()
            contig_b = l.split("\t")[9].rstrip()
            contig_length_b = l.split("\t")[10].rstrip()
            start_b = l.split("\t")[11].rstrip()
            stop_b = l.split("\t")[12].rstrip()
            length_b = l.split("\t")[13].rstrip()
            
            #print "\t".join([img_id_a,
            #                 img_id_b])
            
            if (img_id_a == "637000025" or img_id_a == "646206273" or img_id_a == "647533110" or img_id_a == "637000024") and (img_id_b == "637000025" or img_id_b == "646206273" or img_id_b == "647533110" or img_id_b == "637000024"):
                pass 
                #print "I see you!"
            else:
                print "\t".join([img_id_a,
                                 genome_tree_id_a,
                                 contig_a,
                                 contig_length_a,
                                 start_a,
                                 stop_a,
                                 length_a,
                                 img_id_b,
                                 genome_tree_id_b,
                                 contig_b,
                                 contig_length_b,
                                 start_b,
                                 stop_b,
                                 length_b                                 
                                 ])
                
        #637000025
        #646206273
        #647533110
        #637000024
                
                
            
            
            
            
            
    
    

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
    parser.add_argument('input_file', help="gut oral file")
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
