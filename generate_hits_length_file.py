#!/usr/bin/env python
###############################################################################
#
# __generate_hits_length_file__.py - Generate a hits-length file from a transfers.csv file
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
    # Objects
    transfers_dict = {} 
    hits_accumulative= 0
    
    #parse gut_oral_contigs.csv file
    with open(args.transfers_file,"r") as fh:
        #header
        header=fh.readline().rstrip()
        #Output style
        """genome_tree_id_a        img_id_a        bodysite_a      genome_tree_id_b        img_id_b        bodysite_b      hits    length"""
        #read through file line by line
        #print header
        for l in fh:
            Good_to_print = True
            tabs= l.split("\t")
            img_id_a= tabs[0]
            genome_tree_id_a= tabs[1]
            contig_a= tabs[2]
            contig_length_a= tabs[3]
            start_a= tabs[4]
            stop_a= tabs[5]
            length_a= tabs[6]
            img_id_b= tabs[7]
            genome_tree_id_b= tabs[8]
            contig_b= tabs[9]
            contig_length_b= tabs[10]
            start_b= tabs[11]
            stop_b= tabs[12]
            length_b= tabs[13].rstrip()    
            
            
            length = (int(length_a) + int(length_b))/2                    
            try:
                # transfers_dict[img_id_a] = {img_id_b:[hits,length,genome_tree_id_a,genome_tree_id_b]}
                transfers_dict[img_id_a][img_id_b][0]+=1
                transfers_dict[img_id_a][img_id_b][1]+=((int(length_a)+int(length_b))/2)
            except KeyError:
                try:
                    transfers_dict[img_id_a][img_id_b] = [1,length,genome_tree_id_a,genome_tree_id_b]
                except KeyError:
                    try: 
                        transfers_dict[img_id_a]= {img_id_b:[1,length,genome_tree_id_a,genome_tree_id_b]}
                    except KeyError:
                        pass
    
    #-----
    #header
    #print "\t".join(["genome_tree_id_a",
    #                 "img_id_a",
    #                 "genome_tree_id_b",
    #                 "img_id_b",
    #                 "hits",
    #                 "length"])
    
    
    #print transfers_dict
    for id in transfers_dict:
        for id_b in transfers_dict[id]:
            genome_tree_a= transfers_dict[id][id_b][2]
            genome_tree_b= transfers_dict[id][id_b][3]
            hits= str(transfers_dict[id][id_b][0])
            length= str(transfers_dict[id][id_b][1])
            hits_accumulative= int(hits) + hits_accumulative
            #print "\t".join([genome_tree_a,
            #                 id,
            #                 genome_tree_b,
            #                 id_b,
            #                 hits,
            #                 length
            #                 ])
    print hits_accumulative
                
                
            
            
            
            
    
    

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
    parser.add_argument('-i','--transfers_file', help="Transfers.csv file")
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
