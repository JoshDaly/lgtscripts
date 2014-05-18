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
    parsing = False
    transfer_annotations = {}

    # read in file containing annotations
    with open(args.input_file,"r") as fh:
        # no header
        header_data= {}
        for l in fh:
            if "###" in l:
                parsing = True
                #parsing = not parsing
                #continue.
                tabs = l.split("\t")
                #print "\t".join([tabs[0],tabs[1],tabs[2],tabs[3],tabs[4],tabs[5],tabs[6]])
                id_a = tabs[1].rstrip()
                contig = tabs[2].rstrip() 
                start = tabs[3].rstrip()
                stop = tabs[4].rstrip()
                length = tabs[5].rstrip()
                id_b =  tabs[6].rstrip()
                header_data['id_a'] = id_a
                header_data['id_b'] = id_b
                header_data['contig'] = contig
                
            if parsing and "###" not in l and "@@@" not in l:
                tabs_2 = l.split("\t")
                start = tabs_2[0].rstrip()
                stop = tabs_2[1].rstrip()
                annotation = tabs_2[2].rstrip()
                try:
                    transfer_annotations[header_data['id_a']][header_data['id_b']][header_data['contig']] +=  [[start,stop,annotation]]
                except KeyError:
                    try:
                        transfer_annotations[header_data['id_a']][header_data['id_b']] = {header_data['contig']:[[start,stop,annotation]]}
                    except KeyError:
                        try: 
                            transfer_annotations[header_data['id_a']] += {header_data['id_b']:{header_data['contig']:[[start,stop,annotation]]}}
                        except KeyError:
                            transfer_annotations[header_data['id_a']] = {header_data['id_b']:{header_data['contig']:[[start,stop,annotation]]}}        
            if "@@@" in l:
                parsing = False
                header_data = {}
    #------
    # objects
    phage_trans = {}
    plasmid_trans = {}
    
    for id_a in transfer_annotations:
        for id_b in transfer_annotations[id_a]:
            for contig in transfer_annotations[id_a][id_b]:
                for i in transfer_annotations[id_a][id_b][contig]:
                    #case insensitive search
                    if "phage" in i[2].lower():
                        
                        
                         
    
            
                
        
        
                    
                    
                    
                    
                    
            
            
            
            
            
    
    

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
    parser.add_argument('input_file', help="Provide file containing transfer annotations")
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
