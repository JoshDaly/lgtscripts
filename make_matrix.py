#!/usr/bin/env python
###############################################################################
#
# __make_matrix__.py - make a matrix from a tab delimited file!
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

import numpy as np
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
    
    """ read in the tab delimited file
    img_id_a        contig_a        img_id_b        contig_b        hits    length
    """ 
    ids_dict = {} #dictionary to store ids and position 
    ids_list= [] #list to store 
    matrix = [] #matrix containing both hits and length
    matrix_hits = [] # matrix for hits
    matrix_length = [] # matrix for length
    #read through file line by line

    with open(args.input_file, "r") as fh:
        
        #capture file header
        header = fh.readline()
        
        #read through each line
        for l in fh: 
            # assign each column in line to a variable
            id_a = l.split('\t')[0]
            id_b = l.split('\t')[2]
            contig_a = l.split('\t')[1]
            contig_b = l.split('\t')[3]
            hits = l.split('\t')[4]
            length = l.split('\t')[5]
        
            # check if id_a in dictionary
            if ids_dict.has_key(id_a) is False: 
                ids_dict[id_a] = len(ids_list)
                ids_list.append(id_a)       
        
             #check if id_b in dictionary
            if ids_dict.has_key(id_b) is False: 
                ids_dict[id_b] = len(ids_list)
                ids_list.append(id_b) 
        
        #-----
              
        # create matrix of 0s
        for i in range(0,len(ids_list)):
            matrix.append(['0']*len(ids_list)) 
            matrix_hits.append(['0']*len(ids_list))
            matrix_length.append(['0']*len(ids_list))
                       
        """ Reset to the beginning of the file"""
        #reset file pointer to beginning of file
        fh.seek(0)
    
        #capture file header
        header = fh.readline()
        
        #read through file again from the beginning    
        for l in fh:
             # assign each column in line to a variable
            id_a = l.split('\t')[0]
            id_b = l.split('\t')[2]
            contig_a = l.split('\t')[1]
            contig_b = l.split('\t')[3]
            hits = l.split('\t')[4]
            length = l.split('\t')[5] 
            
            # Check to see which transformation has been chosen
            if args.hits_transformation == "loge":
                hits = np.log(int(hits))
            if args.hits_transformation == "log10":
                hits = np.log10(int(hits))
            if args.hits_transformation == "sqrt":
                hits = np.sqrt(int(hits))
            if args.length_transformation == "loge":
                length = np.log(int(length))
            if args.hits_transformation == "log10":
                hits = np.log10(int(length))
            if args.length_transformation == "sqrt":
                length = np.sqrt(int(length))
            
            # for each ids pair, pull location from ids_dict
            # hits on top half of matrix triangle
            # length on bottom half of matrix triangle            
            matrix[min(ids_dict[id_a],ids_dict[id_b])][max(ids_dict[id_a],ids_dict[id_b])] = str(hits).rstrip() # strip line endings using rstrip()
            matrix[max(ids_dict[id_a],ids_dict[id_b])][min(ids_dict[id_a],ids_dict[id_b])] = str(length).rstrip()
            matrix_hits[min(ids_dict[id_a],ids_dict[id_b])][max(ids_dict[id_a],ids_dict[id_b])] = str(hits).rstrip() # strip line endings using rstrip()
            matrix_length[max(ids_dict[id_a],ids_dict[id_b])][min(ids_dict[id_a],ids_dict[id_b])] = str(length).rstrip()
               
        
        #-----
        
        #open output file for combined data matrix
        with open(args.output_file,'w')as fh2: 
            #write to the output_file, and always use + instead of , when printing/writing to a file. 
            fh2.write("img_ids\t"+"\t".join(ids_list)+"\n")
            for i in range(0,len(matrix)):
                fh2.write(ids_list[i]+"\t"+"\t".join(matrix[i])+"\n")
                
        # open output for hits only        
        with open(args.output_file2,'w')as fh2: 
            #write to the output_file, and always use + instead of , when printing/writing to a file. 
            fh2.write("img_ids\t"+"\t".join(ids_list)+"\n")
            for i in range(0,len(ids_list)):
                fh2.write(ids_list[i]+"\t"+"\t".join(matrix_hits[i])+"\n")
        # open output for length only
        with open(args.output_file3,'w')as fh2: 
            #write to the output_file, and always use + instead of , when printing/writing to a file. 
            fh2.write("img_ids\t"+"\t".join(ids_list)+"\n")
            for i in range(0,len(ids_list)):
                fh2.write(ids_list[i]+"\t"+"\t".join(matrix_length[i])+"\n")
        

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
    parser.add_argument('-o','--output_file', help="Output file 1")
    parser.add_argument('-o2','--output_file2', help="Output file 2 containing hits")
    parser.add_argument('-o3','--output_file3', help="Output file 3 containing length")
    parser.add_argument('-ht','--hits_transformation', help="Data transformation for hits e.g. loge")
    parser.add_argument('-lt','--length_transformation', help="Data transformation for cumulative contig length e.g. loge")   
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
