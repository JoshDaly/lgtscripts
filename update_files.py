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
from Bio import SeqIO
from Bio.Seq import Seq

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
    interacting_genomes={}
    
    # open in interacting genomes list
    with open(args.genome_list,"r") as fh:
        # no header
        for l in fh:
            id = l.rstrip()
            interacting_genomes[id]=0
    genomes = {}
    
    # read in fasta file
    for accession,sequence in SeqIO.to_dict(SeqIO.parse(args.fasta,"fasta")).items():
        print accession
        #else:
            #print accession
            #genomes_dict[accession] = [len(sequence),img_id, sequence.seq]
    
    
    
    
    
    
    # Read in input file for formatting
    with open(args.input_file, "r") as fh:
        # capture header
        header= fh.readline()
        for l in fh:
            """img_id_a        genome_tree_id_a        contig_a        contig_length_a start_a stop_a  length_a        img_id_b        genome_tree_id_b        contig_b        contig_length_b start_b stop_b  length_b"""
            tabs = l.split("\t")
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
            try:
                genomes[img_id_a]+=1
            except KeyError:
                genomes[img_id_a]=1
            try:
                genomes[img_id_b]+=1
            except KeyError:
                genomes[img_id_b]=1
            
            #if img_id_a in interacting_genomes and img_id_b in interacting_genomes:
            #interacting_genomes[img_id_a]+=1
            #interacting_genomes[img_id_b]+=1
            """
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
                """
        #for key in genomes.keys():
        #    print key +"\t"+str(genomes[key])
        #for key in interacting_genomes.keys():
        #    print key + "\t"+ str(interacting_genomes[key])
            
            
            
            
    
    

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
    parser.add_argument('-i','--input_file', help="...")
    parser.add_argument('-g','--genome_list', help="List of interacting genomes")
    parser.add_argument('-f','--fasta', help="File containing transferred sequences in fasta format")
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
