#!/usr/bin/env python
###############################################################################
#
# __get_gff_annotations__.py - Pull genome annotations information from gff files for each genome
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
import glob 
import math

from multiprocessing import Pool
from subprocess import Popen, PIPE

import os
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
    
    """
    Grab annotations for lgt sequences from genome gff files
    """
    # global variables
    transfers_dict = {}
    contig_annotations = {}
    
    
    # input files from directory
    listing = glob.glob('%s/*/*.gff' % args.genome_dir)
    
    #-----
    
    # Read in gut_oral_contigs_UniVec_blast_cleaned_bact_removed.csv
    with open(args.input_file,"r") as fh:
        header = fh.readline().rstrip()
        for l in fh:
            tabs = l.split("\t")
            img_id_a = tabs[0].rstrip()
            genome_tree_id_a = tabs[1].rstrip()
            contig_a = tabs[2].rstrip()
            contig_length_a = tabs[3].rstrip()
            start_a = tabs[4].rstrip()
            stop_a = tabs[5].rstrip()
            length_a = tabs[6].rstrip()
            img_id_b = tabs[7].rstrip()
            genome_tree_id_b =tabs[8].rstrip()
            contig_b = tabs[9].rstrip()
            contig_length_b = tabs[10].rstrip()
            start_b = tabs[11].rstrip()
            stop_b = tabs[12].rstrip()
            length_b = tabs[13].rstrip()
            # Add genome and contig data to dictionary
            try: 
                transfers_dict[img_id_a][contig_a] += [[start_a,stop_a,img_id_b]]
            except KeyError:
                transfers_dict[img_id_a] = {contig_a:[[start_a,stop_a,img_id_b]]}
            try: 
                transfers_dict[img_id_b][contig_b] += [[start_b,stop_b,img_id_a]]
            except KeyError:
                transfers_dict[img_id_b] = {contig_b:[[start_b,stop_b,img_id_a]]}
    #-----
        
    for g_file in listing:
        #grab img id from file name
        img_id = g_file.split("/")[-1].split(".")[0]
        #Is ID in dictionary
        count = 0
        if img_id in transfers_dict:
            # tab delimited
            """ ##gff-version 3
            NZ_ABQU01000001 img_er_v350     CDS     1       1237    .       -       0       ID=643923893;locus_tag=HpulM9_010100000005;product=hypothetical protein"""
            with open(g_file,"r") as fh:
                # capture header
                header = fh.readline()
                #print header
                for l in fh:
                    tabs = l.split("\t")
                    contig = tabs[0].rstrip()
                    start = int(tabs[3].rstrip())
                    stop = int(tabs[4].rstrip())
                    annotation = str(tabs[8].rstrip())
                    try:
                        contig_annotations[img_id][contig] += [[start,stop,annotation]]
                    except KeyError:
                        try:
                            contig_annotations[img_id][contig] = [[start,stop,annotation]]
                        except KeyError:                       
                            contig_annotations[img_id] = {contig:[[start,stop,annotation]]}
                #print contig_annotations
    
    #-----
    
    for id in transfers_dict:
        for contig in transfers_dict[id]:
            for i in range(len(transfers_dict[id][contig])):
                start = int(transfers_dict[id][contig][i][0])
                stop = int(transfers_dict[id][contig][i][1])
                
                diff = math.fabs(start - stop)
                lgt_annotations = []
                if start > stop:
                    start, stop = stop, start
                print "###"+"\t"+id +"\t" + contig + "\t" + str(start) + "\t" + str(stop) +"\t" + str(diff) + "\t"+ transfers_dict[id][contig][i][2]
                # check to see where start and stop finish
                try: 
                    for j in contig_annotations[id][contig]:
                        """ [152, 547, 'ID=650966343;locus_tag=HMPREF0538_22255;product=replication protein']"""
                        #print j
                        orf_start = j[0]
                        orf_end = j[1] 
                        if (orf_end > start and orf_end < stop) or (orf_start < stop and orf_start > start):
                            #print "yoyoyoyo"
                            print "\t".join([str(j[0]),
                                                 str(j[1]),
                                                 j[2]
                                                 ]) 
                except KeyError:
                    pass
                print "@@@"
    
    """                  
    for key in contig_annotations:
        for contig in contig_annotations[key]: 
            for i in contig_annotations[key][contig]:
                print "\t".join([key,
                                 contig,
                                 i[0],
                                 i[1],
                                 i[2]
                                 ])
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
    parser.add_argument('input_file', help="Provide this file: gut_oral_contigs_UniVec_blast_cleaned_bact_removed.csv")
    parser.add_argument('genome_dir', help="Directory containing img genomes")
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
