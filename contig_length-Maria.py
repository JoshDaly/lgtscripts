#!/usr/bin/env python
###############################################################################
#
# __contig_length__.py - Grab sequence information from a directory of fna files.
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
import datetime

from multiprocessing import Pool
from subprocess import Popen, PIPE
import glob 
from Bio import SeqIO
from Bio.Seq import Seq

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
    """Read in a directory containing multiple fna files, then file by file count the length of contigs.
    Then write this information to a file""" 
    #global variables    
    #Parse many fna files using glob
    #print "start", datetime.datetime.now()
    listing = glob.glob('%s/*.fna' % args.genomes_directory)
    genomes_dict = {} #dictionary to store index information of genomes
    #gut_oral_contigs = {} #dictionary to store information from gut_oral_contigs.csv
    
    #gut_oral_contigs = {} # dictionary to store information regarding contig length
    count = 0
    
    for c_file in listing:
        img_id = c_file.split("/")[2].split(".")[0] 
        if count <1:
            #print c_file
            for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
                if accession in genomes_dict:
                    #print "Hello Pilgrim"
                    pass
                else:
                    #print accession
                    genomes_dict[accession] = [len(sequence),img_id, sequence.seq]
        #for i in genomes_dict.keys():
        #    print i + "\t" + str(genomes_dict[i][0])
    
    #print "end", datetime.datetime.now()        
    #print header
    """
    print "\t".join(["img_id_a",
                     "genome_tree_id_a",
                     "contig_a",
                     "contig_length_a",
                     "start_a",
                     "stop_a",
                     "length_a",
                     "img_id_b",
                     "genome_tree_id_b",
                     "contig_b",
                     "contig_length_b",
                     "start_b",
                     "stop_b",
                     "length_b"
                     ])
    """        
    """   
    for key in  gut_oral_dict:
        print "\t".join([key, 
                         str(len(gut_oral_dict[key][0])),
                         gut_oral_dict[key][1]
                         ])
    """    
    """img_id_a        genome_tree_id_a        contig_a        start_a stop_a    length_a img_id_b        genome_tree_id_b        contig_b        start_b stop_b    length_b"""    
    with open(args.contig_file,"r") as fh:
        #capture header
        header = fh.readline()
        #read through line by line
        count = 0
        uniq_id = 1
        if count < 1:
            for l in fh:
                img_id_a = l.split("\t")[0].rstrip()
                genome_tree_id_a = l.split("\t")[1].rstrip() 
                contig_a = l.split("\t")[2].rstrip()
                contig_length_a = l.split("\t")[3].rstrip()
                start_a = l.split("\t")[4].rstrip()
                stop_a = l.split("\t")[5].rstrip()
                length_a = l.split("\t")[6].rstrip()
                img_id_b = l.split("\t")[7].rstrip()
                genome_tree_id_b = l.split("\t")[8].rstrip()
                contig_b = l.split("\t")[9].rstrip()
                contig_length_b = l.split("\t")[10].rstrip()
                start_b = l.split("\t")[11].rstrip()
                stop_b = l.split("\t")[12].rstrip()
                length_b = l.split("\t")[13].rstrip()
                
                #print "\t".join([str(int(start_a)+1),str(int(stop_a)+1)])
                
                try:
                    if int(start_a) > int(stop_a):
                        new_start_a = stop_a
                        new_stop_a = start_a
                        
                        #contig_length_a = str(genomes_dict[contig_a][0])
                        contig_seq_a = genomes_dict[contig_a][2]
                        #contig_length_b = str(genomes_dict[contig_b][0])
                        print ">"+contig_a+"-img_a:"+img_id_a+"-genome_tree_a:"+genome_tree_id_a+"-"+"Start:"+str(start_a)+"-"+"Stop:"+str(stop_a)+"-id_b:"+img_id_b+"-genome_tree_b:"+genome_tree_id_b+"-"+str(uniq_id)
                        #print contig_seq_a[(int(start_a)-1):(int(stop_a)+1)]
                        print contig_seq_a[int(new_start_a)-1:int(new_stop_a)+1]
                    else:
                        contig_seq_a = genomes_dict[contig_a][2]
                        #contig_length_b = str(genomes_dict[contig_b][0])
                        print ">"+contig_a+"-img_a:"+img_id_a+"-genome_tree_a:"+genome_tree_id_a+"-"+"Start:"+str(start_a)+"-"+"Stop:"+str(stop_a)+"-id_b:"+img_id_b+"-genome_tree_b:"+genome_tree_id_b+"-"+str(uniq_id)
                        #print contig_seq_a[(int(start_a)-1):(int(stop_a)+1)]
                        print contig_seq_a[int(start_a)-1:int(stop_a)+1]
                except KeyError:
                    pass
                
                uniq_id += 1
                
                try:
                    if int(start_b) > int(stop_b):
                        new_start_b = stop_b
                        new_stop_b = start_b
                        
                        contig_seq_b = genomes_dict[contig_b][2]
                        print ">"+contig_b+"-img_a:"+img_id_b+"-genome_tree_a:"+genome_tree_id_b+"-"+"Start:"+str(start_b)+"-"+"Stop:"+str(stop_b)+"-id_b:"+img_id_a+"-genome_tree_b:"+genome_tree_id_a+"-"+str(uniq_id)
                        #print contig_seq_b[(int(start_b)-1):(int(stop_b)+1)]
                        print contig_seq_b[int(new_start_b)-1:int(new_stop_b)+1]
                    else:
                        contig_seq_b = genomes_dict[contig_b][2]
                        print ">"+contig_b+"-img_a:"+img_id_b+"-genome_tree_a:"+genome_tree_id_b+"-"+"Start:"+str(start_b)+"-"+"Stop:"+str(stop_b) +"-id_b:"+img_id_a+"-genome_tree_b:"+genome_tree_id_a+"-"+str(uniq_id)
                        #print contig_seq_b[(int(start_b)-1):(int(stop_b)+1)]
                        print contig_seq_b[int(start_b)-1:int(stop_b)+1]
                except KeyError:
                    pass
                uniq_id += 1
                
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
                
                
                #count +=1 
    #print "fin", datetime.datetime.now()
            
                
    """
    for key in gut_oral_dict.keys():
        seq_record = gut_oral_dict[key]
        print seq_record 
       
        print "\t".join([key,
                        len(seq_record)
                        ])
    """
    # read through the file
    """
    for seq_file in SeqIO.parse(listing,"fasta"):
        print (seq_file.id)
        print (repr(seq_file.seq))
        print (len(seq_file))
        break
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
    parser.add_argument('genomes_directory', help="File containing genomes in fna format")
    parser.add_argument('contig_file', help="gut_oral_contigs.csv")
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
