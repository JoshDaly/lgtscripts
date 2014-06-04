#!/usr/bin/env python
###############################################################################
#
# __split-josh_1_per_line.tsv-sep-fasta-files__.py - Split file Sam produced containing multiple 16S sequences into multiple fasta files!
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

import os
import errno

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

class capture16S( object ):
    def __init__(self):
        self.dict_16S = {}
        
    def addGenome16S( self,img_id,seq_16S ):
        self.dict_16S[img_id] = seq_16S
        
    def checkSeqSize(self,seq_16S):
        if len(seq_16S) > 400:
            return True
    def printDict(self):
        for img_id in self.dict_16S.keys():
            print "\t".join([img_id,self.dict_16S[img_id]])
            
    def printToFile(self,):
        for img_id in self.dict_16S.keys():
            output_file = "%s_16S.fna" % img_id
            output_dir = os.path.join(args.output_directory,"%s" % img_id)
            FULL_PATH = os.path.join(args.output_directory,output_file)
            
            if not os.path.exists(output_dir):
                os.system("mkdir %s" % (output_dir))
                
            with open(FULL_PATH,"w") as fh:
                fh.write(">"+img_id+"\n")
                fh.write(self.dict_16S[img_id])
        
        
        
class genomeInfo(object):
    def __init__(self):
        self.dict_genome_tree = {}
        self.dict_genome_name = {}
        
    def addGenomeTree(self,img_id,genome_id):
        self.dict_genome_tree[img_id] = genome_id
    
    def addGenomeName(self, img_id, genome_name):
        self.dict_genome_name[img_id] = genome_name
        

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
    """genome_id    img_id    Genome_name    16Sseq1    seq2    seqn"""
    # objects
    dict_16S = capture16S() # hold 16S sequence
    dict_info = genomeInfo() # hold additional info
    count = 0
    
    # read in file
    with open(args.tsv_file,"r") as fh:
        for l in fh:
            tabs = l.split("\t")
            genome_id = tabs[0]
            img_id = tabs[1]
            genome_name = tabs[2]
            seq_16S = tabs[3:]
            if count < 10:
                if len(seq_16S) > 1:
                    for seq in seq_16S:
                        seq = seq.rstrip()
                        if dict_16S.checkSeqSize(seq):
                            dict_16S.addGenome16S(img_id, seq)
                            break 
            dict_info.addGenomeTree(img_id, genome_id)
            dict_info.addGenomeName(img_id, genome_name)
            #count+=1 
    dict_16S.printToFile() # print out to separate fna files
                
                
            
            
            
            
    """
# parse fasta file using biopython
for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
if accession in genomes_dict:
pass
else:
#print accession
genomes_dict[accession] = [len(sequence),img_id, sequence.seq
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
    parser.add_argument('-tsv','--tsv_file', help="...")
    parser.add_argument('-o','--output_directory', help="...")
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
