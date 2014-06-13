#!/usr/bin/env python
###############################################################################
#
# __add_genome_name_to_kmer_full_info__.py - Add genome species name to kmer_full_info.10.06.14.csv file
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

from multiprocessing import Pool
from subprocess import Popen, PIPE

#from Bio import SeqIO
#from Bio.Seq import Seq

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

class kmerFileParser(object):
    # constants
    _LGT        = 0
    _IMG_ID_1   = 1
    _IMG_ID_2   = 2
    _KMER_SCORE = 3
    _MEAN_DGG   = 4
    _DG_1       = 5
    _DG_2       = 6
    
    def __init__(self):
        self.prepped = False
    
    def readKmerFile(self,fh):
        line = None
        while True:
            if not self.prepped:
                for l in fh:
                    if l[0:3]=="lgt":
                        self.prepped = True
                        break
            for l in fh:
                fields = l.split("\t")
                yield [fields[0],
                       fields[1],
                       fields[2],
                       fields[3],
                       fields[4],
                       fields[5],
                       fields[6]]
            break # done! 
        
class kmerdata(object):
    def __init__(self):
        self.kmer_file_dict = {} # lgt id -> 
        self.genome_db = {}
        
    def addGenome(self,img_id,genome_name):
        self.genome_db[img_id] = genome_name
        
    def addLGTkmer(self,lgt_id,id_a,id_b,kmer_score,dgg,dg1,dg2):
        self.kmer_file_dict[lgt_id] = [id_a,id_b,kmer_score,dgg,dg1,dg2]

    def printHeader(self):
        print "\t".join(["lgt","img_id_a","genome_a","img_id_b","genome_b","kmer_score","mean_dg","dg1","dg2"])
        
    def printINFO(self):
        for lgt in self.kmer_file_dict.keys():
            id_a        = self.kmer_file_dict[lgt][0]
            id_b        = self.kmer_file_dict[lgt][1]
            kmer_score  = self.kmer_file_dict[lgt][2]
            dgg         = self.kmer_file_dict[lgt][3]
            dg1         = self.kmer_file_dict[lgt][4]
            dg2         = self.kmer_file_dict[lgt][5]
            genome_1    = self.genome_db[id_a]
            genome_2    = self.genome_db[id_b]
            print "\t".join(lgt,id_a,genome_1,id_b,genome_2,kmer_score,dgg,dg1,dg2)

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
    KP = kmerFileParser()
    KMER = kmerdata()
    
    #-----
    with open(args.metadata_file,"r") as fh:
        header = fh.readline() # capture header
        for l in fh:
            tabs = l.split("\t")
            img_id = tabs[0]
            genome_name = tabs[4]
            KMER.addGenome(img_id, genome_name)
    
    #-----
    with open(args.kmer_file,"r") as fh:
        for hit in KP.readKmerFile(fh):
            lgt_id      = hit[KP._LGT]
            id_a        = hit[KP._IMG_ID_1]
            id_b        = hit[KP._IMG_ID_2]
            kmer_score  = hit[KP._KMER_SCORE]
            dgg         = hit[KP._MEAN_DGG]
            dg1         = hit[KP._DG_1]
            dg2         = hit[KP._DG_2]
            KMER.addLGTkmer(lgt_id, id_a, id_b, kmer_score, dgg, dg1, dg2)
    KMER.printHeader()
    KMER.printINFO()
            
            
            
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
    parser.add_argument('-metadata','--metadata_file', help="...")
    parser.add_argument('-kmer','--kmer_file', help="...")
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
