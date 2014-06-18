#!/usr/bin/env python
###############################################################################
#
# __parse_coords_methylase__.py - Parse nucmer output file, and grab information regarding methylase gene transfers!
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

class NucMerParser(object):
    """Wrapper class for parsing nucmer output"""
    # constants to make the code more readable
    _START_1  = 0
    _END_1    = 1
    _START_2  = 2
    _END_2    = 3
    _LEN_1    = 4
    _LEN_2    = 5
    _IDENTITY = 6
    _ID_1     = 7
    _ID_2     = 8

    def __init__(self):
        self.prepped = False

    def reset(self):
        self.prepped = False

    def readNuc(self, fp):
        """Read through a nucmer coords file

        this is a generator function
        """
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fp: # search for the first record
                        if l[0] == '=': # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fp:
                fields = l.split('|')
                yield ([int(i) for i in fields[0].split()] +
                       [int(i) for i in fields[1].split()] +
                       [int(i) for i in fields[2].split()] +
                       [float(i) for i in fields[3].split()] +
                       fields[4].split())
            break # done!

class faaParser(object):
    def __init__(self,accession):
        self.readFAA(accession)
    
    def readFAA(self,accession):
        dashes          = accession.rstrip().split("-")
        self.uid        = dashes[-1]
        self.img_id_a   = dashes[1].split(":")[-1]
        self.img_id_b   = dashes[5].split(":")[-1]

class IMGmetadata(object):
    def __init__(self,l):
        self.readMetadata(l)
        
    def readMetadata(self,l):
        """ taxon  = 0, genome_name = 4"""
        tabs = l.rstrip().split("\t")
        self.img_id = tabs[0]
        self.genome_name = tabs[4]
    
class methylaseGenesDB(object):
    def __init__(self):
        self.methylase_dict = {}
        self.lgt_dict       = {}
        self.metadata_dict  = {}
        
    def addLGTdata(self,accession):
        line = faaParser(accession)
        self.lgt_dict[line.uid] = [line.img_id_a,line.img_id_b]
    
    def transfered_methylase_genes(self,lgt_id,rebase_id,id_per,rebase_len_1,lgt_len_):
        self.methylase_dict[lgt_id] = [rebase_id,id_per,rebase_len_1,lgt_len_2]
       
    def addMetadata(self,l):
        line = IMGmetadata(l)
        for lgt_id in self.lgt_dict.keys():
            if self.lgt_dict[lgt_id][0] == line.img_id:
                self.metadata_dict[self.lgt_dict[lgt_id][0]] = line.genome_name
            if self.lgt_dict[lgt_id][1] == line.img_id:
                self.metadata_dict[self.lgt_dict[lgt_id][1]] = line.genome_name
                
    def printOUT(self):
        for lgt in self.methylase_dict.keys():
            try:
                rebase      = self.methylase_dict[lgt][0] 
                img_id_a    = self.lgt_dict[lgt][0]
                img_id_b    = self.lgt_dict[lgt][1] 
                genome_a    = self.metadata_dict[img_id_a]
                genome_b    = self.metadata_dict[img_id_b] 
                id_perc     = self.methylase_dict[lgt][1]
                rebase_len  = self.methylase_dict[lgt][2]
                lgt_len     = self.methylase_dict[lgt][3]
                print "\t".join([lgt,rebase,id_perc,img_id_a,genome_a,img_id_b,genome_b,rebase_len,lgt_len])
            except KeyError:
                pass

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
    NP = NucMerParser() # call class
    METHYL = methylaseGenesDB()
    
    # open .coords file
    with open(args.coords_file,"r") as fh:
        for hit in NP.readNuc(fh):
            lgt_id      = hit[NP._ID_2].split("_")[0] 
            rebase_id   = hit[NP._ID_1]
            id_perc     = hit[NP._IDENTITY]
            lgt_len     = hit[NP._LEN_2]
            rebase_len  = hit[NP._LEN_1]
            METHYL.transfered_methylase_genes(lgt_id, rebase_id,id_perc,rebase_len,lgt_len)
    
    # read in fasta file containing unique id information
    # parse fasta file using biopython
    for accession,sequence in SeqIO.to_dict(SeqIO.parse(args.fasta_file,"fasta")).items():
        METHYL.addLGTdata(accession)
    
    # read in metadata file
    with open(args.metadata,"r") as fh:
        header = fh.readline() # capture header
        for l in fh:
            METHYL.addMetadata(l)
            
    METHYL.printOUT()      
    
    
    

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
    parser.add_argument('-c','--coords_file', help="...")
    parser.add_argument('-f','--fasta_file', help="...")
    parser.add_argument('-m','--metadata', help="...")
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
