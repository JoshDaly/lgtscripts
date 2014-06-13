#!/usr/bin/env python
###############################################################################
#
# __vectors_vs_seq_platform__.py - Compare sequencing vector contamination to sequencing platform used.
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

class TransferParser(object):
    """Wrapper class for parsing transfer files"""
    #constants to make the code readable
    _UID_1          = 0
    _IMG_ID_1       = 1
    _GT_ID_1        = 2
    _CONTIG_1       = 3
    _CONTIG_LEN_1   = 4
    _START_1        = 5
    _STOP_1         = 6
    _LEN_1          = 7
    _SEQ_PLAT_1     = 8
    _UID_2          = 9
    _IMG_ID_2       = 10
    _GT_ID_2        = 11
    _CONTIG_2       = 12
    _CONTIG_LEN_2   = 13
    _START_2        = 14
    _STOP_2         = 15
    _LEN_2          = 16
    _SEQ_PLAT_2     = 17
    _16S_ID         = 18
    
    def __init__(self):
        self.prepped = False
        
    def reset(self):
        self.prepped = False
    
    def readTrans(self,fh):
        
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                for l in fh: # search for the first record
                    if l[0:3] =="uid": # next line is good
                        #print l
                        self.prepped = True
                        break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield [fields[0],
                       fields[1],
                       fields[2],
                       fields[3],
                       fields[4],
                       fields[5],
                       fields[6],
                       fields[7],
                       fields[8],
                       fields[9],
                       fields[10],
                       fields[11],
                       fields[12],
                       fields[13],
                       fields[14],
                       fields[15],
                       fields[16],
                       fields[17],
                       fields[18]]
            break # done!  

class TransferDB(object):
    def __init__(self):
        self.dirty_transfers_dict = {}
        self.dirty_seq_platform = {}
        self.clean_transfers_dict = {}
        self.clean_seq_platform = {}
        self.platform_dirty = {}
        self.platform_clean = {}
        
    def addDirtyTransfer(self,img_id):
        try:
            self.dirty_transfers_dict[img_id] += 1
        except KeyError:
            self.dirty_transfers_dict[img_id] = 1
    
    def addCleanTransfer(self,img_id):
        try:
            self.clean_transfers_dict[img_id] += 1
        except KeyError:
            self.clean_transfers_dict[img_id] = 1
        
    def addDirtyPlatform(self,img_id,platform):
        if "454" in platform: # collate 454 platforms
            platform = "454"
        self.dirty_seq_platform[img_id] = platform
        
    def addCleanPlatform(self,img_id,platform):
        if "454" in platform: # collate 454 platforms
            platform == "454"
        self.clean_seq_platform[img_id] = platform
        
    def compareDicts(self):
        for id in self.dirty_transfers_dict.keys():
            try: 
                diff = 1 - float(self.clean_transfers_dict[id])/self.dirty_transfers_dict[id]
                print "\t".join([id,str(self.dirty_transfers_dict[id]),str(self.clean_transfers_dict[id]),str(diff),self.dirty_seq_platform[id]])
            except KeyError:
                diff = 1
                print "\t".join([id,str(self.dirty_transfers_dict[id]),str(0),str(diff),self.dirty_seq_platform[id]])
        
    def printDict(self):
        for id in self.dirty_transfers_dict.keys():
            print "\t".join([id, str(self.dirty_transfers_dict[id])])
    
    def returnHits(self,id):
        if id in self.dirty_transfers_dict:
            return self.dirty_transfers_dict[id]
        else:
            return 0
        
    def returnIDs(self):
        for id in self.dirty_transfers_dict.keys():
            print id
        
    def collatePlatforms(self):
        for id in self.dirty_transfers_dict:
            print self.dirty_transfers_dict[id]
            try:
                self.platform_dirty[self.dirty_seq_platform[id]] += self.dirty_transfers_dict[id]
            except KeyError:
                self.platform_dirty[self.dirty_seq_platform[id]] = self.dirty_transfers_dict[id]
            try:
                self.platform_clean[self.clean_transfers_dict[id]] += self.clean_seq_platform[id]
            except KeyError:
                self.platform_clean[self.clean_transfers_dict[id]] = self.clean_seq_platform[id]
    
    def printCollatedPlatforms(self):
        for platform in self.platform_dirty:
            try:
                print "\t".join([self.platform_dirty[platform],self.platform_clean[platform]])
            except KeyError:
                print platform
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

def DirtyVSClean(dirty_dict,clean_dict):
    try:
        print "\t".join([id,str(dirty_dict[id]),str(clean_dict)])
    except KeyError:
        print "\t".join([id,str(dirty_dict[id]),str(0)])
        
def printHeader():
    print "\t".join(["img_id","hits_pre_clean","hits_post_clean","ratio","platform"])

def doWork( args ):
    """ Main wrapper"""
    """read in two transfer files, and capture information in dictionaries"""
    """compare the vector contamination with sequencing platform"""
    # objects
    transfers_dict = TransferDB()
    TP = TransferParser()
    
    #-----
    # read in dirty transfers file
    with open(args.dirty,"r") as fh:
        for l in TP.readTrans(fh):
            transfers_dict.addDirtyTransfer(l[TP._IMG_ID_1])
            transfers_dict.addDirtyTransfer(l[TP._IMG_ID_2])
            transfers_dict.addDirtyPlatform(l[TP._IMG_ID_1], l[TP._SEQ_PLAT_1])
            transfers_dict.addDirtyPlatform(l[TP._IMG_ID_2], l[TP._SEQ_PLAT_2])
    #-----
    # read in clean transfers file
    with open(args.clean,"r") as fh: 
        for l in TP.readTrans(fh):
            transfers_dict.addCleanTransfer(l[TP._IMG_ID_1])
            transfers_dict.addCleanTransfer(l[TP._IMG_ID_2])
            transfers_dict.addCleanPlatform(l[TP._IMG_ID_1], l[TP._SEQ_PLAT_1])
            transfers_dict.addCleanPlatform(l[TP._IMG_ID_2], l[TP._SEQ_PLAT_2])
    printHeader()
    transfers_dict.collatePlatforms()
    #transfers_dict.printCollatedPlatforms()
    #transfers_dict.compareDicts()
        
            
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
    parser.add_argument('-dirty','--dirty', help="Transfers file")
    parser.add_argument('-clean','--clean', help="Transfers file")
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
