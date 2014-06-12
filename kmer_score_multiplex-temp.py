#!/usr/bin/env python
###############################################################################
#
# __kmer_score_multiplex__.py - Determine the kmer scores for each LGT event!
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
import datetime

from multiprocessing import Pool
from subprocess import Popen, PIPE

from Bio import SeqIO
from Bio.Seq import Seq

from scipy.spatial.distance import pdist

#import os
#import errno

import numpy as np
np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class LGTInfoStore( object ):
    """Helper class for storing information about directionality of LGT events"""
    # __init__ links to object, and must be present in every Class! 
    def __init__( self ):
        self.lgtTmer = {}
        self.genomeTmers = {}
        self.lgtGenomes = {}
        self.what = None
        self.Dist_dict = {} # store distances
        self.lgtScores = {}

    def addWhat( self, what ):
        self.what = what
            
    def addLGT( self, lgt_id, GID ):
        """populate lgtgenomes dictionary"""
        self.lgtScores[lgt_id] = None
        try:
            self.lgtGenomes[lgt_id] += [GID]
        except KeyError:
            self.lgtGenomes[lgt_id] = [GID] 

    def addLGTTmer( self, lgt_id, tmer ):
        self.lgtTmer[lgt_id] = tmer 
        
    def addGenomeTmer( self, GID, tmer ):
        self.genomeTmers[GID] = tmer

    def getClosestGID( self , rounded):
        """Calculate the kmer score
        
        returns (score, (closestGID, dist), (furthestGID, dist))
        """
        LGTs = self.lgtGenomes.keys()
        #print LGTs
        dgs = []
        for lgt_id in LGTs:
            #print lgt_id
            dg1 = pdist([self.genomeTmers[self.lgtGenomes[lgt_id][0]], self.lgtTmer[lgt_id] ])
            dg2 = pdist([self.genomeTmers[self.lgtGenomes[lgt_id][1]], self.lgtTmer[lgt_id] ]) 
            rounded_score = float(np.round(dg1/(dg1+dg2),decimals=2))
            score = float(dg1/(dg1+dg2))
            #print rounded_score
            self.lgtScores[lgt_id] = [score,float(np.mean([dg1,dg2])),dg1,dg2]
            print self.lgtScores
            if rounded:
                try:
                    self.Dist_dict[rounded_score]+=1
                except KeyError:
                    self.Dist_dict[rounded_score]=1
            else:
                self.Dist_dict[score]=[float(np.mean([dg1,dg2])),dg1,dg2]
            #print self.lgtScores
                
    def printInfoHeader(self):
        print "\t".join(["lgt","img_id_a","img_id_b","kmer_score","mean_dg","dg1","dg2"])
    
    def printInfo(self):
        for lgt in self.lgtScores.keys():
            kmer_score  = str(self.lgtScores[lgt][0])
            mean_dg     = str(self.lgtScores[lgt][1])
            dg1         = str(self.lgtScores[lgt][2])
            dg2         = str(self.lgtScores[lgt][3])
            g1          = self.lgtGenomes[lgt][0]
            g2          = self.lgtGenomes[lgt][1]
            print "\t".join([lgt,g1,g2,kmer_score,mean_dg,dg1,dg2])
        
        
    def getDistHisto(self):
        for score in self.Dist_dict.keys():
            print "\t".join([str(score),str(self.Dist_dict[score][0]),str(self.Dist_dict[score][1]),str(self.Dist_dict[score][2])])
        
    def printDict( self ):
        for uid in self.lgtGenomes:
            print "\t".join([uid,self.lgtGenomes[uid][0],self.lgtGenomes[uid][1]])

    def __str__( self ):
        """print function"""
        return "GID1: %s GID2: %s" % (i for i in self.genomeTmers.keys()) 
            
class TransferParser(object):
    """Wrapper class for parsing transfer files"""
    """img_id_a        genome_tree_id_a        contig_a        contig_length_a start_a stop_a  length_a        img_id_b        genome_tree_id_b        contig_b        contig_length_b start_b stop_b  length_b"""
    #constants to make the code readable
    _UID_1          = 0
    _IMG_ID_1       = 1
    _GT_ID_1        = 2
    _CONTIG_1       = 3
    _CONTIG_LEN_1   = 4
    _START_1        = 5
    _STOP_1         = 6
    _LEN_1          = 7
    _UID_2          = 8 
    _IMG_ID_2       = 9
    _GT_ID_2        = 10
    _CONTIG_2       = 11
    _CONTIG_LEN_2   = 12
    _START_2        = 13
    _STOP_2         = 14
    _LEN_2          = 15
    
    def __init__(self):
        self.prepped = False
    
    def readTrans(self,fh):
        
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fh: # search for the first record
                        if l[0:3] =="uid": # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield [fields[0],
                       fields[1],
                       fields[2],
                       fields[3],
                       int(fields[4]),
                       int(fields[5]),
                       int(fields[6]),
                       int(fields[7]),
                       fields[8],
                       fields[9],
                       fields[10],
                       fields[11],
                       int(fields[12]),
                       int(fields[13]),
                       int(fields[14]),
                       int(fields[15])
                       ]
            break # done!
        
class lgtTransfersDict( object ):
    def __init__( self ):
        self.lgtdict = {}
        
    def addLGT( self,lgt ):
        self.lgtdict[lgt] = 0
    
    def getKeys( self ):
        for uid in self.lgtdict.keys():
            print uid
        
    def checkUID(self, uid):
        if uid in self.lgtdict:
            return True
        
    
 
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

def printDict( dict ):
    for key in dict.keys():
        print "\t".join([str(key),str(dict[key])])

def printHeader():
    print "\t".join(["score","mean_distance","dg1","dg2"])

def getIDs( column ):
    genome_1= column.rstrip().split("-")[2].split(":")[1] 
    genome_2= column.rstrip().split("-")[6].split(":")[1]
    lgt_id= column.rstrip().split("-")[0]
    return (lgt_id,genome_1,genome_2)

def doWork( args ):
    """ Main wrapper"""
    # objects
    kmer_directories = glob.glob('%s/*' % args.kmers_directory)
    LGT_kmers = LGTInfoStore()      # call class
    lgt_dict =  lgtTransfersDict()  # dict of lgt events from transfers file
    TP = TransferParser()           # call class
    count = 0 
    #-----
    """read in transfers file"""
    with open(args.transfer_file,"r") as fh:
        for hit in TP.readTrans(fh):
            lgt_dict.addLGT(hit[TP._UID_1]) # Add uid to dict
            lgt_dict.addLGT(hit[TP._UID_2]) # Add uid to dict
        
    print "Start", datetime.datetime.now()          
    for kmer_dir in kmer_directories:
        kmer_files = glob.glob('%s/*.kmer_counts.csv' % kmer_dir)
        for kmer in kmer_files:
            id = kmer.split("/")[-1].split(".")[0]
            if len(id) < 5: # lgt id
                lgt_tmp_array = []
                lgt_id = id
                if lgt_dict.checkUID(lgt_id):
                    with open(kmer,'r') as lgt_fh:
                        for l in lgt_fh:
                            fields = l.rstrip().split("\t") 
                            if l[0:2] == "ID":
                                pass
                            else:
                                lgt_tmp_array.append([float(i) for i in fields[1:]])
                        lgt_tmer = np.mean(lgt_tmp_array, axis=0)
                        LGT_kmers.addLGTTmer(lgt_id, lgt_tmer)         
            else: # genome file
                if lgt_dict.checkUID(lgt_id):
                    GID = id
                    #print GID
                    LGT_kmers.addLGT(lgt_id,GID)
                    GID_tmp_array = []
                    with open(kmer,'r') as GID_fh:
                        for l in GID_fh:
                            fields = l.rstrip().split("\t")
                            if l[0:2] =="ID":
                                pass
                            else:
                                GID_tmp_array.append([float(i) for i in fields[1:]])
                        GID_tmer = np.mean(GID_tmp_array, axis=0)
                        LGT_kmers.addGenomeTmer(GID, GID_tmer)
            count +=1 
            if count >=1000:
                break
                #pass
        if count >=1000:
                break
                #pass
        
        
                #print id
            #count+=1 # troubleshooting
            #end of for loop
    #printHeader() 
    LGT_kmers.getClosestGID(False)
    #LGT_kmers.getDistHisto()
    LGT_kmers.printInfoHeader()
    LGT_kmers.printInfo()
    #LGT_kmers.printDict()
    print "Stop", datetime.datetime.now()        
         
            
            
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
    parser.add_argument('-t','--transfer_file', help="...")
    parser.add_argument('-k','--kmers_directory', help="...")
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
