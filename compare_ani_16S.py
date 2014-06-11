#!/usr/bin/env python
###############################################################################
#
# __compare_ani_16S__.py - Compare ani scores vs 16S distance!
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

class parse16S(object):
    # constants to make code readable
    _img_id_a = 0
    _img_id_b = 1
    _identity = 2
    
    def __init__(self):
        self.prepped = False

    def reset(self):
        self.prepped = False
        
    def read16S(self,fh):
        """Read through a 16S comparison file

        this is a generator function
        """
        line = None
        while True:
            if not self.prepped:
                # strip header
                for l in fh:
                    if l[0:2] == "ID":
                        self.prepped = True
                        break
            for l in fh:
                fields = l.split("\t")
                yield ([fields[0],
                        fields[1],
                        float(fields[2])])
            break # done

class parseANI(object):
    # constants to make the code readable
    "IMG_Taxon1      Species IMG_taxon2      Species2        ANI1    ANI2    AF1     AF2"
    _img_id_a   = 0
    _species_a  = 1
    _img_id_b   = 2
    _species_a  = 3
    _ANI_1      = 4
    _ANI_2      = 5
    _AF_1       = 6
    _AF_2       = 7
    
    def __init__(self):
        self.prepped = False
        
    def readANI(self,fh):
        """Read through a ANI comparison file

        this is a generator function
        """
        line = None
        while True:
            if not self.prepped:
                # strip header
                for l in fh:
                    if l[0:3] == "IMG":
                        self.prepped = True
                        break
            for l in fh:
                fields = l.split("\t")
                yield ([fields[0],
                        fields[1],
                        fields[2],
                        fields[3],
                        float(fields[4]),
                        float(fields[5]),
                        float(fields[6]),
                        float(fields[7])])
            break # done
 
class store16S(object):
    def __init__(self):
       self.dict_16S = {} 
    
    def add16S(self,g1,g2,ID):
        try:
            self.dict_16S[g1][g2] = [ID]
        except KeyError:
            self.dict_16S[g1] = {g2:[ID]}
        try:
            self.dict_16S[g2][g1] = [ID]
        except KeyError:
            self.dict_16S[g2] = {g1:[ID]}
            
    def addANI(self,g1,g2,ani1,ani2):
        try:
            self.dict_16S[g1][g2]+= [ani1,ani2]
        except KeyError:
            #print "\t".join([g1,g2])
            pass
        try:
            self.dict_16S[g2][g1]+= [ani1,ani2]
        except KeyError:
            #print "\t".join([g2,g1])
            pass 
        
    def printOUT(self):
        for g1 in self.dict_16S.keys():
            for g2 in self.dict_16S[g1]:
                print self.dict_16S[g1][g2]
                
                #ID_16S  = str(self.dict_16S[g1][g2][0])
                #ANI_1   = str(self.dict_16S[g1][g2][1])
                #ANI_2   = str(self.dict_16S[g1][g2][2])
                #print "\t".join([g1,g2,ID_16S,ANI_1,ANI2])

class storeANI(object):
    def __init__(self):
        self.dict_ANI = {}
    
    def addANI(self,g1,g2,ani1,ani2):
        try:
            self.dict_ANI[g1][g2] = [ani1,ani2]
        except KeyError:
            self.dict_ANI[g1]= {g2:[ani1,ani2]}
        try:
            self.dict_ANI[g2][g1] = [ani2,ani1]
        except KeyError:
            self.dict_ANI[g2]= {g1:[ani2,ani1]}
    
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

def compare16S_ANI(ani1,ani2,id_16S):
    pass

def printHeader():
    print "\t".join(["img_id_a","img_id_b","16S_ID","ANI_1","ANI_2"])

def printOut(dict_ani,dict_16S):
    for g1 in dict_16S.keys():
        for g2 in dict_16S[g1]:
            print 


def doWork( args ):
    """ Main wrapper"""
    # objects
    S16_p = parse16S()
    ANI_p = parseANI()
    #ANI_db = storeANI()
    S16_db = store16S()
    count = 0
    counter = 0
            
    # read in 16S file
    with open(args.s16_file,"r") as fh_16S:
        for line in S16_p.read16S(fh_16S):
            S16_db.add16S(line[S16_p._img_id_a], line[S16_p._img_id_b], line[S16_p._identity])
            if count >= 10000:
                break                
            count += 1 
        
    # read in ANI file
    with open(args.ani_file,"r") as fh_ani:
        for line in ANI_p.readANI(fh_ani):
            S16_db.addANI(line[ANI_p._img_id_a], line[ANI_p._img_id_b], line[ANI_p._ANI_1], line[ANI_p._ANI_2])
            if counter >= 10000:
                break
            counter += 1
        
    printHeader() # print out header
    S16_db.printOUT() # print out data
    
     
                
            
            
            
            
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
    parser.add_argument('-ani_file','--ani_file', help="Please provide file containing ANI scores")
    parser.add_argument('-s16_file','--s16_file', help="Please provide file containing 16S scores")
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
