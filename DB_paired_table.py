#!/usr/bin/env python
###############################################################################
#
# __DB_paired_table__.py - DATABASE FILE: create paired genome table for database!
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
        
class METAparser(object):
    def __init__(self,l):
        self.readMetadata(l)
        
    def readMetadata(self,l):
        """ taxon  = 0, genome_name = 4"""
        """/srv/db/img/07042014/metadata/img_metadata_07042014.tsv"""
        tabs = l.rstrip().split("\t")
        self._img_id        = tabs[0]
        self._genome_name   = tabs[4]
        self._body_site     = tabs[56]
        
class ANIparser(object):
    def __init__(self,l):
        self.readANI(l)
        
    def readANI(self,l):
        """IMG_Taxon1      Species IMG_taxon2      Species2        ANI1    ANI2    AF1     AF2"""
        tabs = l.rstrip().split("\t")
        self._img_id_1  = tabs[0]
        self._img_id_2  = tabs[2]
        self._ANI_1     = tabs[4]
        self._ANI_2     = tabs[5]
        
class genomeTreeParser(object):
    def __init__(self,l):
        self.readGT(l)
    
    def readGT(self,l):
        """genome_tree_id  img_id  order   tax_string"""
        tabs = l.rstrip().split("\t")
        self._gt_id  = tabs[0]  
        self._img_id = tabs[1]

class paired_data(object):
    def __init__(self):
        self.img_to_gt_dict = {} # dict to store img-> genome tree ids
        self.img_metadata_dict = {}
        self.ANI_scores = {}
        self.path_to_genome = {} # genome_tree_id -> PATH_TO_FILE
        self.body_site_dict = {} # genome_tree_id -> body_site
         
    def addGT(self,l):
        line = genomeTreeParser(l)
        self.img_to_gt_dict[line._img_id] = line._gt_id
        
    def printGTs(self):
        for key in self.img_to_gt_dict.keys():
            print "\t".join([key,self.img_to_gt_dict[key]])
        
    def addMETA(self,l):
        line = METAparser(l)
        # Simplify IMG bodysite metadata
        # Simplification not needed: Nose, Oral, Airways, Ear, Eye, Gastrointestinal tract, Urogenital tract
        """batches"""
        # 1 = gastrointestinal tract
        # 2 = oral
        # 3 = airways
        # 4 = urogenital tract, nose, ear, eye, nose
        # 5 = skin 
        # 6 = internal organs
        # 7 = plant
        # 8 = other
        
        body_site = "NA"
        # other 11954 genomes 
        batch = 8
        if line._body_site == "Nose": # 30 genomes
            body_site = "nose"
            batch = 4
        if line._body_site == "Oral": # 543 genomes
            body_site = "oral"
            batch = 2
        if line._body_site == "Airways": # 352 genomes
            body_site = "airways"
            batch = 3
        if line._body_site == "Ear": # 15 genomes 
            body_site = "ear"
            batch = 4
        if line._body_site == "Eye": # 10 genomes
            body_site = "eye"
            batch = 4
        if line._body_site == "Gastrointestinal tract": # 1034
            body_site = "gastrointestinal tract"
            batch = 1
        if line._body_site == "Urogenital tract": # 256 
            body_site = "urogenital tract"
            batch = 4
        # Skin # 221
        skin = ["skin","abdomen","ankle","limb","wound"]
        for skin_site in skin:
            if skin_site in line._body_site.lower():
                body_site = "skin"
                batch  = 5
                break
        # Internal organs # 257 genomes
        internal_organs = ["spinal cord","heart","liver","lymph nodes","blood","bladder","bone","brain"]
        for organ in internal_organs:
            if organ in line._body_site.lower():
                body_site = "internal organ"
                batch  = 6
                break
        # Plant # 176 genomes
        plant = ["plant","root"]
        for plant_site in plant:
            if plant_site in line._body_site.lower():
                body_site = "plant"
                batch  = 7
                break
        # add to dictionary
        self.img_metadata_dict[line._img_id] = [line._genome_name,body_site,batch]
    
    def printMETA(self):
        for key in self.img_metadata_dict.keys():
            print "\t".join([key,self.img_metadata_dict[key][0],self.img_metadata_dict[key][1]])
        
    def addANI(self,l,pid):
        line = ANIparser(l)
        self.ANI_scores[pid] = [line._img_id_1,line._img_id_2,line._ANI_1,line._ANI_2] 
        
    def printANIs(self):
        for key in self.ANI_scores.keys():
            print "\t".join([str(key),self.ANI_scores[key][0],self.ANI_scores[key][1],self.ANI_scores[key][2],self.ANI_scores[key][3]])
    
    def checkID(self,l,type): # check if it has genome tree id
        if type == "ANI":
            line =  ANIparser(l)
            if line._img_id_1 in self.img_to_gt_dict and line._img_id_2 in self.img_to_gt_dict:
                return True
        if type == "META":
            line = METAparser(l)
            if line._img_id in self.img_to_gt_dict:
                return True
            
    def BodySiteDict(self):
        for pid in  self.ANI_scores.keys():
            try:
                img_id_1            = self.ANI_scores[pid][0]
                img_id_2            = self.ANI_scores[pid][1]
                ANI_1               = str(self.ANI_scores[pid][2])
                ANI_2               = str(self.ANI_scores[pid][3])
                genome_tree_id_1    = self.img_to_gt_dict[img_id_1]
                genome_tree_id_2    = self.img_to_gt_dict[img_id_2]
                batch_a             = self.img_metadata_dict[img_id_1][2]
                batch_b             = self.img_metadata_dict[img_id_2][2]
                body_site_a         = self.img_metadata_dict[img_id_1][1]
                body_site_b         = self.img_metadata_dict[img_id_2][1]
                # add to dictoinary
                self.body_site_dict[genome_tree_id_1] = body_site_a
                self.body_site_dict[genome_tree_id_2] = body_site_b
            except KeyError:
                pass
            
    def printBodySiteDict(self):
        for id in self.body_site_dict.keys():
            print "\t".join([id,self.body_site_dict[id]])
        
            
    def printPIDtable(self):
        for pid in  self.ANI_scores.keys():
            try:
                img_id_1            = self.ANI_scores[pid][0]
                img_id_2            = self.ANI_scores[pid][1]
                ANI_1               = str(self.ANI_scores[pid][2])
                ANI_2               = str(self.ANI_scores[pid][3])
                genome_tree_id_1    = self.img_to_gt_dict[img_id_1]
                genome_tree_id_2    = self.img_to_gt_dict[img_id_2]
                batch_a             = self.img_metadata_dict[img_id_1][2]
                batch_b             = self.img_metadata_dict[img_id_2][2]
                body_site_a         = self.img_metadata_dict[img_id_1][1]
                body_site_b         = self.img_metadata_dict[img_id_2][1]
                batch = 0
                if batch_a == batch_b:
                    batch = batch_a
                    print "\t".join([str(pid),genome_tree_id_1,genome_tree_id_2,ANI_1,ANI_2,str(batch)])
                if batch_a > batch_b:
                    batch = batch_a
                    print "\t".join([str(pid),genome_tree_id_1,genome_tree_id_2,ANI_1,ANI_2,str(batch)])
                if batch_a < batch_b:
                    batch = batch_b
                    print "\t".join([str(pid),genome_tree_id_1,genome_tree_id_2,ANI_1,ANI_2,str(batch)])
            except KeyError:
                pass
            
    def addPathToGenome(self):
        """/srv/db/img/07042014/genomes/%s/%s.fna"""
        for pid in self.ANI_scores.keys():
            try:
                img_id_1                        = self.ANI_scores[pid][0]
                img_id_2                        = self.ANI_scores[pid][1]
                path_to_file_1                  = "/srv/db/img/07042014/genomes/%s/%s.fna" % (self.ANI_scores[pid][0],self.ANI_scores[pid][0])
                path_to_file_2                  = "/srv/db/img/07042014/genomes/%s/%s.fna" % (self.ANI_scores[pid][1],self.ANI_scores[pid][1])
                self.path_to_genome[img_id_1]   = path_to_file_1
                self.path_to_genome[img_id_2]   = path_to_file_2
            except KeyError:
                pass
        
    def printPathToGenomes(self):
        print "\t".join(["genome_tree_id","PATH_TO_FILE"])
        for id in self.path_to_genome.keys():
            print "\t".join([self.img_to_gt_dict[id], self.path_to_genome[id]])
                 
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

def printHEADER():
    print "\t".join(["pid","genome_1","genome_2","ANI_1","ANI_2","batch"])

def doWork( args ):
    """ Main wrapper"""
    # objects
    PD = paired_data() # call class
    pid = 1 # paired genome ID
    TS = 0
    
    # read in genome tree file
    with open(args.genome_tree_file,"r") as fh:
        header = fh.readline()
        for l in fh:
            PD.addGT(l)
    
    
    # read in ANI file
    with open(args.ANI_file,"r") as fh:
        header = fh.readline()
        for l in fh:
            if PD.checkID(l, "ANI"): # ID in genome tree list
                PD.addANI(l, pid)
                pid += 1
            TS += 1 
            #if TS >= 50000:
            #    break
        
    # read in IMG metadata
    with open(args.img_metadata,"r") as fh:
        header = fh.readline()
        for l in fh:
            if PD.checkID(l, "META"): # ID in genome tree list
                PD.addMETA(l)
    #printHEADER()
    #PD.printPIDtable()            
    #PD.addPathToGenome()
    #PD.printPathToGenomes()    
    PD.BodySiteDict()
    PD.printBodySiteDict()        
            
            
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
    parser.add_argument('-gt','--genome_tree_file', help="...")
    parser.add_argument('-ani','--ANI_file', help="...")
    parser.add_argument('-m','--img_metadata', help="...")
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
