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
import glob

from multiprocessing import Pool
from subprocess import Popen, PIPE

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

class NucMerParser:
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

  
class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readFasta(self, fp): # this is a generator function
        header = None
        seq = None
        while True:
            for l in fp:
                if l[0] == '>': # fasta header line
                    if header is not None:
                        # we have reached a new sequence
                        yield header, "".join(seq)
                    header = l.rstrip()[1:].partition(" ")[0] # save the header we just saw
                    seq = []
                else:
                    seq.append(l.rstrip())
            # anything left in the barrel?
            if header is not None:
                yield header, "".join(seq)
            break  
  
  
  

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
    
    """Read in file containing contig lengths"""
    """img_id_a        genome_tree_id_a        contig_a        contig_length_a start_a stop_a  length_a        img_id_b        genome_tree_id_b        contig_b        contig_length_b start_b stop_b  length_b"""
    #variables
    contigs_dict = {}
    #contig_size_cutoff = args.contig_size_cutoff
    dirty_contigs = {} #dict to store contaminated contigs
    
    #read in dirty contigs and create dict 
    with open(args.dirty_contigs, "r") as fh:
        for l in fh:
            dirty_contigs[l.strip()] = 1
    #print dirty_contigs
    
    with open(args.contigs_file,"r") as fh:
        #capture header
        header = fh.readline()
        #read through line by line
        for l in fh:
            img_id_a = l.split("\t")[0].rstrip()
            genome_tree_id_a =  l.split("\t")[1].rstrip()
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
            
            try:
                contigs_dict[img_id_a][contig_a] = int(contig_length_a) 
            except KeyError:
                contigs_dict[img_id_a] = {contig_a:int(contig_length_a)}
            try:
                contigs_dict[img_id_b][contig_b] = int(contig_length_b) 
            except KeyError:
                contigs_dict[img_id_b] = {contig_b:int(contig_length_b)}


    #-----
    #objects
    transfer_annotations = {}
    
    #read in annotation file
    with open(args.anno_file,"r") as fh:
        # no header 
        for l in fh:
            tabs=l.split("\t")
            lines=tabs[0].split("|")
            contig=lines[0].rstrip()
            id_a=lines[1].split(":")[1].rstrip()
            genome_tree_a=lines[2].split(":")[1].rstrip()
            start=int(lines[3].split(":")[1].rstrip())
            stop=int(lines[4].split(":")[1].rstrip())
            id_b=lines[5].split(":")[1].rstrip()
            genome_tree_b=lines[6].split(":")[1].rstrip()
            unique_id=lines[7].rstrip()
            annotation=tabs[8]
            COG_group=tabs[5]
            COG_annotation=tabs[6]
            try:
                transfer_annotations[id_a][id_b][contig] +=  [[genome_tree_a,genome_tree_b,unique_id,start,stop,COG_group,COG_annotation,annotation]]
            except KeyError:
                try:
                    transfer_annotations[id_a][id_b][contig] = [[genome_tree_a,genome_tree_b,unique_id,start,stop,COG_group,COG_annotation,annotation]]
                except KeyError:
                    try: 
                        transfer_annotations[id_a][id_b] = {contig:[[genome_tree_a,genome_tree_b,unique_id,start,stop,COG_group,COG_annotation,annotation]]}
                    except KeyError:
                        transfer_annotations[id_a] = {id_b:{contig:[[genome_tree_a,genome_tree_b,unique_id,start,stop,COG_group,COG_annotation,annotation]]}}
        
    
    """
    #-----
    # global variables
    parsing = False
    transfer_annotations = {}

    # read in file containing annotations
    with open(args.anno_file,"r") as fh:
        # no header
        header_data= {}
        for l in fh:
            if "###" in l:
                parsing = True
                #parsing = not parsing
                #continue.
                tabs = l.split("\t")
                #print "\t".join([tabs[0],tabs[1],tabs[2],tabs[3],tabs[4],tabs[5],tabs[6]])
                id_a = tabs[1].rstrip()
                contig = tabs[2].rstrip() 
                start = tabs[3].rstrip()
                stop = tabs[4].rstrip()
                length = tabs[5].rstrip()
                id_b =  tabs[6].rstrip()
                header_data['id_a'] = id_a
                header_data['id_b'] = id_b
                header_data['contig'] = contig
                
            if parsing and "###" not in l and "@@@" not in l:
                tabs_2 = l.split("\t")
                start = tabs_2[0].rstrip()
                stop = tabs_2[1].rstrip()
                annotation = tabs_2[2].rstrip()
                try:
                    transfer_annotations[header_data['id_a']][header_data['id_b']][header_data['contig']] +=  [[start,stop,annotation]]
                except KeyError:
                    try:
                        transfer_annotations[header_data['id_a']][header_data['id_b']] = {header_data['contig']:[[start,stop,annotation]]}
                    except KeyError:
                        try: 
                            transfer_annotations[header_data['id_a']] += {header_data['id_b']:{header_data['contig']:[[start,stop,annotation]]}}
                        except KeyError:
                            transfer_annotations[header_data['id_a']] = {header_data['id_b']:{header_data['contig']:[[start,stop,annotation]]}}        
            if "@@@" in l:
                parsing = False
                header_data = {}
    """
    #-----
    
    #global variables
    ids_dict = {} #dictinary to store img_id -> genome_tree_id
    
    #read in genome tree file
    """genome_tree_id    img_id"""
    with open(args.id_file,"r") as fh:
        #no header
        #read file line by line
        for l in fh:
            genome_tree_id = l.split(" ")[0].rstrip() 
            img_id = l.split(" ")[1].rstrip()
            ids_dict[img_id] = genome_tree_id
    #-----
    #variables
    gut_ids = {}
    oral_ids = {}
    both_ids = {}   
    
    with open(args.gut_ids,"r") as fh:
        for l in fh:
            id = l.rstrip()
            gut_ids[id] = 0
            
    with open(args.oral_ids,"r") as fh:
        for l in fh:  
            id = l.rstrip()      
            oral_ids[id] = 0
            
    with open(args.both_ids,"r") as fh:
        for l in fh:
            id = l.rstrip() 
            both_ids[id] = 0
    
    

    #-----
    NP = NucMerParser()
    
    #Need to include glob command to parse MANY coords files
    listing = glob.glob('%s/*.coords' % args.coords_dir)
    # tracking variables
    hits = 0 # total number of hits between the genome
    length = 0 # total length of DNA shared between genomes
    interaction = 0 # variable to indicate whether it has data
       
    #print file header
    print "\t".join(["genome_tree_id_a",
                     "img_id_a",
                     "body_site_a",
                     "genome_tree_id_b",
                     "img_id_b",
                     "body_site_b",
                     "hits", 
                     "length",
                     "interaction",
                     "phage",
                     "plasmid",
                     "transposon"
                     ])   
    
    """    Script for calculating the number of hits between genomes
           and the total length of shared genes*
    """
    for c_file in listing: 
        # split the img IDs from the file name in the form ../gut_oral_all_v_all/2500069000v637000103.gen.coords
        img_id_a = c_file.split('/')[2].split('.')[0].split('v')[0]
        img_id_b = c_file.split('/')[2].split('.')[0].split('v')[1]
        bodysite_a = "n"
        bodysite_b = "n"
        
        if img_id_a in gut_ids:
            body_site_a = "gut"
        elif img_id_a in oral_ids:
            body_site_a = "oral"
        elif img_id_a in both_ids:
            body_site_a = "both"
        if img_id_a in gut_ids:
            body_site_b = "gut"
        elif img_id_a in oral_ids:
            body_site_b = "oral"
        elif img_id_a in both_ids:       
            body_site_b = "both"
        
        # open the file as fh
        with open(c_file, "r") as fh:
            #set interaction to no
            interaction = "n"
            phage = "n"
            plasmid = "n"
            transposon = "n"
            
            # apply NucmerParser (NP) function readNuc on fh
            for hit in NP.readNuc(fh):
                # does it satisfy the criteria?
                try:
                    if (hit[NP._IDENTITY] >= 99 
                        and hit[NP._LEN_1] >= 500 
                        and hit[NP._LEN_2] >= 500
                        and img_id_a in ids_dict
                        and img_id_b in ids_dict
                        ):
                            interaction = "y"
                            if (hit[NP._ID_1] not in dirty_contigs
                                and hit[NP._ID_2] not in dirty_contigs
                                ):
                                # add one to the count
                                # print hit[NP._LEN_1],"\t",hit[NP._LEN_2] 
                                
                                try:
                                    for contig in transfer_annotations[img_id_a][img_id_b]:
                                        for i in transfer_annotations[img_id_a][img_id_b][contig]:
                                            if "phage" in i[-1].lower():
                                                phage = "y"
                                            if "plasmid" in i[-1].lower():
                                                plasmid = "y"
                                            if "transposon" in i[-1].lower():
                                                transposon = "y"
                                            
                                except KeyError:
                                    pass
                                try:
                                    for contig in transfer_annotations[img_id_b][img_id_a]:
                                        for i in transfer_annotations[img_id_b][img_id_a][contig]:
                                            if "phage" in i[-1].lower():
                                                phage="y"
                                            if "plasmid" in i[-1].lower():
                                                plasmid = "y"
                                            if "transposon" in i[-1].lower():
                                                transposon = "y"
                                except KeyError:
                                    pass
                                                                    
                                hits += 1
                                length = ((hit[NP._LEN_1] + hit[NP._LEN_2]) / 2) + length
                                #print "Yar, she be dirty"
                            #else:
                                #print "Yep, theys be dirty!" + [hit[NP._ID_1]] + "\t" + [hit[NP._ID_2]]                           
                except KeyError:
                    pass
            # Only print if it has at least one hit
            
            if interaction == "y":                
                print "\t".join([ids_dict[img_id_a],
                                 img_id_a,
                                 body_site_a,
                                 ids_dict[img_id_b],
                                 img_id_b,
                                 body_site_b,
                                 str(hits),
                                 str(length),
                                 interaction,
                                 phage,
                                 plasmid,
                                 transposon
                                 ])
            
            # reset hits and length to zero
            hits = 0
            length = 0
            # Nucmer parser needs to be reset to read in new nucmer coords file
            NP.reset()
    
    

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
    parser.add_argument('-cd','--coords_dir', help="directory containing nucmer coords files")
    parser.add_argument('-cf','--contigs_file', help="file containing contig length information")
    parser.add_argument('-d','--dirty_contigs', help="file containing dirty contig IDs")
    parser.add_argument('-g','--anno_file',help="File containing annotation information")
    parser.add_argument('-id','--id_file',help="File containing matching genome tree and img IDs")
    parser.add_argument('-gi','--gut_ids', help="gut_img_ids")
    parser.add_argument('-oi','--oral_ids', help="oral_img_ids")
    parser.add_argument('-bi','--both_ids', help="ids_present_gut_and_oral.csv")
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
