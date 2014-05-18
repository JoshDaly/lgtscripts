#!/usr/bin/env python
###############################################################################
#
# nc_parse.py - parse nucmer
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import argparse
import sys

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

def doWork( args ):
    """ Main wrapper"""
    #objects
    contam_ids= {}
    
    """
    #-----
    # parse nucmer file
    NP = NucMerParser()
    with open(args.coord, 'r') as fh:
        for hit in NP.readNuc(fh):
            lines=hit[NP._ID_2].split("-")
            contig=lines[0].rstrip()
            id_a=lines[1].split(":")[1].rstrip()
            genome_tree_a=lines[2].split(":")[1].rstrip()
            start=int(lines[3].split(":")[1].rstrip())
            stop=int(lines[4].split(":")[1].rstrip())
            id_b=lines[5].split(":")[1].rstrip()
            genome_tree_b=lines[6].split(":")[1].rstrip()
            unique_id=lines[7].rstrip()
            contam_ids[unique_id] = [contig,id_a,genome_tree_a,id_b,genome_tree_b,start,stop]
    #for key in contam_ids:
    #    print key
    """
    #-----
    #read in dirty uids lists
    #univec
    with open(args.univec_list, 'r') as fh:
        for l in fh:
            id = l.rstrip()
            contam_ids[id]= 1
    #blast
    with open(args.blast_list, "r") as fh:
        for l in fh:
            id = l.rstrip()
            contam_ids[id]=0    
    #print contam_ids
            
    # read in gut_oral_contigs_uniq.fna   
    with open(args.fasta_file,"r") as fh:
        for l in fh:
            # check if fasta header
            if ">" in l:
                fasta_info = l.split(">")[1].split("-")
                contig= fasta_info[0].rstrip()
                img_a= fasta_info[1].split(":")[1].rstrip()
                genome_tree_a= fasta_info[2].split(":")[1].rstrip()
                start= fasta_info[3].split(":")[1].rstrip()
                stop= fasta_info[4].split(":")[1].rstrip()
                img_b= fasta_info[5].split(":")[1].rstrip()
                genome_tree_b= fasta_info[6].split(":")[1].rstrip()
                uid= fasta_info[7].rstrip()
                try:
                    if uid in contam_ids:
                        contam_ids[uid] = [contig,img_a,genome_tree_a,img_b,genome_tree_b,start,stop]
                except KeyError:
                    pass
        #print contam_ids
                
    
    
    #-----
    """img_id_a        genome_tree_id_a        contig_a        contig_length_a start_a stop_a  length_a        img_id_b        genome_tree_id_b        contig_b        contig_length_b start_b stop_b  length_b"""
    
    #parse gut_oral_contigs.csv file
    with open(args.gut_oral_contigs,"r") as fh:
        #header
        header=fh.readline().rstrip()
        #read through file line by line
        print header
        for l in fh:
            Good_to_print = True
            tabs= l.split("\t")
            img_id_a= tabs[0]
            genome_tree_id_a= tabs[1]
            contig_a= tabs[2]
            contig_length_a= tabs[3]
            start_a= tabs[4]
            stop_a= tabs[5]
            length_a= tabs[6]
            img_id_b= tabs[7]
            genome_tree_id_b= tabs[8]
            contig_b= tabs[9]
            contig_length_b= tabs[10]
            start_b= tabs[11]
            stop_b= tabs[12]
            length_b= tabs[13].rstrip()
            for uid in contam_ids:
                if contig_a == contam_ids[uid][0] or contig_b == contam_ids[uid][0]:
                    if int(start_a) == int(contam_ids[uid][5]) and int(stop_a) == int(contam_ids[uid][6]):
                        Good_to_print = False
                    if int(start_b) == int(contam_ids[uid][5]) and int(stop_b) == int(contam_ids[uid][6]):
                        Good_to_print = False
            if Good_to_print:
                print l.rstrip() 
           
    
    
    
    
    #print len(contam_ids.keys())
                     
                   
                   
    #print test

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

def printSeq(rawSeq, start, stop):
    # if start > stop
    pass

def revComp(seq):
    # see rob edwards
    pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--coord', help="nucmer coords file")
    parser.add_argument('-u','--univec_list', help="List of dirty sequence IDs (UniVec)")
    parser.add_argument('-f','--fasta_file', help="Fasta file containing transfer and uid information")
    parser.add_argument('-b','--blast_list', help="List of dirty sequence IDs (blast)")
    parser.add_argument('-go','--gut_oral_contigs', help="gut_oral_contigs.csv")
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

