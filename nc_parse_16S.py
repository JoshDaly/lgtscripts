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
import glob

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure

class JoshError(Exception): pass

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

def doWork( args ):
    """ Main wrapper"""
    NP = NucMerParser()
        
    print "ID_1", "\t", "ID_2", "\t" ,"%IDENTITY"
    #coords_files = args.coords
    #print len(args.coords)
    
    listing = glob.glob('%s/*.coords' % args.coordsdir)
    
    for c_file in listing:
        with open(c_file, "r") as fh:
            for hit in NP.readNuc(fh):
                if hit[NP._IDENTITY] < args.identity_cutoff:
                    print hit[NP._ID_1], "\t", hit[NP._ID_2], "\t", hit[NP._IDENTITY]
            # required for NucMer parser. Needs to reset otherwise it doesn't know there is a new file
            NP.reset()

    """# parse fasta file
    CP = ContigParser()
    contigs = {}
    for c_file in args.contigs:
        try:
            with open(c_file, "r") as fh:
                for header, seq in CP.readFasta(fh):
                    contigs[header] = seq
        except:
            print "Error opening file:", c_file, exc_info()[0]
            raise
    
    #print contigs.keys()

    # parse nucmer file
    NP = NucMerParser()
    with open(args.coords, 'r') as fh:
        for hit in NP.readNuc(fh):
	    if hit[NP._IDENTITY] >= 99:
	        if hit[NP._LEN_1] >= 500 and hit[NP._LEN_2] >= 500:
                    print "===="
                    print contigs[hit[NP._ID_1]][hit[NP._START_1]:hit[NP._END_1]]
                    print "--"
                    print contigs[hit[NP._ID_2]][hit[NP._START_2]:hit[NP._END_2]]
                    print hit[NP._START_1],"\t",hit[NP._END_1],"\t",hit[NP._START_2],"\t",hit[NP._END_2],"\t",hit[NP._LEN_1],"\t",hit[NP._LEN_2],"\t",hit[NP._ID_1],"\t",hit[NP._ID_2]

    test = contigs[contigs.keys()[0]]
    #print test
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
    parser.add_argument('-c','--coordsdir', help="directory containing nucmer coords files")
    parser.add_argument('-i','--identity_cutoff',type=int,default=97,help="Set 16S identity cutoff default: 97%")
    #parser.add_argument('contigs', nargs='+', help="contigs mapped using nucmer")
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

