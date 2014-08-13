#!/usr/bin/env python
###############################################################################
#
# __hitDataStats__.py - Description!
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

# put classes here 

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
    sequencingPlatform = {}
    sequencingCentre   = {}
    phylum             = {}
    genus              = {}
    status             = {}
    bodysite           = {}
    
    with open(args.hitdata, 'r') as fh:
        header = fh.readline()
        for l in fh:
            tabs = l.split("\t")
            bodysite1 = tabs[1]
            genus1    = tabs[2]
            phylum1   = tabs[3]
            seqplat1  = tabs[6]
            seqcent1  = tabs[5]
            status1   = tabs[7]
            bodysite2 = tabs[9]
            genus2    = tabs[10]
            phylum2   = tabs[11]
            seqplat2  = tabs[14]
            seqcent2  = tabs[13]
            status2   = tabs[15]
            hitspre   = int(tabs[16])
            hitspost  = int(tabs[17])
            try:
                sequencingPlatform[seqplat1][0] += hitspre
                sequencingPlatform[seqplat1][1] += hitspost
            except KeyError:
                sequencingPlatform[seqplat1] = [hitspre,hitspost]
                
            try:
                sequencingPlatform[seqplat2][0] += hitspre
                sequencingPlatform[seqplat2][1] += hitspost
            except KeyError:
                sequencingPlatform[seqplat2] = [hitspre,hitspost]
    
    for key in sequencingPlatform.keys():
        print "\t".join([key,
                         str(sequencingPlatform[key][0]),
                         str(sequencingPlatform[key][1])])
            
    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-hit','--hitdata', help="...")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    #parser.add_argument('-X', '--optional_X', action="store_true", type=int,default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
