#!/usr/bin/env python
###############################################################################
#
# compute_whole_genome_nucmer.py - run nucmer between IMG genomes
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
import datetime

from multiprocessing import Pool
from subprocess import Popen, PIPE

import os
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

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE, stderr=PIPE)
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""

    subgrp = 1000

    # make sure the out dir exists
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)    
    
    """Read through tab-delimited file containing pairwise comparisons to be run"""
    pool = Pool(args.num_threads)
    count = 0 # counter for checkpoint of runtime
    cmds = [[]]
    stdouts = []

    with open(args.genome_comparisons) as fh: 
        """ ID_1         ID_2         %IDENTITY
            2500069000    2501799900    93.58
            2500069000    637000103    83.25
            2500069000    637000104    81.01
        """
        header = fh.readline()
        counter = 0 
        for l in fh:
            fields = l.rstrip().split("\t")
            #cmds[-1].append("nucmer %s.fna %s.fna --mum --coords -p %s.gen > /%s" % (fields[0], fields[1], "%sv%s" %(fields[0],fields[1]), args.out_dir))
            fasta1 = os.path.join(args.fasta_dir, "%s.%s" % (fields[0], args.fasta_suffix))
            fasta2 = os.path.join(args.fasta_dir, "%s.%s" % (fields[1], args.fasta_suffix))
            nucmer_prefix = os.path.join(args.out_dir, "%sv%s.gen" %(fields[0],fields[1]))
            
            cmds[-1].append("nucmer %s %s --mum --coords -p %s" % (fasta1,
                                                                   fasta2,
                                                                   nucmer_prefix)
                            )
            #cmds[-1].append("echo %s" %("%s %s %s" %(fasta1, fasta2, nucmer_prefix)))
            counter += 1
            if counter >= subgrp:
                cmds.append([])
                counter = 0

    print "start", datetime.datetime.now()
    for sub_cmds in cmds:
        stdouts.append(pool.map(runCommand, sub_cmds))            # list of tuples [(stdout, stderr)]
        print "%d done" % subgrp, datetime.datetime.now()
    
    print "finish", datetime.datetime.now()
    
    print "writing stdouts"
    for (out, err) in stdouts[0]:
        err_file = "%s.txt" % err.split('/')[2].split('.')[0]
        with open(os.path.join(args.out_dir, err_file), 'w') as err_fh:
            for line in err:
                err_fh.write(line)
    
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

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_comparisons', help="File containing list of fasta files to be compared")
    parser.add_argument('out_dir', help="Please provide output directory")
    parser.add_argument('fasta_dir', help="Location of fasta files")
    parser.add_argument('--fasta_suffix', default='fna', help="Prefix of fasta files")
    parser.add_argument('--num_threads', type=int, default=3, help="Number of threads to use")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
