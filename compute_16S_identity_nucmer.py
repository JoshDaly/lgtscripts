#!/usr/bin/env python
###############################################################################
#
# compute_16S_identity_nucmer.py - compute 16S identity between IMG genomes
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
import glob 

from multiprocessing import Pool
from subprocess import Popen, PIPE, call,STDOUT

import os
import errno

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
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()
def returnGenomeName(fasta_file):
    """../../img_4.1_all_lgt/16S_fasta_files/651285000/651285000_16S.fna"""
    slashes = fasta_file.rstrip().split("/")
    genome_name = slashes[-1].split(".")[0]
    return genome_name

def doesDirectoryExist(output_dir):
    if not os.path.exists(output_dir):
        os.system("mkdir %s" % (output_dir))

def runJobs((output_directory,fasta_1,fasta_2,genome_1,genome_2)):
    doesDirectoryExist(output_directory)
    os.chdir(output_directory)
    with open("/dev/null","w") as fh:
        call("nucmer %s %s --mum --coords -p %s" % (fasta_1, fasta_2, "%s_v_%s" %(genome_1,genome_2)),stdout=fh,stderr=STDOUT)

def doWork( args ):
    """ Main wrapper"""
    
    # run something external in threads
    counter = 0   # set up a counter
    pool = Pool(args.num_threads, maxtasksperchild=1)                              # 6 threads
    cmds = []
    fasta_files = glob.glob('%s/*/*.fna' % args.fasta_directory)
    jobs = []
    stdouts = []
    subgrp = 1000
    count = 0

    # set stdout to file
    #sys.stderr = open('nucmer.log',"w")
    
    for i in range(len(fasta_files)-1):                    # Mikes example commands for running the script
        for j in range(i+1, len(fasta_files)): # +1 and -1 to the for loops, means that only the bottom half of the triangle will be compared.
            genome_1 = returnGenomeName(fasta_files[i])
            genome_2 = returnGenomeName(fasta_files[j])
            output_directory = "/srv/projects/lgt/img_4.1_all_lgt/16S_nucmer_97/%s_v_%s" % (genome_1,genome_2)              
            fasta_1 = "../../16S_fasta_files/%s/%s.fna" % (genome_1.split("_")[0],genome_1)
            fasta_2 = "../../16S_fasta_files/%s/%s.fna" % (genome_2.split("_")[0],genome_2)
            jobs.append((output_directory,fasta_1,fasta_2,genome_1,genome_2))
            counter += 1
            count +=1
            if counter >= subgrp:
                jobs.append(())
                counter = 0
            if count >=10:
                break
        if count >=10:
                break
    
    print "Start", datetime.datetime.now()
    pool.map(runJobs,jobs)
    print "%d done" % subgrp, datetime.datetime.now()
    print "finish", datetime.datetime.now()
    
    #run commands
    #pool.map(runJobs,jobs)
    #pool.close()
    #pool.join()
    
    
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
    parser.add_argument('-f','--fasta_directory', help="Directory containing 16S fasta file directories")
    parser.add_argument('-so','--std_out_dir', help="Standard out directory")
    parser.add_argument('-num_threads','--num_threads', type=int, default=6, help="Set the number of threads for multiplexing")
    #parser.add_argument('positional_arg', help="")
    #parser.add_argument('positional_arg2', nargs='+' help="")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
