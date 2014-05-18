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

from multiprocessing import Pool
from subprocess import Popen, PIPE

import os
import errno

import numpy as np
np.seterr(all='raise')

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure


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

def doWork( args ):
    """ Main wrapper"""

    # variables
    col_labels = []
    row_labels = []
    data = []
    
    #Stored arguments as variables
    output_file = args.output_image_name + ".%s" % (args.image_format)
    image_format = args.image_format
    
    
    
    # open the file handle
    with open(args.matrix_file,"r") as fh:
        #capture the first line as header
        col_labels = fh.readline().rstrip().split("\t")[1:] # remove img_ids from the header
        # read through the file line by line
        for line in fh:
            # capture first element of line in row labels
            row_labels.append(line.rstrip().split("\t")[0])
            # collected data in an array of arrays
            # mask if zero
            """ Masked array???
            for i in line.rstrip().split("\t")[1:]:
                if int(i) == 0: 
                    
                else:
            """
            data.append(line.rstrip().split("\t")[1:])
        #-----
        #masked array
        mask = []
        for row in data:
            mask.append([1 if x == str(0) else 0 for x in row])
        
        
        
        #-----
        fig, ax = plt.subplots()
        data = np.array(data).astype('float')
        masked_data = np.ma.masked_array(data,mask)
        """
        #loop over data
        for i in range(len(data)):
            for j in range(len(data)):
                if i > j:
                    #length
                    data[i][j] = (data[i][j])/10000
                    #print data[i][j]
                    
                if i < j:
                    #hits
                    data[i][j] = data[i][j]
                    print data[i][j]
        #qprint np.max(data)    
        """
        #Set mask color to white:
        plt.cm.YlOrRd.set_bad(color='white', alpha=None)
        
        # produce heatmap
        ax.pcolormesh(masked_data, cmap=plt.cm.YlOrRd)    
            
        # put the major ticks at the middle of each cell
        ax.set_xticks(np.arange(data.shape[0])+0.1, minor=False)
        ax.set_yticks(np.arange(data.shape[1])+0.1, minor=False)
        
        # turn off the frame
        ax.set_frame_on(False)
        
        #set the limits of the plot to the limits of the data
        plt.axis([0, len(data), 0, len(data)])

        #want a more natural, table-like display
        ax.invert_yaxis() # flips heatmap
        ax.xaxis.tick_top() # places xlabels on top
        
        # Turn off all the ticks
        jax = plt.gca()

        for t in jax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        
        # colour bar
        cax = ax.imshow(data, aspect = 'auto', interpolation = 'nearest',cmap=plt.cm.YlOrRd)
        cbar = fig.colorbar(cax)
        
        
        
        # plot labels
        ax.set_xticklabels(col_labels, minor=False,size=1,rotation='vertical')
        ax.set_yticklabels(row_labels, minor=False,size=1)
        #plt.savefig(output_file,dpi=args.dpi,format=str(image_format))
        plt.show()

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
    parser.add_argument('-m','--matrix_file', help="Matrix file required")
    parser.add_argument('-o','--output_image_name', default='heatmap', help="Output file name; default: out.txt")
    parser.add_argument('-i','--image_format', default='png', help="Image format; default: png")
    parser.add_argument('-d','--dpi', type=int, default=400, help="Set image dpi")
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
