#!/usr/bin/env python
###############################################################################
#
# __rank_interactions__.py - rank genomes by its number of interacting genomes
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
    
    # Required files
    # 1. gut_oral_genomes_total_hits_length_clean_contig_removed.csv
    """genome_tree_id_a        img_id_a        bodysite_a      genome_tree_id_b        img_id_b        bodysite_b      hits    length"""
    # 2. img_ids_ordered_by_genome_tree.Mar14.csv
    """genome_tree_id    img_id    order    tax_string"""
    # 3. gut_oral_contigs.csv
    """img_id_a    genome_tree_id_a    contig_a    contig_length_a    start_a    stop_a    length_a    img_id_b    genome_tree_id_b    contig_b    contig_length_b    start_b    stop_b    length_b"""
    
    #variables
    tree_id_genus_taxa = {} # dictionary to store tree_id -> taxa at genus level
    gut_oral_contigs = {} #dictionary to store tree_id -> [img_id, contig_id, contig_length]
     
    # read in img_ids_ordered_by_genome_tree.Mar14.csv
    with open(args.ids_taxa_file,"r") as fh:
        # capture header
        header = fh.readline()
        
        for l in fh:
            tree_id = l.split("\t")[0]
            img_id = l.split("\t")[1]
            tree_order = l.split("\t")[2] 
            taxa_string = l.split("\t")[3]
            #print taxa_string
            if taxa_string.split(";")[0] == "d__Archaea":
                g = taxa_string.split(";")[6] # genus
                f = taxa_string.split(";")[5] # family
                o = taxa_string.split(";")[4] # order
                c = taxa_string.split(";")[3] # class
                p = taxa_string.split(";")[2] # phylum
            else:
                g = taxa_string.split(";")[5] # genus
                f = taxa_string.split(";")[4] # family
                o = taxa_string.split(";")[3] # order
                c = taxa_string.split(";")[2] # class
                p = taxa_string.split(";")[1] # phylum
            # lowest taxa level
            if len(g) == 4:
                if len(f) == 4:
                    if len(o) == 4:
                        if len(c) == 4:
                            if len(p) == 4:
                                pass # removed 4 genomes which have very little taxonomical information
                            else:
                                #print p
                                tree_id_genus_taxa[tree_id] = [p,tree_order] 
                        else:
                            #print c
                            tree_id_genus_taxa[tree_id] = [c,tree_order] 
                    else:
                        #print o
                        tree_id_genus_taxa[tree_id] = [o,tree_order] 
                else:
                    #print f
                    tree_id_genus_taxa[tree_id] = [f,tree_order]   
            else:
                #print g
                tree_id_genus_taxa[tree_id] = [g,tree_order]

    #variables
    counts = {} # dictionary to store counts and img_ids  
    
    with open(args.contigs, "r") as fh:
        # capture header
        header = fh.readline()
        # read through line by line
        for l in fh:
            img_id_a = l.split("\t")[0].rstrip()
            genome_tree_id_a = l.split("\t")[1].rstrip()
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
            # add to dictionary
            if genome_tree_id_a in gut_oral_contigs:
                pass
            else:
                try:
                    gut_oral_contigs[genome_tree_id_a] = [img_id_a,
                                                          contig_a,
                                                          int(contig_length_a)]
                except KeyError:
                    pass
            if genome_tree_id_b in gut_oral_contigs:
                pass
            else:    
                try:
                    gut_oral_contigs[genome_tree_id_b] = [img_id_b,
                                                          contig_b,
                                                          int(contig_length_b)]
                except KeyError:
                    pass
            
    
    with open(args.lgt_events,"r") as fh:
        #capture header
        header = fh.readline()
        for l in fh:
            tree_id_a = l.split("\t")[0]
            img_id_a = l.split("\t")[1]
            tree_id_b = l.split("\t")[3]
            img_id_b = l.split("\t")[4]
            contig_length_cutoff = 10000
            
            # only count interactions between genomes from different genera
            if gut_oral_contigs[tree_id_a][2] > contig_length_cutoff and gut_oral_contigs[tree_id_a][2] > contig_length_cutoff: 
                if "g__" in tree_id_genus_taxa[tree_id_a][0] and "g__" in tree_id_genus_taxa[tree_id_b][0]: #both contain genus level taxa
                    if tree_id_genus_taxa[tree_id_a][0] == tree_id_genus_taxa[tree_id_b][0]:
                        pass
                    else:
                        try:
                            counts[tree_id_a][0] += 1
                        except KeyError:
                            counts[tree_id_a] = [1, img_id_a, tree_id_genus_taxa[tree_id_a][1]]
                        try:
                            counts[tree_id_b][0] += 1
                        except KeyError:
                            counts[tree_id_b] = [1, img_id_b, tree_id_genus_taxa[tree_id_b][1]]
                else:
                    try:
                        counts[tree_id_a][0] += 1
                    except KeyError:
                        counts[tree_id_a] = [1, img_id_a, tree_id_genus_taxa[tree_id_a][1]]
                    try:
                        counts[tree_id_b][0] += 1
                    except KeyError:
                        counts[tree_id_b] = [1, img_id_b, tree_id_genus_taxa[tree_id_b][1]]
    #print header
    print "\t".join(["genome_tree_id",
                     "img_id",
                     "tree_order",
                     "#_lgt_partners"
                     ])
    for key in counts:
        print "\t".join([key,
                         counts[key][1],
                         counts[key][2],
                         str(counts[key][0])
                         ])
                
            

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
    parser.add_argument('ids_taxa_file', help="img_ids_ordered_by_genome_tree.Mar14.csv")
    parser.add_argument('lgt_events', help="gut_oral_genomes_total_hits_length_clean_contig_removed.csv")
    parser.add_argument('contigs', help="gut_oral_contigs.csv")
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
