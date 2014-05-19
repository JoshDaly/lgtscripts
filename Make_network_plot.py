#!/usr/bin/env python
###############################################################################
#
# __make_network_plot__.py - Produce network plot from LGT interaction file
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

import networkx as nx 

#import os
#import errno

import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
import matplotlib.pyplot as plt
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
    #objects
    ids_dict = {}
    unordered_ids_list = []
    working_ids_list = []
    ids_list = []
    genome_tree_img_dict = {}
    genome_interactions = []
    phylum_cols = {}
    
    #-----
    #variables to store phyla lists
    intra_phyla = 0 
    hits_accumulative = 0
    phylums = {}
    
    
    with open(args.phyla,"r") as fh:
        for l in fh:
            tabs= l.split()
            id= tabs[0].rstrip()
            phylum= tabs[1].rstrip()
            if phylum=="firmicutes":
                phylums[id] = 0
            elif phylum=="actinobacteria":
                phylums[id] = 1
            elif phylum=="synergistes":
                phylums[id] = 2
            elif phylum=="proteobacteria":
                phylums[id] = 3
            elif phylum=="fusobacteria":
                phylums[id] = 4
            elif phylum=="bacteroidetes":
                phylums[id] = 5
    
    #-----
    with open(args.genome_tree_file,"r") as fh:
           for l in fh:
               genome_tree_id = l.split('\t')[0].rstrip()
               img_id = l.split('\t')[1].rstrip()
               ladder_position = l.split('\t')[2].rstrip()
               taxa_string = l.split('\t')[3].rstrip()
               genome_tree_img_dict[img_id] = [ladder_position,genome_tree_id,taxa_string] 
               
    #-----           
    """"""
    # read in file containing annotation data
    #objects
    transfer_annotations = {}
    
    #read in annotation file
    with open(args.annotation_file,"r") as fh:
        # no header 
        for l in fh:
            tabs = l.split("\t")
            lines = tabs[0].split("-")
            contig = lines[0].rstrip()
            id_a = lines[1].split(":")[1].rstrip()
            genome_tree_a = lines[2].split(":")[1].rstrip()
            start = int(lines[3].split(":")[1].rstrip())
            stop = int(lines[4].split(":")[1].rstrip())
            id_b = lines[5].split(":")[1].rstrip()
            genome_tree_b = lines[6].split(":")[1].rstrip()
            unique_id = lines[7].rstrip()
            annotation = tabs[8]
            try:
                transfer_annotations[id_a][id_b][contig] +=  [[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation]]
            except KeyError:
                try:
                    transfer_annotations[id_a][id_b][contig] = [[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation]]
                except KeyError:
                    try: 
                        transfer_annotations[id_a][id_b] = {contig:[[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation]]}
                    except KeyError:
                        transfer_annotations[id_a] = {id_b:{contig:[[genome_tree_a,genome_tree_b,unique_id,start,stop,annotation]]}} 
    
    #-----
    # To print out the annotations
    """
    for id_a in transfer_annotations:
        for id_b in transfer_annotations[id_a]:
            for contig in transfer_annotations[id_a][id_b]:
                for i in transfer_annotations[id_a][id_b][contig]:
                    if "phage" in i[-1].lower():
                        print id_a+ "\t" + id_b + "\t" + contig + "\t" + i[-1] + "\t"+i[2]
    
    """
    #-----
    """Opening the hits-length file"""
    with open(args.input_file, "r") as fh:
          # capture the file header
          header = fh.readline()
                         
          #read through the file line by line
          for l in fh:
              # assign each column in line to a variable
              tabs = l.split("\t")
              genome_tree_a = tabs[0]
              id_a = tabs[1]
              body_site_a = tabs[2]
              genome_tree_b = tabs[3]
              id_b = tabs[4]
              body_site_b = tabs[5]
              hits = int(tabs[6].rstrip())
              length = int(tabs[7].rstrip())
              #interaction = tabs[8].rstrip()
              #phage = tabs[9].rstrip()
              #plasmid = tabs[10].rstrip()
              #transposon = tabs[11].rstrip()
              rRNA_id = float(tabs[8].rstrip())
              #Hits.append(int(hits))
              #Lengths.append(int(length))
              if hits>0:
                  genome_interactions.append((id_a,id_b))
    
                  try:
                      ids_dict[id_a][id_b] = [int(hits),int(length),body_site_a,body_site_b,rRNA_id]
                  except KeyError:
                      ids_dict[id_a] = {id_b:[int(hits),int(length),body_site_a,body_site_b,rRNA_id]}
                  try:                
                      ids_dict[id_b][id_a] = [int(hits),int(length),body_site_b,body_site_a,rRNA_id]  
                  except KeyError: 
                      ids_dict[id_b] = {id_a:[int(hits),int(length),body_site_b,body_site_a,rRNA_id]}
                
                 
          # make ids_list a numpy array
          ids_list= np.array(ids_dict.keys())
          for key in ids_dict.keys():
              unordered_ids_list.append(int(genome_tree_img_dict[key][0]))    
          # ordered img_ids  
          working_ids_list = ids_list[np.argsort(unordered_ids_list)]
          
    for i in working_ids_list:
        if genome_tree_img_dict[i][1]==0: #firmicutes
            phylum_cols[i] ="blue"
        elif genome_tree_img_dict[i][1]==1: #actinobacteria:
            phylum_cols[i]="red"
        elif genome_tree_img_dict[i][1]==2: #synergistes
            phylum_cols[i]="yellow"
        elif genome_tree_img_dict[i][1]==3: #proteobacteria
            phylum_cols[i]="grey"
        elif genome_tree_img_dict[i][1]==4: #fusobacteria
            phylum_cols[i]="purple"
        elif genome_tree_img_dict[i][1]==5: #bacteroidetes
            phylum_cols[i]="orange"
              
    #print phylum_cols   
    
    #-----
    """Creating the network graph"""
    G=nx.Graph()
    G.add_nodes_from(working_ids_list)
    # edge width is proportional to the number of hits between genomes
    edgewidth=[]
    edgecolour=[]
    for id_a in ids_dict:
        for id_b in ids_dict[id_a]:
            #print ids_dict[id_a][id_b][0]
            #print id_a+"\t"+id_b
            #print id_a+"\t"+id_b+"\t"+str(ids_dict[id_a][id_b][0])
            G.add_edge(id_a,id_b)
            #print G.edges()[-1]
            #edgewidth.append(int(ids_dict[id_a][id_b][0]))
            #G.add_edge(id_a,id_b,weight= 1)
    #print G.edges()
    #G.add_edges_from(genome_interactions)
    #i = 0
    # Change properties of edges
    edges = []
    count = 0
    for edge in G.edges(): 
        
        #-----
        """
        # Colour edges based on vectors
        already_seen = False
        try:
            for contig in transfer_annotations[edge[0]][edge[1]]:
                for i in transfer_annotations[edge[0]][edge[1]][contig]:
                    if "phage" in i[-1].lower():
                        if already_seen == False:
                            edgecolour.append("#FF0000")
                            already_seen = True
        except KeyError:
            try:
                for contig in transfer_annotations[edge[1]][edge[0]]:
                    for i in transfer_annotations[edge[1]][edge[0]][contig]:
                        if "phage" in i[-1].lower():
                            if already_seen == False:
                                edgecolour.append("#FF0000")
                                already_seen = True
            except KeyError: 
                pass
        
        if already_seen == False:
                edgecolour.append("#D1D1D1")
        """
        #-----
        """
            try:
                for contig in transfer_annotations[edge[1]][edge[0]]:
                    for i in transfer_annotations[edge[1]][edge[0]][contig]:
                        if "phage" in i[-1].lower():
                            edgecolour.append("#FF0000")
                        else:
                            edgecolour.append("#FFFFFF")
            except Error: 
                pass
        """  
        edgewidth.append(int(ids_dict[edge[0]][edge[1]][0]))
        edges.append(edge)
        #print edge[0]+"\t"+edge[1]+"\t"+edgecolour[-1]+"\t"+str(len(edgecolour))
    
    
    values = [phylum_cols.get(node,0.25) for node in G.nodes()]
    pos= nx.spring_layout(G)
    nx.draw(G,pos,node_color=values,with_labels=args.labels,width=edgewidth,font_size=12,font_color='#006400',alpha=args.alpha)
    #nx.draw(G,pos)
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
    parser.add_argument('-i','--input_file', help="gut_oral_genomes_total_hits_length_clean.csv")
    parser.add_argument('-g','--genome_tree_file', help="Genome tree file")
    parser.add_argument('-p','--phyla', help="phyla file")
    parser.add_argument('-a','--annotation_file', help="File containing annotated IMG blast output")
    parser.add_argument('-alpha','--alpha', type=float,help="Set transparency of edges")
    parser.add_argument('-labels','--labels', type=float,default=False,help="Labels set True or False: default False")
    #parser.add_argument('-p1','--phyla_1', help="phyla file")
    #parser.add_argument('-p2','--phyla_2', help="phyla file")
    #parser.add_argument('-p3','--phyla_3', help="phyla file")
    #parser.add_argument('-p4','--phyla_4', help="phyla file")
    #parser.add_argument('-p5','--phyla_5', help="phyla file")
    #parser.add_argument('-p6','--phyla_6', help="phyla file")
    #parser.add_argument('-p7','--phyla_7', help="phyla file")
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
