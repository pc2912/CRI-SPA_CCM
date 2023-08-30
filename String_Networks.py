
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import numpy as np
from matplotlib import pyplot as plt
import re


def return_color_dict():
    '''define the network's edge colors + draw legend'''
    color_dict = {
                 'textmining': [0.4,0.3,0.6],
                  'neighborhood' : [0.9,0,0],
                 'fusion': [0,0.8,0],
                 'cooccurence': [0,0,0.9],
                 'homology': [0.95,0.7,0],
                 'coexpression': [0.9,0,0.8],
                 'database': [0.5,0.9,0.8]
                 }
    fig, axs= plt.subplots(1,7, figsize=(70,30))
    axs = axs.flatten()
    i=0
    for link, color in color_dict.items():
        axs[i].imshow(np.array(color_dict[link]).reshape(1,1,3).astype(float))
        axs[i].set_title(link, size=100)
        axs[i].axis('off')
        i+=1
    return color_dict
color_dict = return_color_dict()


def draw_network( G, pos , node_color,font_color=[0,0,0],figsize=(22,18) ):
    fig, ax= plt.subplots(figsize=figsize)
    shift =0.1#0.04
    label_pos = {node:(x*(1.05+shift*(i%2)),y*(1.05+shift*(i%2))) for i, (node, (x,y)) in zip(range(len(pos)),pos.items())}
    #we first draw our nodes
    nx.draw_networkx_labels(G, pos=label_pos,
                            ax=ax,font_color = font_color )
    
    node_pos = {node:(x*(1+0.03*(i%2)),y*(1+0.03*(i%2))) for i, (node, (x,y)) in zip(range(len(pos)),pos.items())}
    nx.draw_networkx_nodes(G, pos=pos, ax=ax, node_color= node_color, node_size=100)
    all_rads={1:[0],2:[0.1,-0.1],3:[0,0.1,-0.1],4:[0.08,-0.08,0.14,-.14],
         5:[0, 0.08,-0.08,0.14,-.14], 6:[ 0.03,-0.03,0.08,-0.08,0.14,-.14],
          7:[0,0.03,-0.03, 0.08,-0.08,0.14,-.14]}
    
    for edge in set(G.edges()):
        edge_attributes = G.get_edge_data(edge[0],edge[1])
        n_edges= len(edge_attributes )
        for n, (metric, attributes) in enumerate(edge_attributes.items()):
            
            rad = all_rads[n_edges][n]
            nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(edge[0],edge[1],metric)],
                           node_size=100,
                           edge_color= attributes['color'] ,
                           connectionstyle='arc3,rad='+str(rad),
                           width =  attributes['width'],
                           arrowstyle='-')
    plt.grid(b=None)
    plt.axis('off')
    return fig, ax


def draw_tonsof_network( G, pos , node_color,fig,ax, font_color=[0,0,0], node_size=100):
    shift =0.1#0.04
    label_pos = {node:(x*(1.05+shift*(i%2)),y*(1.05+shift*(i%2))) for i, (node, (x,y)) in zip(range(len(pos)),pos.items())}
    #we first draw our nodes
   # nx.draw_networkx_labels(G, pos=label_pos,
      #                      ax=ax,font_color = font_color )
    
    node_pos = {node:(x*(1+0.03*(i%2)),y*(1+0.03*(i%2))) for i, (node, (x,y)) in zip(range(len(pos)),pos.items())}
    nx.draw_networkx_nodes(G, pos=pos, ax=ax, node_color= node_color, node_size=node_size)
    all_rads={1:[0],2:[0.1,-0.1],3:[0,0.1,-0.1],4:[0.08,-0.08,0.14,-.14],
         5:[0, 0.08,-0.08,0.14,-.14], 6:[ 0.03,-0.03,0.08,-0.08,0.14,-.14],
          7:[0,0.03,-0.03, 0.08,-0.08,0.14,-.14]}
    
    for edge in set(G.edges()):
        edge_attributes = G.get_edge_data(edge[0],edge[1])
        n_edges= len(edge_attributes )
        for n, (metric, attributes) in enumerate(edge_attributes.items()):
            
            rad = all_rads[n_edges][n]
            nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(edge[0],edge[1],metric)],
                           node_size=100,
                           edge_color= attributes['color'] ,
                           connectionstyle='arc3,rad='+str(rad),
                           width =  attributes['width'],
                           arrowstyle='-')
    return fig, ax

def draw_network_spring( G, pos , node_color,font_color=[0,0,0],figsize=(22,18) ):
    fig, ax= plt.subplots(figsize=figsize)
    plt.grid(b=None)
    label_pos = {node:(x+figsize[0]*0.001,y+figsize[1]*0.001) for i, (node, (x,y)) in zip(range(len(pos)),pos.items())}
     
        #we first draw our nodes
    nx.draw_networkx_labels(G, pos=label_pos,
                            ax=ax,font_color = font_color )
    

    nx.draw_networkx_nodes(G, pos=pos, ax=ax, node_color= node_color, node_size=100)
    all_rads={1:[0],2:[0.1,-0.1],3:[0,0.1,-0.1],4:[0.08,-0.08,0.14,-.14],
         5:[0, 0.08,-0.08,0.14,-.14], 6:[ 0.03,-0.03,0.08,-0.08,0.14,-.14],
          7:[0,0.03,-0.03, 0.08,-0.08,0.14,-.14]}
    
    for edge in set(G.edges()):
        edge_attributes = G.get_edge_data(edge[0],edge[1])
        n_edges= len(edge_attributes )
        for n, (metric, attributes) in enumerate(edge_attributes.items()):
            
            rad = all_rads[n_edges][n]
            nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(edge[0],edge[1],metric)],
                           node_size=100,
                           edge_color= attributes['color'] ,
                           connectionstyle='arc3,rad='+str(rad),
                           width =  attributes['width'],
                           arrowstyle='-')
            
    plt.axis('off')
    
            

def build_network(gene_group,gene2systematic,systematic2gene, String, G= nx.MultiGraph(directed=False)):

    i1 = String.protein1.isin(['4932.'+gene2systematic[g] for g in gene_group])

    i2 =  String.protein2.isin(['4932.'+gene2systematic[g] for g in gene_group])

    (i1 &i2).sum()

    subString= String[i1 &i2]
    #we go through all A,B B,A pairs in the group of gene
    # we created (A,B, Metric) for all linking metrics in string
    # we created (B,A, Metric) for all linking metrics in string
    # if a link is only in A,B or only in B,A it is recorded
    #if link is both in A,B and B,A the properties of B,A overwrites the link properties
    for i in range(len(subString)):
        gene1= subString.iloc[i].protein1[5:]
        gene1 = systematic2gene[gene1]
        gene2= subString.iloc[i].protein2[5:]
        gene2 = systematic2gene[gene2]

        G.add_node(gene1)
        G.add_node(gene2)
        for  metric, color in color_dict.items():
            # set confindence score higher for textmining which tends to link everything
            #see  https://string-db.org/cgi/help?sessionId=bQhE0Cc56XQA
            confidence = 700 if metric == 'textmining' else 700
            if subString.iloc[i][metric]> confidence: #check that their is a link>0 between the 2 genes at a given metric
                #we check that B,A,metric not an edge before adding A,B, metric
                # if A,B, metric an edge we just overwrite it
                if G.has_edge(gene2, gene1, metric) == False:
                    G.add_edge(gene1, gene2, metric) #B,A will overwrite A,B with the same properties (assuming String is symmetric)
                    G.edges[(gene1, gene2, metric)]['width'] =subString.iloc[i][metric]/400
                    G.edges[(gene1, gene2, metric)]['color'] = color
                
     
    return G