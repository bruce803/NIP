import networkx as nx
import numpy as np
from numpy import array, empty
from scipy.stats import hypergeom as hg
import pandas as pd
import pylab as plt
import argparse
import itertools
import sys
import csv
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle


# merge

def node_inter_g(dic, gNodeList):
    '''get the intersection between pathway and node_layer. return the intersection as dict, 
    key:pathway name, value:common genelist.
    '''
    mergeNode = []
    for k, v in dic.items():
        intersec = geneListInter(gNodeList, list(v))
        if len(intersec) > 1:
            mergeNode.append((k, intersec))
    return dict(mergeNode)


def node_both_inter_g(dic, list1, list2):
    '''return the nodes from both dic&list1, and dic&list2'''
    mergeNode = []
    for k, v in dic.items():
        intersec1 = geneListInter(list1, list(v))
        intersec2 = geneListInter(list2, list(v))
        if intersec1 and intersec2:
            intersec = intersec1 + intersec2
            if len(intersec) > 3:
                mergeNode.append((k, intersec))
    return dict(mergeNode)


def noninter_patition(sorted_funInterNode):
    '''find the non-intersected partitions of nodes set to merge'''
    geneNode = set(sorted_funInterNode[0][1])
    noninterList = [sorted_funInterNode[0]]
    for item in sorted_funInterNode:
        if not (set(item[1]) & set(geneNode)):
            geneNode = geneNode.union(item[1])
            noninterList.append(item)
    return noninterList


def merge_nodes(G, nodes, new_node, attr_dict=None, **attr):
    """
    Merges the selected `nodes` of the graph G into one `new_node`,
    meaning that all the edges that pointed to or from one of these
    `nodes` will point to or from the `new_node`.
    attr_dict and **attr are defined as in `G.add_node`.
    """

    G.add_node(new_node, attr_dict, **attr)  # Add the 'merged' node

    for n1, n2, data in G.edges(data=True):
        # For all edges related to one of the nodes to merge,
        # make an edge going to or coming from the `new gene`.
        if n1 in nodes:
            G.add_edge(new_node, n2, data)
        elif n2 in nodes:
            G.add_edge(n1, new_node, data)


def copy_merge_nodes(G, G_copy, nodes, new_node, attr_dict=None, **attr):
    """
    Merge the nodes to new_node, but copy this merge to a new graph G_copy
    attr_dict and **attr are defined as in `G.add_node`.
    """

    #     G_copy.add_node(new_node, attr_dict, **attr) # Add the 'merged' node

    for n1, n2 in G.edges():
        # For all edges related to one of the nodes to merge,
        # make an edge going to or coming from the `new gene`.
        if n1 in nodes:
            G_copy.add_edge(new_node, n2)
        elif n2 in nodes:
            G_copy.add_edge(n1, new_node)


# load gene list
def readGeneFromFile(file_path):
    gene_list = []
    with open(file_path, "r") as fp:
        for line in fp:
            tmp = line.rstrip('\n')
            try:
                gene_list.append(tmp)
            except:pass
    return gene_list


def geneListInter(list1, list2):
    '''get the instersection of two sets of nodes or edges. The nodes (char) or edges(tuple) should be stored as List.
    For example list1=geneList, and list2=netNodes are genes' names'''
    return list(set(list1)&set(list2))


# convert from NAME to id number ('key' in nodes_dict)
def name2id(inter_immune,nodes_inv_dict):
    inter_immune_id = []
    for i in range(len(inter_immune)):
        inter_immune_id.append(nodes_inv_dict[inter_immune[i]])
    return inter_immune_id


def complete_bipartite(list1,list2):
    '''generate a complete bipartite graph between list1 and list2'''
    c12 = list(itertools.product(list1, list2)) #考虑(1,4)和（4，1）两种情况
    c21 = list(itertools.product(list2, list1))
    return c12 + c21


def read_dict_from_folder(indir):
    '''read the function or pathway into dictionary, key:pathway name, value:pathway'''
    pathFunDict = {}
    for root, dirs, filenames in os.walk(indir):
        for fn in filenames:
            path = os.path.join(indir,fn)
            geneList = readGeneFromFile(path)
            pathFunDict[fn]= geneList
    return pathFunDict


def str2bool(v):
    '''parse boolen values with argparse'''
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def correct_pvalues_for_multiple_testing(pvalues, correction_type="Benjamini-Hochberg"):
    """ 
    p values to fdrs
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """

    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n - rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n / rank) * pvalue)
        for i in range(0, int(n) - 1):
            if new_values[i] < new_values[i + 1]:
                new_values[i + 1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


def allSimplePathkHop(graph, source, target, k):
    ''' return all the k hops simple paths from nodes in the source list to the target list.
    Parameters
    ----------
    graph : NetworkX graph

    source : node list
       Starting node for path

    target : node list
       Ending node for path. If provided only predecessors between
       source and target are returned

    cutoff : integer, optional
        Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    pred : list of paths
    '''
    hop_k = []
    for node_s in source:
        for node_t in target:
            path_kHop = list(nx.all_simple_paths(graph, node_s, node_t, cutoff=k))
            hop_k = hop_k + path_kHop

    hop_k = [item for item in hop_k if len(item)==k+1] # only keep the paths that their length==k+1
    return hop_k


def get_edges(pathList,k):
    '''extract the all the edges from pathList'''
    edge_layer_all = []
    for i in range(len(pathList)):
        for j in range(k):
            edge_layer_all.append(tuple([np.array(pathList)[i,j],np.array(pathList)[i,j+1]]))
    return edge_layer_all


def extract_node_edges_from_file(fpath1, fpath2, graph, k):  ###update
    '''output the results of subgraph analysis, including the edges, nodes in different layer of the subgraph.
    k: length of path'''
    gene_list1 = readGeneFromFile(fpath1)
    gene_list2 = readGeneFromFile(fpath2)

    # invert the nodes mapping
    key = list(graph.nodes)  # len(key)=1886, load the nodes of immune net

    list1_inter_immune = geneListInter(gene_list1, key)
    list2_inter_immune = geneListInter(gene_list2, key)
    # list1_inter_immune_id = name2id(list1_inter_immune, nodes_inv_dict)
    # list2_inter_immune_id = name2id(list2_inter_immune, nodes_inv_dict)

    if list1_inter_immune and list2_inter_immune:
        pathList = allSimplePathkHop(graph, list1_inter_immune, list2_inter_immune, k)
        # pathList = allSimplePathkHop(graph, list1_inter_immune_id, list2_inter_immune_id, k)
        edge_layer_all = get_edges(pathList, k)

        if pathList:
            if k == 1:
                return list(set(np.array(pathList)[:, 0])), list(set(np.array(pathList)[:, 1])), edge_layer_all
            elif k == 2:
                return list(set(np.array(pathList)[:, 0])), list(set(np.array(pathList)[:, 1])), \
                       list(set(np.array(pathList)[:, 2])), edge_layer_all
            elif k == 3:
                return list(set(np.array(pathList)[:, 0])), list(set(np.array(pathList)[:, 1])), \
                       list(set(np.array(pathList)[:, 2])), list(set(np.array(pathList)[:, 3])), edge_layer_all
        else:
            print("There are no interaction between the two groups!!!")
            return
    else:
        sys.exit("There are no intersections between gene list and global network")


def node_inter_g_fdr(dic, gNodeList, fdrThreshByUser):
    '''有fdr选项的合并函数'''
    mergeNode = []
    len_gNodeList = len(gNodeList)
    for k, v in dic.items():
        intersec = geneListInter(gNodeList, list(v))
        if len(intersec) >= 3:  # step 1, filtering with length
            pVal = hg.sf(len(intersec), 19718, len(v), len_gNodeList)
            mergeNode.append([k, intersec, pVal])
    mergeNode_df = pd.DataFrame(mergeNode, columns=['PwID', 'nodeList', 'pVal'])

    mergeNode_df['fdr'] = pd.Series(correct_pvalues_for_multiple_testing(list(mergeNode_df.pVal)),
                                        index=mergeNode_df.index)
    mergeNode_df = mergeNode_df[mergeNode_df.fdr > fdrThreshByUser]  # step2: filtering with fdr

    return dict(zip(mergeNode_df.PwID, mergeNode_df.nodeList))


def node_both_inter_g_fdr(dic, list1, list2, fdrThreshByUser):
    '''return the nodes from both dic&list1, and dic&list2
    filtering with fdr threshold'''
    mergeNode = []
    len_gNodeList = len(list1) + len(list2)
    for k, v in dic.items():
        intersec1 = geneListInter(list1, list(v))
        intersec2 = geneListInter(list2, list(v))
        if intersec1 and intersec2:
            intersec = intersec1 + intersec2
            if len(intersec) > 3:
                pVal = hg.sf(len(intersec), 19718, len(v), len_gNodeList)
                mergeNode.append([k, intersec, pVal])

    mergeNode_df = pd.DataFrame(mergeNode, columns=['PwID', 'nodeList', 'pVal'])
    mergeNode_df['fdr'] = pd.Series(correct_pvalues_for_multiple_testing(list(mergeNode_df.pVal)),
                                        index=mergeNode_df.index)
    mergeNode_df = mergeNode_df[mergeNode_df.fdr > fdrThreshByUser]  # step2: filtering with fdr
    return dict(zip(mergeNode_df.PwID, mergeNode_df.nodeList))

#load network
g = nx.read_gpickle('.\\data\\graph.gpickle')

# load nodes dictionary
# with open('.\\data\\nodeDict.data','rb') as f1:
#     nodes_dict = pickle.load(f1)

# load inverse nodes dictionary
# with open('.\\data\\nodeDictInv.data','rb') as f2:
#     nodes_inv_dict = pickle.load(f2)


with open('.\\data\\inListPath.data','rb') as f3:
    cle_up = pickle.load(f3)

with open('.\\data\\outListPath.data','rb') as f4:
    ple_up = pickle.load(f4)

with open('.\\data\\funpathDir.data','rb') as f5:
    funpathDir = pickle.load(f5)



# read pathway or function , we union them together
# funDict = read_dict_from_folder(funDir) #396
pathwayDict = read_dict_from_folder(funpathDir) #1807


def experiment1(src, dst, graph, nodeslist):
    '''1 hop'''

    gene_list1 = readGeneFromFile(src)
    gene_list2 = readGeneFromFile(dst)

    graph_nodes = list(graph.nodes)
    list1_inter_immune = geneListInter(gene_list1, graph_nodes)
    list2_inter_immune = geneListInter(gene_list2, graph_nodes)
    if list1_inter_immune and list2_inter_immune:
        pathList1hops = allSimplePathkHop(graph, list1_inter_immune, list2_inter_immune, 1)
        edge_layer_all1hops = get_edges(pathList1hops, 1)
        print(edge_layer_all1hops)

    if edge_layer_all1hops:
        n_layer0 = list(set(np.array(pathList1hops)[:, 0]))
        n_layer1 = list(set(np.array(pathList1hops)[:, 1]))
    else:
        print('Gene List 1 is:', list1_inter_immune)
        print('Gene List 2 is:', list2_inter_immune)
        print('There are no 1 hop connections.')
        sys.exit()

    g_1hop= nx.Graph()
    g_1hop.add_edges_from(edge_layer_all1hops)

    if not n_layer1:
        print('There are no interactions between the two gene groups!')
        return

    else:
        if nodeslist:  # act delete
            g_1hop.remove_nodes_from(nodeslist)
            n_layer0 = geneListInter(n_layer0, g_1hop.nodes())
            n_layer1 = geneListInter(n_layer1, g_1hop.nodes())

        x0 = np.repeat(0, len(n_layer0))
        x1 = np.repeat(6, len(n_layer1))

        y0 = np.linspace(0, 100, len(n_layer0))
        y1 = np.linspace(0, 100, len(n_layer1))

        pos_0 = np.column_stack([x0, y0])
        pos_1 = np.column_stack([x1, y1])

        pos_0 = dict(zip(n_layer0, pos_0))
        pos_1 = dict(zip(n_layer1, pos_1))

        nx.draw_networkx_nodes(g_1hop, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_1hop, pos_1, nodelist=n_layer1, node_color='c', node_size=300, alpha=0.8)
        # nx.draw_networkx_nodes(g_pleuucle_1hop,pos_3,nodelist=n_layer3,node_color='c',node_size=300,alpha=0.8)

        # edges
        pos = {}
        pos.update(pos_0)
        pos.update(pos_1)
        nx.draw_networkx_edges(g_1hop, pos, edgelist=nx.edges(g_1hop), width=1, alpha=0.8, edge_color='k')
        nx.draw_networkx_labels(g_1hop, pos, font_size=16, font_family='sans-serif')

        #     fig = plt.gcf()
        #     plt.draw()
        #     fig.set_size_inches(8, 5)
        plt.axis('off')

        # save graph and statistics of graph into file
        nx.write_graphml(g_1hop, ".\\data\\g_1hop.graphml")
        # nx.write_gml(g_1hop, ".\\data\\g_1hop.gml")
        nx.write_gpickle(g_1hop, ".\\data\\g_1hop.gpickle")

        with open('.\\data\\g_1hop_stat.csv', 'w+') as f:
            writer = csv.DictWriter(f, fieldnames=['nodes_layer0', 'nodes_layer1'])
            writer.writeheader()
            csv.writer(f).writerows(itertools.zip_longest(n_layer0, n_layer1))

        fig = plt.gcf()
        plt.draw()
        fig.set_size_inches(8, 5)
        plt.axis('off')
        # plt.title('The merged subnetwork from inList up to OutList')
        plt.show()
        # fig.savefig('g_pathway_merge.png', dpi=900)
        # fig.savefig('g_function_merge.png', dpi=900)


def experiment2(src, dst, graph, fdr, nodeslist, merge):
    '''2 hops'''
    gene_list1 = readGeneFromFile(src)
    gene_list2 = readGeneFromFile(dst)

    graph_nodes = list(graph.nodes)
    list1_inter_immune = geneListInter(gene_list1, graph_nodes)
    list2_inter_immune = geneListInter(gene_list2, graph_nodes)

    if list1_inter_immune and list2_inter_immune:
        pathList2hops = allSimplePathkHop(graph, list1_inter_immune, list2_inter_immune, 2)
        edge_layer_all2hops = get_edges(pathList2hops, 2)

    if pathList2hops:
        n_layer0 = list(set(np.array(pathList2hops)[:, 0]))
        n_layer1 = list(set(np.array(pathList2hops)[:, 1]))
        n_layer2 = list(set(np.array(pathList2hops)[:, 2]))
    else:
        print('Gene List 1 is:', list1_inter_immune)
        print('Gene List 2 is:', list2_inter_immune)
        print('There are no 2 hops connections.')
        sys.exit()


    g_2hop= nx.Graph()
    # g_2hop.add_nodes_from(node_all)
    g_2hop.add_edges_from(edge_layer_all2hops)

    node_02_copy = n_layer0 + n_layer2  # 37
    corssNodes = geneListInter(n_layer1, node_02_copy)  # 中间层和src、dst层重复的nodes
    n_layer1 = [node for node in n_layer1 if node not in corssNodes]  # 把重复的notes从中间层踢掉

    g_2hop_copy = nx.Graph()
    g_2hop_copy.add_nodes_from(node_02_copy)

    if not merge: # update, no merge, no fdr
        if nodeslist:
            g_2hop.remove_nodes_from(nodeslist)

        x0 = np.repeat(0, len(n_layer0))
        x1 = np.repeat(3, len(n_layer1))
        x2 = np.repeat(6, len(n_layer2))

        y0 = np.linspace(0, 100, len(n_layer0))
        y1 = np.linspace(0, 100, len(n_layer1))
        y2 = np.linspace(0, 100, len(n_layer2))

        pos_0 = np.column_stack([x0, y0])
        pos_1 = np.column_stack([x1, y1])
        pos_2 = np.column_stack([x2, y2])

        #     pos_0 = rescale_layout(pos_0, scale=1)
        #     pos_1 = rescale_layout(pos_1, scale=1)
        #     pos_2 = rescale_layout(pos_2, scale=1)

        pos_0 = dict(zip(n_layer0, pos_0))
        pos_1 = dict(zip(n_layer1, pos_1))
        pos_2 = dict(zip(n_layer2, pos_2))

        nx.draw_networkx_nodes(g_2hop, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_2hop, pos_1, nodelist=n_layer1, node_color='r', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_2hop, pos_2, nodelist=n_layer2, node_color='c', node_size=300, alpha=0.8)
        # nx.draw_networkx_nodes(g_pleuucle_1hop,pos_3,nodelist=n_layer3,node_color='c',node_size=300,alpha=0.8)

        # edges
        pos = {}
        pos.update(pos_0)
        pos.update(pos_1)
        pos.update(pos_2)
        # pos.update(pos_3)
        nx.draw_networkx_edges(g_2hop, pos, edgelist=nx.edges(g_2hop), width=1, alpha=0.8, edge_color='k')
        nx.draw_networkx_labels(g_2hop, pos, font_size=16, font_family='sans-serif')

        plt.axis('off')
        plt.show()


        # plt.title('The merged subnetwork from inList up to outList up')
        # fig.savefig('g_pathway_merge.png', dpi=900)
        # fig.savefig('g_function_merge.png', dpi=900)

        # save graph and statistics of graph into file
        nx.write_graphml(g_2hop, ".\\data\\g_2hop.graphml")
        nx.write_gml(g_2hop, ".\\data\\g_2hop.gml")
        nx.write_gpickle(g_2hop, ".\\data\\g_2hop.gpickle")

        with open('.\\data\\g_2hop_stas.csv', 'w+') as f:
            writer = csv.DictWriter(f, fieldnames=['nodes_layer0', 'nodes_layer1', 'nodes_layer2'])
            writer.writeheader()
            csv.writer(f).writerows(itertools.zip_longest(n_layer0, n_layer1, n_layer2))

    elif merge:
        if not nodeslist:  # not act delete

            if fdr is None:
                inter_node = node_inter_g(pathwayDict, n_layer1) # pathway
                # inter_node = node_inter_g(funDict, n_layer1) #function
            else:
                inter_node = node_inter_g_fdr(pathwayDict, n_layer1, fdr)  # pathway

            for k, v in inter_node.items():
                copy_merge_nodes(g_2hop, g_2hop_copy, v, k)

        else:             # act delete
            # g_2hop.remove_nodes_from(nodeslist)
            # n_layer1 = geneListInter(g_2hop.edges(), n_layer1)

            if fdr is None:
                inter_node = node_inter_g(pathwayDict, n_layer1)

            else:
                inter_node = node_inter_g_fdr(pathwayDict, n_layer1, fdr)

            for k, v in inter_node.items():
                copy_merge_nodes(g_2hop, g_2hop_copy, n_layer0, n_layer2, v, k)

            g_2hop_copy.remove_nodes_from(nodeslist)  # delete nodes in nodeslist

        node_layer1_copy = [item for item in g_2hop_copy.nodes() if item not in node_02_copy] # g_3hop_copy.nodes - node_copy

        lens = [len(n_layer0), len(node_layer1_copy), len(n_layer2)]

        pos_0 = {}
        x0 = 1
        # const = 2
        y = 100
        for i in range(len(n_layer0)):
            pos_0[n_layer0[i]] = [x0, y - i * (max(lens) / len(n_layer0))]

        x1 = 3
        pos_1 = {}
        for j in range(len(node_layer1_copy)):
            pos_1[node_layer1_copy[j]] = [x1, y - j * (max(lens) / len(node_layer1_copy))]

        x2 = 5
        pos_2 = {}
        for n in range(len(n_layer2)):
            pos_2[n_layer2[n]] = [x2, y - n * (max(lens) / len(n_layer2))]

        nx.draw_networkx_nodes(g_2hop_copy, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_2hop_copy, pos_1, nodelist=node_layer1_copy, node_color='r', node_size=600, alpha=0.8)
        nx.draw_networkx_nodes(g_2hop_copy, pos_2, nodelist=n_layer2, node_color='c', node_size=300, alpha=0.8)
        # nx.draw_networkx_nodes(g_pleuucle_1hop,pos_3,nodelist=n_layer3,node_color='c',node_size=300,alpha=0.8)

        # edges
        pos = {}
        pos.update(pos_0)
        pos.update(pos_1)
        pos.update(pos_2)
        # pos.update(pos_3)
        nx.draw_networkx_edges(g_2hop_copy, pos, edgelist=nx.edges(g_2hop_copy), width=1, alpha=0.8, edge_color='k')
        nx.draw_networkx_labels(g_2hop_copy, pos, font_size=16, font_family='sans-serif')


        fig = plt.gcf()
        plt.draw()
        fig.set_size_inches(8, 5)
        plt.axis('off')
        plt.show()
        # plt.title('The merged subnetwork from inList up to outList')
        # fig.savefig('g_pathway_merge.png', dpi=900)
        # fig.savefig('g_function_merge.png', dpi=900)

        # save graph and statistics of graph into file
        nx.write_graphml(g_2hop_copy, ".\\data\\g_2hop.graphml")
        nx.write_gml(g_2hop_copy, ".\\data\\g_2hop.gml")
        nx.write_gpickle(g_2hop_copy, ".\\data\\g_2hop.gpickle")

        with open('.\\data\\g_2hop_stas.csv', 'w+') as f:
            writer = csv.DictWriter(f, fieldnames=['nodes_layer0', 'nodes_layer1', 'nodes_layer2'])
            writer.writeheader()
            csv.writer(f).writerows(itertools.zip_longest(n_layer0, node_layer1_copy, n_layer2))


def experiment3(src, dst, graph, fdr, nodeslist, merge):
    ''' 3 hops'''
    gene_list1 = readGeneFromFile(src)
    gene_list2 = readGeneFromFile(dst)

    graph_nodes = list(graph.nodes)
    list1_inter_immune = geneListInter(gene_list1, graph_nodes)
    list2_inter_immune = geneListInter(gene_list2, graph_nodes)

    if list1_inter_immune and list2_inter_immune:
        pathList3hops = allSimplePathkHop(graph, list1_inter_immune, list2_inter_immune, 3)
        # pathList = allSimplePathkHop(graph, list1_inter_immune_id, list2_inter_immune_id, k)
        edge_layer_all3hops = get_edges(pathList3hops, 3)

    if pathList3hops:
        n_layer0 = list(set(np.array(pathList3hops)[:, 0]))
        n_layer1 = list(set(np.array(pathList3hops)[:, 1]))
        n_layer2 = list(set(np.array(pathList3hops)[:, 2]))
        n_layer3 = list(set(np.array(pathList3hops)[:, 3]))
    else:
        print('There are no 3 hops connections.')
        sys.exit()


    g_3hop= nx.Graph()
    # g_3hop.add_nodes_from(node_all)
    g_3hop.add_edges_from(edge_layer_all3hops)

    node_03_copy = n_layer0 + n_layer3  # 最外面两层的节点
    node_12_copy = n_layer1 + n_layer2  # nodes of middle two layers
    corssNodes = geneListInter(node_12_copy, node_03_copy)  # 中间层和src、dst层重复的nodes
    n_layer1 = [node for node in n_layer1 if node not in corssNodes]  # 把重复的notes从中间层踢掉
    n_layer2 = [node for node in n_layer2 if node not in corssNodes]  # 把重复的notes从中间层踢掉

    crossNodes_12 = geneListInter(n_layer1, n_layer2)  # 中间两层的交集

    if (len(n_layer1) >= len(n_layer2)):
        n_layer1 = [node for node in n_layer1 if node not in crossNodes_12]
    else:
        n_layer2 = [node for node in n_layer2 if node not in crossNodes_12]

    if not merge:
        if nodeslist:
            g_3hop.remove_nodes_from(nodeslist)


        x0 = np.repeat(0, len(n_layer0))
        x1 = np.repeat(2, len(n_layer1))
        x2 = np.repeat(4, len(n_layer2))
        x3 = np.repeat(6, len(n_layer3))

        y0 = np.linspace(0, 100, len(n_layer0))
        y1 = np.linspace(0, 100, len(n_layer1))
        y2 = np.linspace(0, 100, len(n_layer2))
        y3 = np.linspace(0, 100, len(n_layer3))

        pos_0 = np.column_stack([x0, y0])
        pos_1 = np.column_stack([x1, y1])
        pos_2 = np.column_stack([x2, y2])
        pos_3 = np.column_stack([x3, y3])

        #     pos_0 = rescale_layout(pos_0, scale=1)
        #     pos_1 = rescale_layout(pos_1, scale=1)
        #     pos_2 = rescale_layout(pos_2, scale=1)

        pos_0 = dict(zip(n_layer0, pos_0))
        pos_1 = dict(zip(n_layer1, pos_1))
        pos_2 = dict(zip(n_layer2, pos_2))
        pos_3 = dict(zip(n_layer3, pos_3))

        nx.draw_networkx_nodes(g_3hop, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_3hop, pos_1, nodelist=n_layer1, node_color='r', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_3hop, pos_2, nodelist=n_layer2, node_color='r', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_3hop, pos_3, nodelist=n_layer3, node_color='c', node_size=300, alpha=0.8)

        #         edges
        pos = {}
        pos.update(pos_0)
        pos.update(pos_1)
        pos.update(pos_2)
        pos.update(pos_3)
        nx.draw_networkx_edges(g_3hop, pos, edgelist=nx.edges(g_3hop), width=1, alpha=0.8,
                               edge_color='k')
        nx.draw_networkx_labels(g_3hop, pos, font_size=10, font_family='sans-serif')

        #         fig = plt.gcf()
        #         plt.draw()
        #         fig.set_size_inches(8, 5)
        plt.axis('off')
        plt.show()
        # plt.title('The merged subnetwork from inList up to outList')
        # fig.savefig('g_pathway_merge.png', dpi=900)
        # fig.savefig('g_function_merge.png', dpi=900)

        # save graph and statistics of graph into file
        # nx.write_graphml(g_3hop, ".\\data\\g_3hop.graphml")
        # nx.write_gml(g_3hop, ".\\data\\g_3hop.gml")
        # nx.write_gpickle(g_3hop, ".\\data\\g_3hop.gpickle")

        # with open('.\\data\\g_3hop_stas.csv', 'w+') as f:
        #     writer = csv.DictWriter(f, fieldnames=['nodes_layer0', 'nodes_layer1', 'nodes_layer2' 'nodes_layer3'])
        #     writer.writeheader()
        #     csv.writer(f).writerows(itertools.zip_longest(n_layer0, n_layer1, n_layer2, n_layer3))


    # return n_layer0, n_layer1, n_layer2, n_layer3, edge_all, g_3hop_copy,
    # inter_node = node_both_inter_g_fdr(pathwayDict, n_layer1, n_layer2, fdr)

    elif merge:
        node_03_copy = n_layer0 + n_layer3  # 37
        g_3hop_copy = nx.Graph()
        g_3hop_copy.add_nodes_from(node_03_copy)
        # True
        if fdr is None:
            inter_node = node_both_inter_g(pathwayDict, n_layer1, n_layer2) #198, pathway
            # inter_node = node_both_inter_g(funDict, n_layer1, n_layer2)  # 28, function
        else:
            inter_node = node_both_inter_g_fdr(pathwayDict, n_layer1, n_layer2, fdr)


        for k, v in inter_node.items():
            #     merge_id = [nodes_inv_dict[gene_name] for gene_name in v] # display with id
            copy_merge_nodes(g_3hop, g_3hop_copy, n_layer0, n_layer3, v, k)

        if nodeslist:
            g_3hop_copy.remove_nodes_from(nodeslist)

        node_layer1_copy = [item for item in g_3hop_copy.nodes() if item not in node_03_copy] # g_3hop_copy.nodes - node_copy

        lens = [len(n_layer0), len(node_layer1_copy), len(n_layer3)]

        pos_0 = {}
        x0 = 1
        const = 2
        y = 100
        for i in range(len(n_layer0)):
            pos_0[n_layer0[i]] = [x0, y - i * (max(lens) / len(n_layer0))]

        x1 = 3
        pos_1 = {}
        for j in range(len(node_layer1_copy)):
            pos_1[node_layer1_copy[j]] = [x1, y - j * (max(lens) / len(node_layer1_copy))]

        x2 = 5
        pos_2 = {}
        for n in range(len(n_layer3)):
            pos_2[n_layer3[n]] = [x2, y - n * (max(lens) / len(n_layer3))]

        nx.draw_networkx_nodes(g_3hop_copy, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_3hop_copy, pos_1, nodelist=node_layer1_copy, node_color='r', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(g_3hop_copy, pos_2, nodelist=n_layer3, node_color='c', node_size=300, alpha=0.8)
        # nx.draw_networkx_nodes(g_pleuucle_1hop,pos_3,nodelist=n_layer3,node_color='c',node_size=300,alpha=0.8)

        # edges
        pos = {}
        pos.update(pos_0)
        pos.update(pos_1)
        pos.update(pos_2)
        # pos.update(pos_3)
        nx.draw_networkx_edges(g_3hop_copy, pos, edgelist=nx.edges(g_3hop_copy), width=1, alpha=0.8,
                               edge_color='k')
        nx.draw_networkx_labels(g_3hop_copy, pos, font_size=8, font_family='sans-serif')

        fig = plt.gcf()
        plt.draw()
        fig.set_size_inches(8, 5)
        plt.axis('off')
        plt.show()
        # plt.title('The merged subnetwork from inList to outList')
        # fig.savefig('g_pathway_merge.png', dpi=900)
        # fig.savefig('g_function_merge.png', dpi=900)

        # save graph and statistics of graph into file
        nx.write_graphml(g_3hop_copy, ".\\data\\g_3hop.graphml")
        nx.write_gml(g_3hop_copy, ".\\data\\g_3hop.gml")
        nx.write_gpickle(g_3hop_copy, ".\\data\\g_3hop.gpickle")

        with open('.\\data\\g_3hop_stas.csv', 'w+') as f:
            writer = csv.DictWriter(f, fieldnames=['nodes_layer0', 'nodes_layer1', 'nodes_layer2'])
            writer.writeheader()
            csv.writer(f).writerows(itertools.zip_longest(n_layer0, node_layer1_copy, n_layer2,))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Network Interactor Parser')

    parser.add_argument('--layer', '-l', type=int, choices=[1, 2, 3], default=2, help="Numbers of hops" )
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument("--merge", "-m", action="store_true", default=False,
    #                    help="Merge the intersected nodes accroding to the pathways or functions")
    parser.add_argument("--merge", "-m", action="store_true", default=False,
                       help="Merge the intersected nodes accroding to the pathways or functions")
    parser.add_argument('--delete', '-d', nargs='+', type=str, help="Delete nodes or edges")
    parser.add_argument('--pvalue', '-p', default=None, nargs='?', dest="fdr", type=float,
                        help="Filtering the merged nodes with FDR")

    # args = parser.parse_args()
    # args = parser.parse_args('-l 3 -m -p 0.05'.split())
    args = parser.parse_args('-l 2 -m '.split())
    #

    if args.merge:
        if args.layer == 1:
            experiment1(cle_up, ple_up, g, args.delete)
        elif args.layer == 2:
            experiment2(cle_up, ple_up, g, args.fdr, args.delete, merge=True)
        elif args.layer == 3:
            experiment3(cle_up, ple_up, g, args.fdr, args.delete, merge=True)
            print('merge the middle nodes with fdr {0}'.format(args.fdr))
    elif not args.merge:
        if args.layer == 1:
            experiment1(cle_up, ple_up, g, args.delete)
        elif args.layer == 2:
            experiment2(cle_up, ple_up, g, args.fdr, args.delete, merge=False)
        elif args.layer == 3:
            experiment3(cle_up, ple_up, g, args.fdr, args.delete, merge=False)
            print('merge both the two layer node with fdr {0}'.format(args.fdr))