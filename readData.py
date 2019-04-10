import networkx as nx
import itertools
import argparse
import os.path
try:
    import cPickle as pickle
except ImportError:
    import pickle


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Data Preprocessing')
    # parser.add_argument("-f", "--filename", required=True, type=argparse.FileType(mode='r', bufsize=-1, encoding=None, errors=None))
    parser.add_argument('-n', '--network', required=True, type=extant_file, help='global network')
    parser.add_argument('-i', '--inlist', required=True, type=str, help='in gene list')
    parser.add_argument('-o', '--outlist', required=True, type=str, help='out gene list')
    parser.add_argument('-f', '--funpath', required=True, type=str, help='pathway or function directory')
    args = parser.parse_args()

    # nodes = []
    edges = []
    with open(args.network, "r") as fp:
        #     next(fp)
        for line in fp:

            tmp = line.rstrip('\n')  # delete change line
            tmp = tmp.split("\t")  # split with backspace
            try:
                edges.append((tmp[0], tmp[1]))
            except:
                pass

    # covert (4,33) to ('UBC', 'ITCH')
    # edges_id = [(nodes_dict[inNode], nodes_dict[outNode]) for inNode, outNode in edges] #('UBC', 'ITCH')
    # whole graph
    # load nodes from Dictionary,g.size()=g.number_of_edges()=10821, g.number_of_nodes()=1886
    # node_id = nodes_dict.values()  # the id are like this, 'MAPK8','MAP2K7', 'UBC', 'ITCH'.....
    # node_key = nodes_dict.keys() # the nodes are likes: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    g = nx.Graph()


    g.add_edges_from(edges)
    # g.add_edges_from(edges_id)#('UBC', 'ITCH')
    # g.add_nodes_from(node_id)

    # write the graph to file

    nx.write_graphml(g, ".\\data\\graph.graphml")
    nx.write_gml(g, ".\\data\\graph.gml")
    nx.write_gpickle(g, ".\\data\\graph.gpickle")

    # invert the nodes mapping
    # nodes_inv_dict = dict((v, k) for k, v in nodes_dict.items())


    with open('.\\data\\funpathDir.data', 'wb') as funpathFile:
        pickle.dump(args.funpath, funpathFile)


    # with open('.\\data\\nodeDict.data', 'wb') as nodeDictFile:
    #     pickle.dump(nodes_dict, nodeDictFile)

    # with open('.\\data\\nodeDictInv.data', 'wb') as nodeDictInvFile:
    #     pickle.dump(nodes_inv_dict, nodeDictInvFile)

    # store the gene list file path directly
    with open('.\\data\\inListPath.data', 'wb') as inFile:
        pickle.dump(args.inlist, inFile)

    with open('.\\data\\outListPath.data', 'wb') as outFile:
        pickle.dump(args.outlist, outFile)
