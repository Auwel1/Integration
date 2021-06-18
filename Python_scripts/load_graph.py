import argparse
from xml.dom import minidom as mdom
from xml.etree import ElementTree as ET
import xml
import pandas as pd
import igraph as ig
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import HMA_connect as hma
import utils
from tqdm import tqdm
import re

# parser = argparse.ArgumentParser
#
# parser.add_argument('-xml', '--xmldoc', help='xml metabolic file')
# parser.add_argument('-t', '--trim', help='path to a node list to eliminate from graph')
# parser.add_argument('-o', '--output', help='output directory path')
# parser.add_argument('-n', '--name', help='file name')
# parser.add_argument('-l', '--list', help = 'list of nodes to trim from the graph')
#
# args = parser.parse_args()


def get_RNAseq_csv(csvpath):
    data = pd.read_csv(csvpath, index_col= 0, sep= "\t", header = 0)
    if len(data.columns) == 0:
        data = pd.read_csv(csvpath, index_col=0, sep=",", header=0)
    return data


def get_node_list(df_path):
    df = pd.read_excel(df_path, index_col=0, engine= 'openpyxl')
    nodelist = df['Metabolite Id']
    return nodelist


def xml_doc_parsing(path):
    """

    :param path: xml species/reaction  list pathway
    :return:
    """

    xml = mdom.parse(path)
    xml_nodes = nodes = xml.getElementsByTagName('species')
    xml_edges = xml.getElementsByTagName('reaction')

    return xml_nodes, xml_edges


def set_metabo_dict(xml_nodes):
    metabo_dict = {}
    for s in xml_nodes:
        attr = s.attributes
        species = dict(attr.items())
        n = species.pop('id')
        metabo_dict[n] = species

    return metabo_dict


def extract_reactant_products(xml_listof):
    mol_stack = []
    if xml_listof is not None:
        for i, mol in enumerate(xml_listof):
            tmp = dict(mol.attributes.items())
            if len(tmp) == 0:
                for child in mol.childNodes:
                    tmp = dict(child.attributes.items())
                    mol_stack.append(tmp.get('species'))
            else:
                mol_stack.append(tmp.get('species'))

    return mol_stack


def set_reaction_dict(xml_edges):
    reaction_dict = {}
    for r in xml_edges:
        ids = dict(r.attributes.items())
        notes = r.getElementsByTagName('notes')
        listOfReactants = r.getElementsByTagName('listOfReactants')[0].childNodes
        listOfProducts = r.getElementsByTagName('listOfProducts')
        n = ids.pop('id')
        ids.pop('name')

        reactants = extract_reactant_products(listOfReactants)
        products = extract_reactant_products(listOfProducts)
        supp = notes[0].lastChild.childNodes[0].childNodes[0].data

        edge_list = []

        for o in reactants:
            for t in products:
                w = (o, t)
                edge_list.append(w)

            reaction_dict[n + '_' + ids.get('metaid')] = {'id': ids, 'reaction': edge_list, 'notes': supp}

    return reaction_dict


def load_graph(path, csvpath, weights = None, weighting = True):
    """

    :param path: path to xml file
    :return: a graph in networkx format
    """
    xml_nodes, xml_edges = xml_doc_parsing(path)
    nodes = set_metabo_dict(xml_nodes)
    edges = set_reaction_dict(xml_edges)
    if weights == None :
        w_edges, gene_list = extract_RNA_weights(edges, csvpath)
    else :
        w_edges = pd.read_csv(weights, sep = ',')



    G = nx.Graph()
    G.add_nodes_from(nodes)
    for k in edges.keys():
        raw_metaid = edges.get(k).get('id').get('metaid')
        metaid = raw_metaid.split("_", maxsplit=2)[2]
        if w_edges[w_edges['reactions'] == metaid]['weight'].values[0] != 0 :
            for e in edges.get(k).get('reaction'):
                G.add_edge(*e)
                if weighting == True:
                    G.edges[e[0],e[1]]['weight'] = w_edges[w_edges['reactions'] == metaid]['weight'].values[0]


    # plt.subplot()
    # nx.draw(G, node_size = 1)
    # plt.show()

    dict_graph = {'nodes': nodes, 'edges': edges}

    return G , gene_list


def graph_trimming(g: nx.Graph, nodelist=None):
    nodelist =  ['M_'+val for val in nodelist]
    print("Before trimming : ",len(g), " nodes and ", len(g.edges), " edges")
    if nodelist is not None:
        g.remove_nodes_from(nodelist)

    g.remove_nodes_from(list(nx.isolates(g)))
    print("After trimming : ",len(g), " nodes and ", len(g.edges), " edges")
    return g


def load_RNA_data(dfpath):
    """

    :param dfpath: pathway to dataframe (csv file)
    :return: RNA-seq dataframe
    """
    rna_data = pd.read_csv(dfpath, sep='\t')

    return rna_data


def getlogFC(gene_id, file):
    row = file[file[utils.gene_id_col] == gene_id]

    log2FC = None
    if not row.empty:
        log2FC = row[utils.logFCcol].values

    return log2FC


def extract_RNA_weights(edges, csvpath):
    weight_table = pd.DataFrame()
    weight_col, reaction_col, = list(), list()
    gene_list = list()
    rna = get_RNAseq_csv(csvpath)
    hma.tempting_connection()
    print("Weight Attribution")
    for r in tqdm(edges.keys(), desc="reaction", total= len(edges.keys())):
        raw_metaid = edges.get(r).get('id').get('metaid')
        metaid = raw_metaid.split("_", maxsplit = 2)[2]
        reaction_col.append(metaid)
        try :
            s = hma.automatic_request_to_MA("reactions", metaid, "HumanGem")
            as_genes = hma.get_ensembl_geneid_list(hma.get_genes(s))
            weight_list = list()
            weight = 0

            for tup in as_genes:
                if tup[0] not in gene_list:
                    gene_list.append(tup[0])
                log2FC = getlogFC(tup[0], rna)
                print(tup[0], log2FC)
                if log2FC is not None:
                    weight_list.append(float(pow(2,log2FC)))

            if len(weight_list) > 0 :
                weight = np.mean(weight_list)
            weight_col.append(weight)

        except :
            print("!! WARNING : An error occured during the weight attribution !!")
            weight_col.append(0)
    weight_table['reactions'] = reaction_col
    weight_table['weight'] = weight_col

    weight_table.to_csv(utils.out_weight_table)
    return weight_table ,gene_list


def saveNxtogml(g:nx.Graph, out):
    nx.readwrite.generate_gmlexternal_(g, out)


def genesFromMulti(path):
    all = pd.read_csv(path, header =0, sep = '\t' )
    trim = all.drop_duplicates('external_gene_name',keep='first')

    trim.to_csv(path)


def load_global_graph(path, csvpath, weights = None):

    give_weight = False
    common_graph = load_graph(path, csvpath, weights, give_weight)
    return common_graph

def metabolic_weights(path, ids_table):
    table = get_RNAseq_csv(path)
    w_list = table[utils.metabo_column]
    for m in table.index:
        print(m)
        m_pat = re.compile(r'.*'+m+'.*')
        for m_node in ids_table['metabolite']:
            match = re.match(m_pat, m_node)
            if match is not None:
                print(m,': ',match.string)






if __name__ == "__main__":
    # ids = get_RNAseq_csv("/run/media/aurelien/ACOFFE/Stage/integration_job/clusters/AB_clusters_weights.csv")
    # metabolic_weights(utils.metabo_path, ids)
    nl = get_node_list("/run/media/aurelien/ACOFFE/Stage/integration_job/12859_2020_3564_MOESM1_ESM_modified.xlsx")

    g, genes = load_graph(path = "/run/media/aurelien/ACOFFE/Stage/integration_job/data_HMA/"
                   "U-251MG.xml-91080a939b86d903928cb7e2c321c2ff/U-251 MG.xml",
                   csvpath= "/run/media/aurelien/ACOFFE/Stage/integration_job/new_result_part3/new_result_part3/tables_deg/"
                            "control_hypoxie_ldha_ldhb_hypoxie_DEG_significant.tsv")

    G = graph_trimming(g, nl)
    with open('/run/media/aurelien/ACOFFE/Stage/integration_job/gene_list.txt', 'w') as genelist :
        for i in genes :
            if i is not None:
                genelist.write(i + '\r')
        genelist.close()

    # nx.write_gml(G,"/run/media/aurelien/ACOFFE/Stage/integration_job/G.gml")
    # Ge = ig.Graph.from_networkx(G)
    # Ge.write_graphml('/run/media/aurelien/ACOFFE/Stage/integration_job/common.graphml')

    # p_g = load_graph(args.xml)
    # g = graph_trimming(g, args.trim)
    # G = ig.Graph.from_networkx(g)
    #
    # output = args.output + args.name
    # G.write_graphml(output)

    # nodes, edges = xml_doc_parsing("/run/media/aurelien/ACOFFE/Stage/integration_job/data_HMA/"
    #                                "U-251MG.xml-91080a939b86d903928cb7e2c321c2ff/U-251 MG.xml")

    # extract_RNA_weights("/run/media/aurelien/ACOFFE/Stage/integration_job/data_HMA/"
    #                      "U-251MG.xml-91080a939b86d903928cb7e2c321c2ff/U-251 MG.xml",
    #                      "/run/media/aurelien/ACOFFE/Stage/integration_job/new_result_part3/new_result_part3/tables_deg/"
    #                      "control_hypoxie_ldha_ldhb_hypoxie_DEG_significant.tsv")
    #

    #
    # t = load_graph("/run/media/aurelien/ACOFFE/Stage/integration_job/data_HMA/"
    #                "U-251MG.xml-91080a939b86d903928cb7e2c321c2ff/U-251 MG.xml")
    # nl = get_node_list("/run/media/aurelien/ACOFFE/Stage/integration_job/12859_2020_3564_MOESM1_ESM_modified.xlsx")
    # g = graph_trimming(t, nl)
    # G = ig.Graph.from_networkx(g)
    #
    # G.write_graphml('/run/media/aurelien/ACOFFE/Stage/integration_job/Metabo_brain_graph.graphml')
    # plt.subplot()
    # nx.draw(g, node_size=1)
    # plt.show()

    # with open('/run/media/aurelien/ACOFFE/Stage/integration_job/Metabo_brain_graph.json', 'w') as out:
    #     out.write(json.dumps(json_graph.node_link_data(g)))
