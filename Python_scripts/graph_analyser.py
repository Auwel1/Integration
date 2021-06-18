import igraph as ig
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
import networkx.algorithms.community as com
import community as lou
import numpy as np
import pandas as pd
import HMA_connect as hma
import utils as ut
from tqdm import tqdm
from sklearn.metrics import silhouette_samples, silhouette_score
import warnings
import load_graph as lg


warnings.filterwarnings('ignore')

def load_graph(graph_path):
    g = ig.load(graph_path)

    return g


def nx_load_graph(graph_path):
    g = nx.readwrite.read_gml(graph_path)
    return g


def plot_Igraph(g: ig.Graph):
    plt.subplot()
    ig.drawing.plot(g)
    plt.show()


def plot_nxGraph(g: nx.Graph, size):
    plt.subplot()
    nx.draw(g, node_size=size, with_labels=False, )
    plt.show()


def largest_component(g: nx.Graph):
    cc = nx.connected_components(g)
    lgcc = max(cc, key=len)
    lgcc_graph = g.subgraph(lgcc)

    return lgcc_graph


def spectral_clustering(g: nx.Graph):
    adj_matrix = nx.to_numpy_matrix(g)
    nodelist = list(g.nodes())
    scores = list()
    for n in range(10,48,1):
            mc = SpectralClustering(n_clusters=n, affinity= 'nearest_neighbors',
                                    assign_labels= 'discretize').fit_predict(adj_matrix)
            silhouette_avg = silhouette_score(adj_matrix, mc)
            print("For n_clusters =", n,
          "The average silhouette_score is :", silhouette_avg)
            scores.append(silhouette_avg)
    val = max(scores)
    n = scores.index(val)+10
    print("k = ", n , " with a score of ", max(scores))
    mc = SpectralClustering(n_clusters=n, affinity='nearest_neighbors',
                            assign_labels='discretize').fit_predict(adj_matrix)

    mc = np.transpose(pd.DataFrame(mc))
    mc.columns = nodelist
    return mc


def get_indicator_matrix(memb_cluster):
    x = list(range(0,max(memb_cluster.values[0])+1))
    indic_matrix = pd.DataFrame( columns= x, index= memb_cluster.columns)
    for c in indic_matrix.columns:
        for m in memb_cluster:
            if int(memb_cluster[m].values[0]) == c:
                indic_matrix.loc[m][c] = 1
            else :
                indic_matrix.loc[m][c] = 0

    return indic_matrix


def summarize_graph(g:nx.Graph, indicator_matrix):
    adj_matrix = nx.to_numpy_matrix(g)
    t_ind_matrix = np.transpose(indicator_matrix)
    summ_adj_matrix = t_ind_matrix.dot(adj_matrix.dot(indicator_matrix))
    summ_graph = nx.from_numpy_matrix(summ_adj_matrix.to_numpy(dtype= 'float'))

    return summ_graph


def get_cluster_members(indicator_matrix):
    id_per_cluster = pd.DataFrame()
    for cluster in indicator_matrix:
        id_list = indicator_matrix[indicator_matrix[cluster] == 1].index
        id_per_cluster = pd.concat([id_per_cluster, pd.Series(id_list)], ignore_index=True, axis=1)

    return id_per_cluster


def metabofromclustHMA(ids_cluster, ids_dataframe= None):
    hma.tempting_connection()

    ids_dat = pd.DataFrame()
    metabo_list = list()
    KEGG_list = list()
    MetaNetX_list = list()
    cluster_list = list()
    id_list = list()
    print('working on clusters')
    for cluster in ids_cluster:
        for id in tqdm(ids_cluster[cluster], desc = 'in cluster '+str(cluster), total= len(ids_cluster[cluster])) :

            if type(id) is str:
                req = id.split(sep = '_')[1]
                if ids_dataframe == None:
                    try:
                        search = hma.automatic_request_to_MA(request_type= 'metabolites',
                                                request= req, model= ut.HMA_model)
                        metabolite = search.get('name')
                    except:
                        print('!!WARNING : metabolite ' +req+ ' seems not be found !!')
                        continue
                    cluster_list.append("cluster " + str(cluster))
                    KEGG_id = ""
                    Meta_id = ""

                    if search.get('externalDbs').get('KEGG') is not None:
                        KEGG_id =  search.get('externalDbs').get('KEGG')[0].get('id')
                    if search.get('externalDbs').get('MetaNetX') is not None :
                        Meta_id = search.get('externalDbs').get('MetaNetX')[0].get('id')
                else:
                    ids_doc = pd.read_csv(ids_dataframe, sep = ',', header= 0)
                    row = ids_doc.loc[ids_doc['id'] == req]
                    if req not in ids_doc['id'].values:
                        continue
                    KEGG_id = ""
                    Meta_id = ""

                    metabolite = row['metabolite'].values[0]
                    if row['KEGG'].values[0] is not None:
                        KEGG_id =  row['KEGG'].values[0]
                    if row['MetaNetX'] is not None :
                        Meta_id = row['MetaNetX'].values[0]
                cluster_list.append("cluster " + str(cluster))
                metabo_list.append(metabolite)
                KEGG_list.append(KEGG_id)
                MetaNetX_list.append(Meta_id)
                id_list.append(id)

    ids_dat['cluster'] = cluster_list
    ids_dat['id'] = id_list
    ids_dat['metabolite'] = metabo_list
    ids_dat['KEGG'] = KEGG_list
    ids_dat['MetaNetX'] = MetaNetX_list
    return ids_dat

def ids_from_xml(xmldoc):
    species, reactions = lg.xml_doc_parsing(xmldoc)
    ids_dat = pd.DataFrame()
    metabo_list = list()
    KEGG_list = list()
    MetaNetX_list = list()
    id_list = list()
    hma.tempting_connection()
    for s in tqdm(species, desc = "species data loading", total = len(species)) :
        id = s.attributes.items()[1][1]
        spl = id.split('_')[1]
        try:
            search = hma.automatic_request_to_MA(request_type='metabolites',
                                                 request=spl, model=ut.HMA_model)
            metabolite = search.get('name')
        except:
            print('!!WARNING : metabolite ' + spl + ' seems not be found !!')
            continue

        id_list.append(id)
        metabo_list.append(metabolite)
        KEGG_id = ""
        Meta_id = ""

        if search.get('externalDbs').get('KEGG') is not None:
            KEGG_id = search.get('externalDbs').get('KEGG')[0].get('id')
        if search.get('externalDbs').get('MetaNetX') is not None:
            Meta_id = search.get('externalDbs').get('MetaNetX')[0].get('id')

        KEGG_list.append(KEGG_id)
        MetaNetX_list.append(Meta_id)

    ids_dat['id'] = id_list
    ids_dat['metabolite'] = metabo_list
    ids_dat['KEGG'] = KEGG_list
    ids_dat['MetaNetX'] = MetaNetX_list
    return ids_dat



if __name__ == '__main__':
    g = nx_load_graph("/run/media/aurelien/ACOFFE/Stage/integration_job/G.gml")
    g = g.to_undirected()
    lgcc = largest_component(g)

    mc = spectral_clustering(lgcc)
    summ_graph = summarize_graph(lgcc, get_indicator_matrix(mc))
    plot_nxGraph(summ_graph, 150)
    #
    ids = get_cluster_members(get_indicator_matrix(mc))
    frame = metabofromclustHMA(ids, "/run/media/aurelien/ACOFFE/Stage/integration_job/metabo_ids.csv")
    frame.to_csv("/run/media/aurelien/ACOFFE/Stage/integration_job/clusters/AB_clusters_noweights.csv")


    #  g = load_graph("/run/media/aurelien/ACOFFE/Stage/integration_job/Test.graphml")
