import requests
import argparse
import pandas as pd
import json

parser = argparse.ArgumentParser()

parser.add_argument('-t', '--request_type', help='"reactions", "metabolites", "genes", "compartments", '
                                                 '"integrated models", "subsystems", "interaction-partners", ...')
parser.add_argument('-r', '--request', help='request to metabolic atlas')
parser.add_argument('-s', '--stack', help = 'request list')
parser.add_argument('-m', '--model', help='model organism ')

args = parser.parse_args()


def tempting_connection():
    print('Attempting to connect Metabolic Atlas API : ')
    response = requests.get('https://metabolicatlas.org/api/v2/#/')
    if response.status_code == 200:
        print('Connected with success')
    else:
        print('Trouble in connecting, see the error code : ', response.status_code)


def manual_request_to_MA():
    request = 'https://metabolicatlas.org/api/v2/' + args.request_type + '/' + args.request + '/' + '?model=' + args.model
    search = requests.get(request).json()

    return search


def automatic_request_to_MA(request_type, request, model):

    try :
        req = 'https://metabolicatlas.org/api/v2/' + request_type + '/' + request + '/' + '?model=' + model
        search = requests.get(req).json()
    except  json.JSONDecodeError :
        print("Request error, maybe ",request, " doesn't exist in HMA" )
        return

    return search


def get_genes(search):
    genes_dict = search.get('genes')
    return genes_dict


def get_ensembl_geneid_list(gene_list):
    id_list = list()
    for k in gene_list:
        id_list.append((k.get('name'), k.get('id')))

    return id_list

def geneid_from_list():
    tempting_connection()


    s = request_to_MA()
    get_ensembl_geneid_list(get_genes(s))



if __name__ == "__main__":
    tempting_connection()