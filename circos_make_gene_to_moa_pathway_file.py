# Created by woochanghwang at 29/05/2020
'''
Modified by Woochang Hwang at 21/07/2020
'''
import pandas as pd
import toolbox.data_handler as dh
from collections import OrderedDict


def main():
    '''
    Reactome enriched result of HIDDEN proteins
    Pathways
    - Metabolism of RNA
    - Immune System
    - Cell Cycle
    :return:
    '''
    # moa_to_pathway_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/MoA/200drug/MoA_Reactome_pathways_lowlevel_SOM_Ref_v1.xlsx"
    reactome_gmt_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/Reactome/ReactomePathways_for_gProfiler.gmt"

    reactome_gene_dict = dict()

    with open(reactome_gmt_addr) as reactome_gmt_f:
        reactome_path_gene = [x.strip().split('\t') for x in reactome_gmt_f.readlines()]

    for path in reactome_path_gene:
        reactome_gene_dict[path[1]] = path[2:]


    moa_to_gene_dict = dict()
    hidden_paths = ['Metabolism of RNA','Immune System', 'Cell Cycle','Metabolism of proteins']
    path_labels = ['Met.RNA','Immune','CellCycle','Met.Proteins']
    for label,path in zip(path_labels, hidden_paths):
        print(label, path)
        path_genes = reactome_gene_dict[path]

        moa_to_gene_dict[label] = list(set(path_genes))


    for key, values in moa_to_gene_dict.items():
        print(key,len(values))


    dh.save_obj(moa_to_gene_dict, "../Data/hidden_to_keyPaths_dict")

if __name__ == '__main__':

    main()