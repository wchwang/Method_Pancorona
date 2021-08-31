# Created by woochanghwang at 08/07/2021

import pandas as pd
from toolbox.data_handler import save_obj, load_obj
import networkx as nx

def id_mapping(protein, mapping_dict):
    return mapping_dict[protein]

def save_graph(threshold):
    # threshold = 700
    human_interaction_df = load_obj(f"/Users/woochanghwang/PycharmProjects/MTIProject/General/data/STRING/v11/human_id_symbol_interaction_{threshold}")
    print(list(human_interaction_df))
    human_Graph = nx.from_pandas_edgelist(human_interaction_df,source='protein1_name',target='protein2_name')
    print(nx.info(human_Graph))
    save_obj(human_Graph,f"../Data/Network/human_string_symbol_graph_{threshold}")

def save_interactions(threshold):
    mouse_interaction_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/STRING/v11/9606.protein.links.v11.0.txt"
    network_info_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/STRING/v11/9606.protein.info.v11.0.txt"

    mouse_interaction_df = pd.read_csv(mouse_interaction_addr, sep=' ')
    network_info_df = pd.read_csv(network_info_addr, sep='\t')

    print(mouse_interaction_df.head())
    # print(list(network_info_df))
    network_info_df = network_info_df[['protein_external_id','preferred_name']]
    network_info_df = network_info_df.set_index('protein_external_id')
    # network_info_df = network_info_df[['preferred_name']]
    network_info_df.transpose()
    # print(list(network_info_df))
    protein_gene_map_dict = network_info_df['preferred_name'].to_dict()
    # threshold = 700
    mouse_interaction_df = mouse_interaction_df[mouse_interaction_df['combined_score']>=threshold]
    mouse_interaction_df['protein1_name'] = mouse_interaction_df.apply(lambda x:id_mapping(x['protein1'],protein_gene_map_dict), axis=1)
    mouse_interaction_df['protein2_name'] = mouse_interaction_df.apply(lambda x:id_mapping(x['protein2'],protein_gene_map_dict), axis=1)
    print(mouse_interaction_df)
    save_obj(mouse_interaction_df,f"/Users/woochanghwang/PycharmProjects/MTIProject/General/data/STRING/v11/human_id_symbol_interaction_{threshold}")

def main():
    threshold = 700
    save_interactions(threshold)
    save_graph(threshold)


if __name__ == '__main__':
    main()
