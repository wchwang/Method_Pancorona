# Created by woochanghwang at 17/07/2021

# Created by woochanghwang at 17/06/2020

import pandas as pd
import toolbox.data_handler as dh
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns

'''
Calculate all shotest paths between DIP and DEP 
get ratio 
length 1
length 2
length 3
'''

def remove_DIP_DEP_self_edges(covid_graph, dip_in_graph, dep_in_graph):
    covid_graph_edges = covid_graph.edges()

    covid_graph_edges_no_DIP_DEP_self = []

    for edge in covid_graph_edges:
        if len(set(edge)&set(dip_in_graph))==2: continue
        elif len(set(edge)&set(dep_in_graph))==2: continue
        else:
            covid_graph_edges_no_DIP_DEP_self.append(edge)

    return nx.Graph(covid_graph_edges_no_DIP_DEP_self)

def find_pair_shortest_path_length(covid_graph, dip_proteins, dep_proteins, virus):

    covid_graph_nodes = list(covid_graph.nodes())

    dip_in_graph = list(set(dip_proteins)&set(covid_graph_nodes))
    dep_in_graph = list((set(dep_proteins)&set(covid_graph_nodes))-set(dip_in_graph))
    hidden_in_graph = list(set(covid_graph_nodes) - set(dip_in_graph) - set(dep_in_graph))

    print(len(dip_in_graph), len(dep_in_graph), len(hidden_in_graph), virus)
    dip_dep_pair_list = list(itertools.product(dip_in_graph, dep_in_graph))

    print(len(dip_dep_pair_list))
    # print(dip_dep_pair_list[:5])

    print(len(covid_graph.edges()))
    covid_graph = remove_DIP_DEP_self_edges(covid_graph, dip_in_graph, dep_in_graph)

    print(len(covid_graph.edges()))

    dip_dep_pair_path_length_list = []
    dip_dep_pair_no_path = []
    for pair in dip_dep_pair_list:
        pair = list(pair)
        pair_length = pair[:]
        try:
            length = nx.shortest_path_length(covid_graph,source=pair[0],target=pair[1])
            pair_length.append(length)
            dip_dep_pair_path_length_list.append(pair_length)

        except:
            print(pair)
            dip_dep_pair_no_path.append(pair)

    dip_dep_pair_path_length_df = pd.DataFrame(dip_dep_pair_path_length_list,columns=['Protein','Protein','Path_length'])
    dip_dep_pair_path_length_addr = f"../result/{virus}/network/{virus}_dip_to_dep_path_length_no_DIP_DEP_self.csv"
    dip_dep_pair_path_length_df.to_csv(dip_dep_pair_path_length_addr,index=False)

    dip_dep_pair_no_path_df = pd.DataFrame(dip_dep_pair_no_path,columns=['DIP','DEP'])
    dip_dep_pair_no_path_df.to_csv(f"../result/{virus}/network/{virus}_dip_to_dep_no_path.csv")

    ################################################

    # dip_dep_pair_path_length_df.plot.hist(grid=True, bins=20, rwidth=0.9,color='#607c8e')
    # dip_dep_pair_path_length_df.plot.hist(grid=True,alpha=0.5, color='#607c8e')
    # plt.title('COVID19 network DIP to DEP path length')
    # plt.xlabel('Length')
    # plt.ylabel('DIP,DEP pair')
    # # plt.grid(axis='y', alpha=0.75)
    # plt.show()

def draw_path_length_histogram():
    dip_dep_pair_path_length_SARS_CoV_addr = "../result/SARS-CoV/network/SARS-CoV_dip_to_dep_path_length_no_DIP_DEP_self.csv"
    dip_dep_pair_path_length_SARS_CoV_2_addr = "../result/SARS-CoV-2/network/SARS-CoV-2_dip_to_dep_path_length_no_DIP_DEP_self.csv"


    pair_path_length_SARS_CoV_df = pd.read_csv(dip_dep_pair_path_length_SARS_CoV_addr)
    pair_path_length_SARS_CoV_2_df = pd.read_csv(dip_dep_pair_path_length_SARS_CoV_2_addr)

    # length_6, pair_numbers_6 = zip(*Counter(pair_path_length_SARS_CoV_df['Path_length'].to_list()).items())
    # length_24, pair_numbers_24 = zip(*Counter(pair_path_length_SARS_CoV_2_df['Path_length'].to_list()).items())
    #
    # print(length_6,pair_numbers_6)

    pair_length_SARS_CoV_counter=  Counter(pair_path_length_SARS_CoV_df['Path_length'].to_list())
    pair_length_SARS_CoV_2_counter = Counter(pair_path_length_SARS_CoV_2_df['Path_length'].to_list())

    pair_length_SARS_CoV_counter_df = pd.DataFrame.from_dict(pair_length_SARS_CoV_counter,orient='index').reset_index()
    pair_length_SARS_CoV_counter_df = pair_length_SARS_CoV_counter_df.rename(columns={'index':'length', 0:'DIP-DEP_pairs'})
    pair_length_SARS_CoV_counter_df['network'] = 'SARS-CoV'

    total_pairs = pair_length_SARS_CoV_counter_df['DIP-DEP_pairs'].sum()
    pair_length_SARS_CoV_counter_df['percentage'] = pair_length_SARS_CoV_counter_df['DIP-DEP_pairs']/total_pairs * 100

    pair_length_SARS_CoV_2_counter_df = pd.DataFrame.from_dict(pair_length_SARS_CoV_2_counter, orient='index').reset_index()
    pair_length_SARS_CoV_2_counter_df = pair_length_SARS_CoV_2_counter_df.rename(columns={'index': 'length', 0: 'DIP-DEP_pairs'})
    pair_length_SARS_CoV_2_counter_df['network'] = 'SARS-CoV-2'
    total_pairs = pair_length_SARS_CoV_2_counter_df['DIP-DEP_pairs'].sum()
    pair_length_SARS_CoV_2_counter_df['percentage'] = pair_length_SARS_CoV_2_counter_df['DIP-DEP_pairs']/total_pairs * 100

    pair_length_counter_df = pd.concat([pair_length_SARS_CoV_counter_df, pair_length_SARS_CoV_2_counter_df])
    # pair_length_counter_df = pair_length_counter_df.sort_values(by='length',ascending=False)
    print(pair_length_counter_df)

    ## sns.catplot(x='length',y='DIP-DEP_pairs',hue='network', data=pair_length_counter_df, kind='bar')
    ax = sns.histplot(data=pair_length_counter_df,
                      x = 'network',
                      weights='percentage',
                      hue='length',
                      multiple='stack',
                      edgecolor='white',
                      # palette='tab20c',
                      palette = "Set3",
                      hue_order=[9,8,7,6,5,4,3,2,1],
                      shrink=0.8)
    # print(pair_length_SARS_CoV_counter)
    # n,bins,patches = plt.hist([pair_path_length_SARS_CoV_df['Path_length'],pair_path_length_SARS_CoV_2_df['Path_length']])
    # plt.grid()
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Virus')
    legend = ax.get_legend()
    legend.set_bbox_to_anchor((1, 1))
    plt.tight_layout()

    # plt.savefig("../result/covid19_dip_to_dep_path_barplot_no_DIP_DEP_self.png")
    # plt.savefig("../result/covid19_dip_to_dep_path_stack_barplot_no_DIP_DEP_self.pdf")
    plt.show()
    # sns.barplot()



def main():
    viruses = ['SARS-CoV', 'SARS-CoV-2']
    virus = viruses[1]
    graph_addr = f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph"
    covid_graph = dh.load_obj(graph_addr)


    dip_proteins = pd.read_csv(f"../Data/DIP/{virus}_DIP.csv")['gene_name'].tolist()
    dep_proteins = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()

    # find_pair_shortest_path_length(covid_graph, dip_proteins, dep_proteins, virus)

    draw_path_length_histogram()


if __name__ == '__main__':
    main()