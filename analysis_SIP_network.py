# Created by woochanghwang at 17/07/2021

import pandas as pd
import toolbox.data_handler as dh
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from venn import venn

def draw_venn_diagram():

    virus_list = ['SARS-CoV', 'SARS-CoV-2']


    #######################
    # CoV , DIP vs DEP
    #######################
    for virus in virus_list:
        graph_addr = f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph"
        covid_graph = dh.load_obj(graph_addr)
        sip_nodes = covid_graph.nodes()
        DIP = pd.read_csv(f"../Data/DIP/{virus}_DIP.csv")['gene_name'].tolist()
        DEP = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()

        DIPnSIP = set(DIP)&set(sip_nodes)
        DEPnSIP = set(DEP)&set(sip_nodes)
        hidden = set(sip_nodes) - DIPnSIP - DEPnSIP
        virus_dip_dep_dict = {
            f"{virus}_DIP":DIPnSIP,
            f"{virus}_DEP":DEPnSIP,
            f"{virus}_Hidden":hidden,
        }
        venn(virus_dip_dep_dict)

        plt.savefig((f"../result/{virus}/{virus}_DIP_DEP_Hidde_venn.pdf"))
        plt.show()



def save_dip_dep_hidden():
    virus_list = ['SARS-CoV', 'SARS-CoV-2']


    #######################
    # CoV , DIP vs DEP
    #######################
    for virus in virus_list:
        graph_addr = f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph"
        covid_graph = dh.load_obj(graph_addr)
        sip_nodes = covid_graph.nodes()
        DIP = pd.read_csv(f"../Data/DIP/{virus}_DIP.csv")['gene_name'].tolist()
        DEP = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()

        DIPnSIP = set(DIP)&set(sip_nodes)
        DEPnSIP = set(DEP)&set(sip_nodes)
        hidden = set(sip_nodes) - DIPnSIP - DEPnSIP
        hidden_df = pd.DataFrame(list(hidden), columns=['Gene_name'])
        hidden_df.to_csv(f"../result/{virus}/{virus}_Hidden.csv", index=False)

def draw_venn_hidden():
    virus_list = ['SARS-CoV', 'SARS-CoV-2']
    #######################
    # DIP CoV vs CoV-2
    ######################
    dip_dict = dict()
    for virus in virus_list:
        DIP = pd.read_csv(f"../result/{virus}/Hidden/{virus}_Hidden.csv")['Gene_name'].tolist()
        dip_dict[f"{virus}_DIP"] = set(DIP)

    venn(dip_dict)
    plt.savefig((f"../result/{virus_list[0]}/Hidden/Hidden_virus_venn.pdf"))
    plt.show()

def main():


    # draw_venn_diagram()
    # save_dip_dep_hidden()
    draw_venn_hidden()

if __name__ == '__main__':
    main()
