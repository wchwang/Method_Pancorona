# Created by woochanghwang at 08/07/2021

import pandas as pd
from toolbox.data_handler import save_obj, load_obj
import networkx as nx
import toolbox.visual_utilities as vu

def is_in_STRING(gene, graph):
    if gene in graph:
        return "IN"
    else:
        return "Not_IN"

def make_DIP_DEP():
    human_Graph = load_obj(f"../Data/Network/human_string_symbol_graph_700")
    proteins_in_Graph = human_Graph.nodes()
    print(len(proteins_in_Graph))

    DEP_both_addr = "../Data/DEP/DEP_two_virus.xlsx"
    DIP_both_addr = "../Data/DIP/DIP_two_virus.xlsx"
    virus_list = ["SARS-CoV","SARS-CoV-2"]

    #####################
    ## Check DIP or DEP are in STRING or not
    ######################
    # for virus in virus_list:
    #     dip_df = pd.read_excel(DIP_both_addr,sheet_name=virus)
    #     # dip_genes = dip_df['gene_name'].tolist()
    #     # print(len(dip_genes))
    #     # print(len(set(dip_genes)&set(proteins_in_Graph)))
    #     dip_df['In_STRING'] = dip_df.apply(lambda row: is_in_STRING(row['gene_name'],proteins_in_Graph),axis=1)
    #     dip_in_STRING_df = dip_df[dip_df['In_STRING']=='IN']
    #     dip_in_STRING_df = dip_in_STRING_df.drop(columns='In_STRING')
    #     dip_in_STRING_df.to_csv(f"../Data/DIP/{virus}_DIP.csv",index=False)

    for virus in virus_list:
        dep_df = pd.read_excel(DEP_both_addr,sheet_name=virus)
        print(dep_df)
        dep_genes = dep_df['Gene_name'].tolist()
        print(len(dep_genes))
        print(len(set(dep_genes)&set(proteins_in_Graph)))
        dep_df['In_STRING'] = dep_df.apply(lambda row: is_in_STRING(row['Gene_name'],proteins_in_Graph),axis=1)
        dep_in_STRING_df = dep_df[dep_df['In_STRING']=='IN']
        dep_in_STRING_df = dep_in_STRING_df.drop(columns='In_STRING')
        dep_in_STRING_df.to_csv(f"../Data/DEP/{virus}_DEP.csv",index=False)

def draw_venn_diagram():
    from venn import venn
    import matplotlib.pyplot as plt

    virus_list = ["SARS-CoV","SARS-CoV-2"]

    ########################
    ## CoV , DIP vs DEP
    ########################
    # for virus in virus_list:
    #     DIP = pd.read_csv(f"../Data/DIP/{virus}_DIP.csv")['gene_name'].tolist()
    #     DEP = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()
    #     virus_dip_dep_dict = {
    #         f"{virus}_DIP":set(DIP),
    #         f"{virus}_DEP":set(DEP)
    #     }
    #     venn(virus_dip_dep_dict)
    #
    #     plt.savefig((f"../result/{virus}/{virus}_DIP_DEP_venn.pdf"))
    #     plt.show()

    ########################
    ## DIP CoV vs CoV-2
    #######################
    # dip_dict = dict()
    # for virus in virus_list:
    #     DIP = pd.read_csv(f"../Data/DIP/{virus}_DIP.csv")['gene_name'].tolist()
    #     dip_dict[f"{virus}_DIP"] = set(DIP)
    #
    # venn(dip_dict)
    # plt.savefig((f"../result/DIP_DEP/DIP_virus_venn.pdf"))
    # plt.show()

    ########################
    ## DEP CoV vs CoV-2
    #######################
    dep_dict = dict()
    for virus in virus_list:
        DIP = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()
        dep_dict[f"{virus}_DEP"] = set(DIP)

    venn(dep_dict)
    plt.savefig((f"../result/DIP_DEP/DEP_virus_venn.pdf"))
    plt.show()

def draw_venn_diagram_v2():
    virus_list = ["SARS-CoV","SARS-CoV-2"]


    DIP_SARS = pd.read_csv(f"../Data/DIP/{virus_list[0]}_DIP.csv")['gene_name'].tolist()
    DIP_Cov2 = pd.read_csv(f"../Data/DIP/{virus_list[1]}_DIP.csv")['gene_name'].tolist()

    vu.draw_venn_2group(group1=set(DIP_SARS), group2=set(DIP_Cov2), group_labels=['SARS','COVID-19'],
                        save_addr=f"../result/DIP_DEP/DIP_virus_venn_v2.pdf", font_size=20)

    DEP_SARS = pd.read_csv(f"../Data/DEP/{virus_list[0]}_DEP.csv")['Gene_name'].tolist()
    DEP_Cov2 = pd.read_csv(f"../Data/DEP/{virus_list[1]}_DEP.csv")['Gene_name'].tolist()

    vu.draw_venn_2group(group1=set(DEP_SARS), group2=set(DEP_Cov2), group_labels=['SARS','COVID-19'],
                        save_addr=f"../result/DIP_DEP/DEP_virus_venn_v2.pdf", font_size=20)

    DIP_SARS = pd.read_csv(f"../Data/DIP/{virus_list[0]}_DIP.csv")['gene_name'].tolist()
    DIP_Cov2 = pd.read_csv(f"../Data/DIP/{virus_list[1]}_DIP.csv")['gene_name'].tolist()

    vu.draw_venn_2group(group1=set(DIP_SARS), group2=set(DEP_SARS), group_labels=['DIP','DEP'],
                        save_addr=f"../result/DIP_DEP/cov_cov2_dip_v2.pdf", font_size=20, title="SARS")

    DEP_SARS = pd.read_csv(f"../Data/DEP/{virus_list[0]}_DEP.csv")['Gene_name'].tolist()
    DEP_Cov2 = pd.read_csv(f"../Data/DEP/{virus_list[1]}_DEP.csv")['Gene_name'].tolist()

    vu.draw_venn_2group(group1=set(DIP_Cov2), group2=set(DEP_Cov2), group_labels=['DIP','DEP'],
                        save_addr=f"../result/DIP_DEP/cov_cov2_dep_v2.pdf", font_size=20, title="COVID-19")


def main():

    # make_DIP_DEP()
    # draw_venn_diagram()
    draw_venn_diagram_v2()


if __name__ == '__main__':
    main()
