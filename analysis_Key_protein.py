# Created by woochanghwang at 18/07/2021


import pandas as pd
import toolbox.data_handler as dh
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from venn import venn
import toolbox.visual_utilities as vu

def draw_venn_diagram():

    virus_list = ['SARS-CoV', 'SARS-CoV-2']


    #######################
    # CoV , DIP vs DEP vs Hidden
    #######################
    for virus in virus_list:
        graph_addr = f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph"
        covid_graph = dh.load_obj(graph_addr)
        sip_nodes = covid_graph.nodes()
        key_protein = pd.read_csv(f"../result/{virus}/Sig.Genes/{virus}_key_protein_every.txt")['Gene'].tolist()
        DIP = pd.read_csv(f"../Data/DIP/{virus}_DIP.csv")['gene_name'].tolist()
        DEP = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()

        DIPnSIP = set(DIP)&set(sip_nodes)
        DEPnSIP = set(DEP)&set(sip_nodes)
        # hidden = set(sip_nodes) - DIPnSIP - DEPnSIP
        virus_dip_dep_dict = {
            f"{virus}_DIP":DIPnSIP,
            f"{virus}_DEP":DEPnSIP,
            f"{virus}_Key_Proteins":set(key_protein),
        }
        # venn(virus_dip_dep_dict)  # version 1
        if virus == 'SARS-CoV':
            virus_color = "#1f77b4"
        elif virus =="SARS-CoV-2":
            virus_color = "#ff7f0e"
        vu.draw_venn_3group(DIPnSIP,DEPnSIP,key_protein, group_labels=["DIP","DEP","Key Proteins"],
                            save_addr=f"../result/{virus}/{virus}_DIP_DEP_Hidden_Keyprotein_venn_v2.pdf",font_size=20,
                            set_colors= ('r','g',virus_color),
                            title=f"{virus}_DIP_DEP_Hidden")

        plt.savefig((f"../result/{virus}/{virus}_DIP_DEP_Hidden_Keyprotein_venn_v3.pdf"))
        plt.show()

def draw_venn_each_centrality():
    virus_list = ['SARS-CoV', 'SARS-CoV-2']

    #######################
    # CoV , DIP vs DEP vs Hidden
    #######################
    for virus in virus_list:
        network_anaysis_df = pd.read_csv(f"../result/{virus}/network_analysis/{virus}_A549_24h_centrality_RWR_result_pvalue.csv")

        eigen_list = network_anaysis_df[network_anaysis_df['Eigen_pvalue']< 0.01]['Gene'].tolist()
        degree_list = network_anaysis_df[network_anaysis_df['Degree_pvalue']< 0.01]['Gene'].tolist()
        bw_list = network_anaysis_df[network_anaysis_df['Between_plvaue']< 0.01]['Gene'].tolist()
        rwr_list = network_anaysis_df[network_anaysis_df['RWR_pvalue']< 0.01]['Gene'].tolist()


        virus_centrality_dict = {
            "Eigen":set(eigen_list),
            "Degree":set(degree_list),
            "Betweeneess":set(bw_list),
            "RWR":set(rwr_list)
        }
        venn(virus_centrality_dict)

        plt.savefig((f"../result/{virus}/{virus}_Centrality_venn.pdf"))
        plt.show()

def draw_venn_each_centraltiy_bw_corona():
    virus_list = ['SARS-CoV', 'SARS-CoV-2']


    centrality_list = ["Eigen", "Degree", "Between","RWR"]

    #######################
    # CoV , DIP vs DEP vs Hidden
    #######################
    for centrality in centrality_list:

        cov1_network_df = pd.read_csv(f"../result/SARS-CoV/network_analysis/SARS-CoV_A549_24h_centrality_RWR_result_pvalue.csv")
        cov2_network_df = pd.read_csv(f"../result/SARS-CoV-2/network_analysis/SARS-CoV-2_A549_24h_centrality_RWR_result_pvalue.csv")

        cov1_list = cov1_network_df[cov1_network_df[f"{centrality}_pvalue"]< 0.01]['Gene'].tolist()
        cov2_list = cov2_network_df[cov2_network_df[f"{centrality}_pvalue"]< 0.01]['Gene'].tolist()

        virus_centrality_dict = {
            "SARS-CoV" : set(cov1_list),
            "SARS-CoV-2" : set(cov2_list)
        }
        # venn(virus_centrality_dict)
        #
        # plt.savefig((f"../Figure/CoV1_CoV2_keyprotein_venn_{centrality}.pdf"))
        # plt.show()

        vu.draw_venn_2group(set(cov1_list), set(cov2_list), group_labels=['SARS','COVID-19'],font_size=20, set_colors= ("#1f77b4","#ff7f0e"),alpha=0.6,
                            save_addr=f"../Figure/CoV1_CoV2_keyprotein_venn_{centrality}_v4.pdf")

def draw_venn_key_proteins():
    virus_list = ['SARS-CoV', 'SARS-CoV-2']

    cov1_keyproteins = pd.read_csv("../result/SARS-CoV/Sig.Genes/SARS-CoV_key_protein_every.txt")['Gene'].tolist()
    cov2_keyproteins = pd.read_csv("../result/SARS-CoV-2/Sig.Genes/SARS-CoV-2_key_protein_every.txt")['Gene'].tolist()

    virus_centrality_dict = {
        "SARS-CoV" : set(cov1_keyproteins),
        "SARS-CoV-2" : set(cov2_keyproteins)
    }

    # venn(virus_centrality_dict)

    # plt.savefig((f"../Figure/CoV1_CoV2_keyprotein_venn.pdf"))
    # plt.show()
    vu.draw_venn_2group(set(cov1_keyproteins),set(cov2_keyproteins),group_labels=['SARS','COVID-19'], set_colors= ("#1f77b4","#ff7f0e"),alpha=0.6,
                        save_addr=f"../Figure/CoV1_CoV2_keyprotein_venn_v4.pdf", font_size=20)
    # #######################
    # # CoV , DIP vs DEP vs Hidden
    # #######################
    # for virus in virus_list:
    #     network_anaysis_df = pd.read_csv(f"../result/{virus}/network_analysis/{virus}_A549_24h_centrality_RWR_result_pvalue.csv")
    #
    #     eigen_list = network_anaysis_df[network_anaysis_df['Eigen_pvalue']< 0.01]['Gene'].tolist()
    #     degree_list = network_anaysis_df[network_anaysis_df['Degree_pvalue']< 0.01]['Gene'].tolist()
    #     bw_list = network_anaysis_df[network_anaysis_df['Between_plvaue']< 0.01]['Gene'].tolist()
    #     rwr_list = network_anaysis_df[network_anaysis_df['RWR_pvalue']< 0.01]['Gene'].tolist()
    #
    #
    #     virus_centrality_dict = {
    #         "Eigen":set(eigen_list),
    #         "Degree":set(degree_list),
    #         "Betweeneess":set(bw_list),
    #         "RWR":set(rwr_list)
    #     }
    #     venn(virus_centrality_dict)
    #
    #     plt.savefig((f"../result/{virus}/{virus}_Centrality_venn.pdf"))
    #     plt.show()


def overlap_similarlity(set1, set2):
    set1 = set(set1)
    set2 = set(set2)
    return len(set1&set2) / min(len(set1), len(set2))


def draw_bar_graph_keyprotein_functions():
    # centrality_list = ["every","eigen","degree","bw","rwr"]
    centrality_list = ["every"]

    similarity_dict = dict()

    for centrality in centrality_list:
    # centrality = centrality_list[3]
        cov1_enriched_result_df = pd.read_csv(f"../result/SARS-CoV/Sig.Genes/gProfiler_hsapiens_keyprotein_{centrality}.csv")
        cov2_enriched_result_df = pd.read_csv(f"../result/SARS-CoV-2/Sig.Genes/gProfiler_hsapiens_keyprotein_{centrality}_cov2.csv")

        cov1_enriched_result_df = cov1_enriched_result_df[cov1_enriched_result_df['source']=='GO:BP']
        cov2_enriched_result_df = cov2_enriched_result_df[cov2_enriched_result_df['source']=='GO:BP']

        cov1_enriched_selected_df = cov1_enriched_result_df[['term_name','negative_log10_of_adjusted_p_value']]
        cov2_enriched_selected_df = cov2_enriched_result_df[['term_name','negative_log10_of_adjusted_p_value']]

        cov1_bp_terms = cov1_enriched_selected_df['term_name'].tolist()
        cov2_bp_terms = cov2_enriched_selected_df['term_name'].tolist()

        cov1_cov2_similarity = overlap_similarlity(cov1_bp_terms, cov2_bp_terms)
        similarity_dict[centrality] = cov1_cov2_similarity

        cov2_enriched_selected_df = cov2_enriched_result_df.head(20)

        # cov1_enriched_selected_df['virus'] = "SARS-CoV"
        # cov2_enriched_selected_df['virus'] = "SARS-CoV-2"

        cov1_enriched_selected_df['virus'] = "SARS"
        cov2_enriched_selected_df['virus'] = "COVID-19"

        path_terms = cov2_enriched_selected_df['term_name'].tolist()
        cov1_enriched_selected_df = cov1_enriched_selected_df[cov1_enriched_selected_df['term_name'].isin(path_terms)]

        cov1_cov2_enriched_selected_df = pd.concat([cov1_enriched_selected_df, cov2_enriched_selected_df])

        # cov1_cov2_enriched_selected_df = pd.merge(left=cov2_enriched_selected_df,
        #                                           right = cov1_enriched_selected_df,
        #                                           how='left',
        #                                           left_on='term_name',
        #                                           right_on='term_name')

        # print(cov1_cov2_enriched_selected_df)
        print(cov2_enriched_selected_df)
        print(len(path_terms))
        print(cov1_enriched_selected_df)
        print(cov1_cov2_enriched_selected_df)
        plt.subplots(figsize=(8, 10))
        ax = sns.barplot(x='negative_log10_of_adjusted_p_value', y='term_name', hue='virus', data=cov1_cov2_enriched_selected_df,
                         saturation=0.8)
        plt.xlabel('-log(adjPvalue)')
        plt.ylabel('GO:BP')
        plt.legend(loc ='lower right')
        plt.tight_layout()

        plt.savefig(f"../Figure/CoV1_CoV2_keyprotein_{centrality}_v2.pdf")
        plt.show()

    print(similarity_dict)
    #
    # cov_similairty_df = pd.DataFrame.from_dict(similarity_dict,orient='index')
    # cov_similairty_df= cov_similairty_df.reset_index()
    # cov_similairty_df = cov_similairty_df.rename(columns={
    #     'index' : 'centrality',
    #     0 : 'similarity'
    # })
    #
    # print(cov_similairty_df)
    # cov_similairty_df.to_excel("../result/Cov1_Cov2_keyproteins_BP_similarity.xlsx")
    # ax = sns.barplot(x='centrality',y='similarity',data=cov_similairty_df)
    # # plt.savefig(f"../Figure/CoV1_CoV2_keyprotein_BP_similarity.pdf")
    # plt.show()





def main():
    # draw_venn_diagram()
    # draw_venn_each_centrality()
    # draw_bar_graph_keyprotein_functions()
    draw_venn_key_proteins()
    draw_venn_each_centraltiy_bw_corona()
if __name__ == '__main__':
    main()
