# Created by woochanghwang at 16/07/2021
import pandas as pd
import numpy as np
import networkx as nx

def calc_F1_score(recall, precision):
    f1_score = 2 * ((precision * recall) / (precision + recall))
    return f1_score

def make_drug_patway():
    virus = "SARS-CoV"
    candidate_drug_df= pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.xlsx")
    candidate_drugs = candidate_drug_df['source'].tolist()
    print(candidate_drugs)

    candidate_drugs_enriched_pathways = pd.DataFrame(columns=['term_name'])

    for drug in candidate_drugs:
        enriched_pathway_addr = f"../result/{virus}/Drug/enriched_pathway/enriched_reactome_pathways_{drug}.csv"
        enriched_pathway_df = pd.read_csv(enriched_pathway_addr)
        print(enriched_pathway_df)
        if len(enriched_pathway_df) > 0:
            enriched_pathway_df[drug] = enriched_pathway_df.apply(lambda row: calc_F1_score(row['recall'],row['precision']),axis=1)
            enriched_pathway_df = enriched_pathway_df[['term_name',drug]]
            candidate_drugs_enriched_pathways = pd.merge(left=candidate_drugs_enriched_pathways,
                                                         right = enriched_pathway_df,
                                                         how='outer')

        else:
            candidate_drugs_enriched_pathways[drug]=np.nan

    candidate_drugs_enriched_pathways['Count'] = candidate_drugs_enriched_pathways.count(axis=1)

    candidate_drugs_enriched_pathways.to_excel("../result/SARS-CoV/Drug/candidate_drug_enriched_pathways.xlsx",index=False)
    print(candidate_drugs_enriched_pathways)

def get_whole_enriched_pathways(virus,drug_result_df, threshold):
    # threshold = 0.001
    effective_drugs = drug_result_df['source'].tolist()
    # pa_addr_pref = "../result/{}/Drug/enriched_pathway_extended_no_path/enriched_reactome_pathways_{}.csv"
    pa_addr_pref = "../result/{}/Drug/extended_all_enriched_pathway/enriched_reactome_pathways_{}.csv"
    whole_pa_list = []
    for drug in effective_drugs:
        pa_df = pd.read_csv(pa_addr_pref.format(virus,drug))
        if len(pa_df) > 0:
            pa_df = pa_df[pa_df['p_value'] < threshold]  #filter
            pa_list = pa_df['term_name'].to_list()
            if len(pa_list) ==0:
                print("no path drug: ", drug)
            whole_pa_list+=pa_list
        else:
            print(drug)

    whole_pa_list = list(set(whole_pa_list))
    print(len(whole_pa_list))

    return whole_pa_list

def make_reactome_hierachi():
    reactome_hierachi_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/Reactome/ReactomePathwaysRelation_HSA.txt"
    reactome_path_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/Reactome/ReactomePathways.txt"

    with open(reactome_hierachi_addr) as r_h_f:
        reactome_pair = [x.strip().split('\t') for x in r_h_f.readlines()]

    with open(reactome_path_addr) as r_p_f:
        reactome_path_list = [x.strip().split('\t') for x in r_p_f.readlines()]

    reactome_path_map_dict = dict()
    reactome_name_to_id_dict = dict()
    for path in reactome_path_list:
        if path[2] == 'Homo sapiens':
            reactome_path_map_dict[path[0]] = path[1]
            reactome_name_to_id_dict[path[1]] = path[0]
    reactome_diGraph = nx.DiGraph(reactome_pair)

    # print(reactome_diGraph.nodes)

    return reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict


def find_reactome_parents(path_name):
    reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict = make_reactome_hierachi()
    reactome_id = reactome_name_to_id_dict.get(path_name,'NA')
    if reactome_id != 'NA':
        ancestor_temr_id = list(nx.ancestors(reactome_diGraph,reactome_id))
        parents = [reactome_path_map_dict[id] for id in ancestor_temr_id]
        parents = ' | '.join(parents)
        return parents
    else:   return 'NA'


def find_reactome_hierachi_top2(path_name):
    reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict = make_reactome_hierachi()
    reactome_id = reactome_name_to_id_dict.get(path_name,'NA')
    # print(path_name, reactome_id)
    if reactome_id != 'NA':
        ancestor_temr_id = list(nx.ancestors(reactome_diGraph,reactome_id))
        # print(ancestor_temr_id, len(ancestor_temr_id))
        if len(ancestor_temr_id) == 1 :
            # print("return:", reactome_path_map_dict[reactome_id])
            return_path = reactome_path_map_dict[reactome_id]
            return return_path

        elif len(ancestor_temr_id) > 1:
            # predecessors = nx.predecessor(reactome_diGraph,reactome_id)
            # predecessor_id = reactome_name_to_id_dict.get(predecessors[0], 'NA')
            predecessors = list(reactome_diGraph.predecessors(reactome_id))
            # print("pre:",predecessors)
            # predecessor_id = reactome_name_to_id_dict.get(predecessors[0],'NA')
            predecessor_name = reactome_path_map_dict[predecessors[0]]
            return find_reactome_hierachi_top2(predecessor_name)
            # return 'NA'

        # parents = ' | '.join(parents)
    else:   return 'NA'

def find_reactome_hierachi_top3(path_name):
    reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict = make_reactome_hierachi()
    reactome_id = reactome_name_to_id_dict.get(path_name,'NA')
    # print(path_name, reactome_id)
    if reactome_id != 'NA':
        ancestor_temr_id = list(nx.ancestors(reactome_diGraph,reactome_id))
        # print(ancestor_temr_id, len(ancestor_temr_id))
        if len(ancestor_temr_id) == 1 :
            # print("return:", reactome_path_map_dict[reactome_id])
            return_path = reactome_path_map_dict[reactome_id]
            return return_path
        if len(ancestor_temr_id) == 2 :
            # print("return:", reactome_path_map_dict[reactome_id])
            return_path = reactome_path_map_dict[reactome_id]
            return return_path

        # elif len(ancestor_temr_id) > 1:
        #     # predecessors = nx.predecessor(reactome_diGraph,reactome_id)
        #     # predecessor_id = reactome_name_to_id_dict.get(predecessors[0], 'NA')
        #     predecessors = list(reactome_diGraph.predecessors(reactome_id))
        #     print("pre:",predecessors)
        #     # predecessor_id = reactome_name_to_id_dict.get(predecessors[0],'NA')
        #     if len(ancestor_temr_id) > 3:
        #         predecessor_name = reactome_path_map_dict[ancestor_temr_id[-3]]
        #     else:
        #         predecessor_name = reactome_path_map_dict[ancestor_temr_id[0]]
        #
        #     return predecessor_name
        elif len(ancestor_temr_id) > 2:
            # predecessors = nx.predecessor(reactome_diGraph,reactome_id)
            # predecessor_id = reactome_name_to_id_dict.get(predecessors[0], 'NA')
            predecessors = list(reactome_diGraph.predecessors(reactome_id))
            # print("pre:",predecessors)
            # predecessor_id = reactome_name_to_id_dict.get(predecessors[0],'NA')
            predecessor_name = reactome_path_map_dict[predecessors[0]]
            return find_reactome_hierachi_top3(predecessor_name)

        # parents = ' | '.join(parents)
    else:   return 'NA'

def find_ancestor_enriched_pathways(virus, whole_pathways, threshold):

    # whole_pathways=whole_pathways[:10]
    enriched_pathway_df= pd.DataFrame(whole_pathways,columns=['Pathways'])
    # print(enriched_pathway_df)
    # enriched_pathway_df = enriched_pathway_df[:4]
    enriched_pathway_df['Parents'] = enriched_pathway_df['Pathways'].apply(find_reactome_parents)
    enriched_pathway_df['Top2'] = enriched_pathway_df['Pathways'].apply(find_reactome_hierachi_top2)

    enriched_pathway_df['Top3'] = enriched_pathway_df['Pathways'].apply(find_reactome_hierachi_top3)
    print(enriched_pathway_df)
    # enriched_pathway_df.to_excel(
    #     f"../result/{virus}/Drug/{virus}_extendAll_Reactome_pathways_hierachi_{threshold}_top2.xlsx", index=False)
    enriched_pathway_df.to_excel(
        f"../result/{virus}/Drug/{virus}_extendAll_Reactome_pathways_hierachi_{threshold}_top3.xlsx", index=False)


def find_low_level_pathways(virus, moa_to_pathway_df, threshold, low_level_pathway_addr):
    # print(moa_to_pathway_df)
    top2_groupby = moa_to_pathway_df.groupby('Top2')['Pathways'].agg(list).reset_index(name='Pathways')
    top2_groupby_size = moa_to_pathway_df.groupby('Top2')['Pathways'].size().reset_index(name='counts')
    top2_pathways_size_df = top2_groupby.merge(top2_groupby_size)

    high_level_top2_df = top2_pathways_size_df[top2_pathways_size_df['counts']>1]
    high_level_top2 = high_level_top2_df['Top2'].to_list()
    high_level_top2 = [x.strip() for x in high_level_top2]
    # print(high_level_top2)

    # moa_to_pathway_df = moa_to_pathway_df[['Pathways','Parents','Top2','SOM','C10','SOM_MoA','Range']]
    # moa_to_pathway_df['Parents_counts'] = moa_to_pathway_df['Parents'].str.split('|').str.len()

    low_level_pathways_df = moa_to_pathway_df.dropna()
    # print(low_level_pathways_df)

    low_level_pathways_df = low_level_pathways_df[~low_level_pathways_df['Pathways'].isin(high_level_top2)]
    print(low_level_pathways_df)
    # low_level_pathways_df.to_excel(f"../result/{virus}/Drug/{virus}_extendAll_Reactome_low_level_{threshold}_top2.xlsx")
    low_level_pathways_df.to_excel(low_level_pathway_addr)

    return low_level_pathways_df['Pathways'].to_list()

def make_matrix(virus,candidate_drug_addr, whole_pa_list, threshold, fi_matrix_addr):

    # pathway_addr_prefix = "../result/{}/Drug/enriched_pathway_extended_no_path/enriched_reactome_pathways_{}.csv"
    pathway_addr_prefix = "../result/{}/Drug/extended_all_enriched_pathway/enriched_reactome_pathways_{}.csv"


    pathway = whole_pa_list
    print(len(pathway), pathway[:3])

    sorted_pathway = sorted(pathway)
    print(sorted_pathway[:3])

    candidate_drug_df = pd.read_excel(candidate_drug_addr)

    drug_names = candidate_drug_df['source'].tolist()

    # print(drug_names)
    print(sorted_pathway)
    print(len(sorted_pathway))
    drug_pathway_matrix = []

    for drug in drug_names:
        drug_pathway_enriched_df = pd.read_csv(pathway_addr_prefix.format(virus, drug))
        drug_pathway_values = [drug]
        # print(drug)
        # print(drug_pathway_enriched_df.head())
        for path in sorted_pathway:
            # path_df = drug_pathway_enriched_df[drug_pathway_enriched_df['term_name'] == path][['precision', 'recall']]
            path_df = drug_pathway_enriched_df[(drug_pathway_enriched_df['term_name'] == path) & (drug_pathway_enriched_df['p_value'] < 0.01)][['precision', 'recall']]
            # print(path_df)
            if len(path_df) > 0:
                recall = path_df.iloc[0]['recall']
                precision = path_df.iloc[0]['precision']
                f1_score = 2 * ((precision * recall) / (precision + recall))
                drug_pathway_values.append(f1_score)
            else:
                drug_pathway_values.append(0.0)
        drug_pathway_matrix.append(drug_pathway_values)

    drug_pathway_df = pd.DataFrame(drug_pathway_matrix, columns=['Drug_name'] + sorted_pathway)
    print(drug_pathway_df)

    # sns.clustermap(drug_pathway_df)
    #
    # plt.show()

    # drug_pathway_df.to_csv(
    #     f"../result/SARS-CoV/Drug/F1/{virus}_less5_drug_to_low_level_pathways_{threshold}_matrix.csv",
    #     index=False)
    drug_pathway_df.to_csv(fi_matrix_addr,index=False)

def make_matrix_merge_to_top3(virus,low_level_pathway_addr,f1_matrix_addr,threshold,f1_matrix_merged_addr):
    f1_matrix_df = pd.read_csv(f1_matrix_addr,index_col=0)
    f1_matrix_df_T = f1_matrix_df.T
    print(f1_matrix_df_T)

    low_level_pathway_df = pd.read_excel(low_level_pathway_addr)
    print(list(low_level_pathway_df))
    low_to_top3_df = low_level_pathway_df[['Pathways','Top3']]
    f1_matrix_top3_df = pd.merge(left = low_to_top3_df,
                                 right = f1_matrix_df_T,
                                 how='right',
                                 left_on='Pathways',
                                 right_index=True)
    f1_matrix_top3_df_T = f1_matrix_top3_df.drop(columns=['Pathways'])
    f1_matrix_top3_merged_df_T = f1_matrix_top3_df.groupby(['Top3']).mean()
    # f1_matrix_top3_merged_df = f1_matrix_top3_merged_df.reset_index()
    # f1_matrix_top3_merged_df = f1_matrix_top3_merged_df.rename(columns={
    #     'Top3' : 'Drugs'
    # })

    f1_matrix_top3_merged_df_T['count_drugs'] = (f1_matrix_top3_merged_df_T>0).sum(1)
    f1_matrix_top3_merged_df_T.to_csv(f"../result/{virus}/Drug/F1/{virus}_extendAll_drug_to_low_level_pathways_{threshold}_merged_top3_matrix_T.csv")

    f1_matrix_top3_merged_df_T = f1_matrix_top3_merged_df_T.drop(columns=['count_drugs'])
    f1_matrix_top3_merged_df = f1_matrix_top3_merged_df_T.T
    print(f1_matrix_top3_merged_df)
    candidate_drug_df = pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP.xlsx")
    print(list(candidate_drug_df))
    candidate_drug_df = candidate_drug_df[['source','name']]
    f1_matrix_df = pd.merge(left = candidate_drug_df,
                            right = f1_matrix_top3_merged_df,
                            how='left',
                            left_on='source',
                            right_index=True)
    f1_matrix_df = f1_matrix_df.drop(columns=['source'])
    f1_matrix_df.to_csv(f1_matrix_merged_addr,index=False)

    f1_pathway = list(f1_matrix_df)
    print(f1_pathway)
    f1_pathway_df = pd.DataFrame(f1_pathway[1:],columns=['pathway'])
    f1_pathway_df.to_csv(f"../result/{virus}/Drug/F1/{virus}_f1_pahtways.csv", index=False)

def main():
    # make_drug_patway()

    virus = "SARS-CoV"
    # virus = "SARS-CoV-2"
    candidate_drug_df= pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.xlsx")
    candidate_drugs = candidate_drug_df['source'].tolist()
    print(candidate_drugs)

    threshold = 1.0e-13     # SARS-CoV
    # threshold = 1.0e-12     # SARS-CoV-2
    whole_pa_list = get_whole_enriched_pathways(virus, candidate_drug_df, threshold)
    # find_ancestor_enriched_pathways(virus,whole_pa_list, threshold)    # one time
    #
    # #####################
    # ## Step 2
    # #####################
    #
    # drug_to_pathway_df= pd.read_excel( f"../result/{virus}/Drug/{virus}_less5_Reactome_pathways_hierachi_{threshold}_top2.xlsx")
    drug_to_pathway_df= pd.read_excel( f"../result/{virus}/Drug/{virus}_extendAll_Reactome_pathways_hierachi_{threshold}_top3.xlsx")
    low_level_pathway_addr = f"../result/{virus}/Drug/{virus}_extendAll_Reactome_low_level_{threshold}_top3.xlsx"
    # low_level_pathways = find_low_level_pathways(virus,drug_to_pathway_df, threshold,low_level_pathway_addr)

    f1_matrix_addr = f"../result/{virus}/Drug/F1/{virus}_extendAll_drug_to_low_level_pathways_{threshold}_matrix.csv"
    # make_matrix(virus,candidate_drug_df,low_level_pathway_addr,threshold,f1_matrix_addr)
    f1_matrix_merged_addr=f"../result/{virus}/Drug/F1/{virus}_extendAll_drug_to_low_level_pathways_{threshold}_merged_top3_matrix.csv"
    make_matrix_merge_to_top3(virus,low_level_pathway_addr,f1_matrix_addr,threshold,f1_matrix_merged_addr)


if __name__ == '__main__':
    main()
