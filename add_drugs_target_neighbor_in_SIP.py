# Created by woochanghwang at 15/07/2021
# Created by woochanghwang at 18/12/2020
import toolbox.data_handler as dh
import pandas as pd
import numpy as np


def get_neighbor(key_proteins, SIP_G):
    key_proteins = str(key_proteins)
    if key_proteins == 'nan':
        return "NA"
    else:
        key_proteins = key_proteins.split(',')
        key_proteins_neighbor = []
        for key in key_proteins:
            key_proteins_neighbor.append(key)
            key_proteins_neighbor += [n for n in SIP_G.neighbors(key)]

        return ','.join(key_proteins_neighbor)

def get_neighbor_less_5(key_proteins, SIP_G):
    key_proteins = str(key_proteins)
    if key_proteins == 'nan':
        return "NA"
    else:
        key_proteins = key_proteins.split(',')
        if len(key_proteins)>5:
            return ','.join(key_proteins)
        else:
            key_proteins_neighbor = []
            for key in key_proteins:
                key_proteins_neighbor.append(key)
                key_proteins_neighbor += [n for n in SIP_G.neighbors(key)]

            return ','.join(key_proteins_neighbor)

def target_in_network(drug_targets,virus):
    # virus = "SARS-CoV-2"
    print(drug_targets)
    sip_nodes = pd.read_csv(f"../result/{virus}/network/{virus}_network_nodes.csv")['ID'].to_list()

    drug_targets = drug_targets.split(', ')
    print(drug_targets)
    print(set(drug_targets)&set(sip_nodes))

    return ','.join(list(set(drug_targets)&set(sip_nodes)))




def make_extended_target_list(virus):

    candidate_drug_target_df = pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP.xlsx")
    SIP_Graph = dh.load_obj(f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph")
    mode = "all"    # less_5
    if mode == "all":
        candidate_drug_target_df['Target_extended'] = candidate_drug_target_df.apply(lambda row:
                                                                                     get_neighbor(row['Drug_Target_SIP'],SIP_Graph),
                                                                                     axis=1)
    elif mode == 'less_5':
        candidate_drug_target_df['Target_extended'] = candidate_drug_target_df.apply(lambda row:
                                                                                     get_neighbor_less_5(row['Drug_Target_SIP'],SIP_Graph),
                                                                                     axis=1)
    candidate_drug_target_df.to_csv(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP_extended_{mode}.tsv", sep='\t', index=False)

def add_target_info_on_candidate_drugs(virus):


    drug_np_target_df = pd.read_csv(f"../result/{virus}/Drug/{virus}_drug_network_proximity_target.tsv", sep='\t',index_col=0)
    drug_np_info_df = pd.read_csv(f"../result/{virus}/Drug/{virus}_drug_network_proximity_druginfo.csv")
    drug_np_info_df = drug_np_info_df.set_index('source')

    drug_np_target_df = drug_np_target_df[['Drug_Target']]
    drug_np_info_df = drug_np_info_df[['Disease', 'n.source', 'n.target', 'd', 'z', 'name', 'type', 'groups', 'atc_codes', 'categories', 'indication']]
    drug_np_info_target_df = pd.concat([drug_np_info_df,drug_np_target_df], axis=1)
    drug_np_info_target_df = drug_np_info_target_df.reset_index()
    print(drug_np_info_target_df)

    candidate_drug_target_df = drug_np_info_target_df[drug_np_info_target_df['z']<=-2]
    candidate_drug_target_df.to_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.xlsx")
    print(candidate_drug_target_df)

    # candidate_drug_target_df['Drug_Target_SIP'] = candidate_drug_target_df['Drug_Target'].apply(target_in_network)
    candidate_drug_target_df['Drug_Target_SIP'] = candidate_drug_target_df.apply(lambda row: target_in_network(row['Drug_Target'],virus), axis=1)
    # # candidate_drug_target_df['Drug_Target_SIP'] = candidate_drug_target_df['Drug_Targetlist'].apply(target_in_network)
    candidate_drug_target_df['Drug_Target_SIPlist'] = candidate_drug_target_df['Drug_Target_SIP'].str.split(',')
    candidate_drug_target_df['n.Drug_Target_SIP'] = candidate_drug_target_df['Drug_Target_SIPlist'].str.len()

    candidate_drug_target_df = candidate_drug_target_df.drop(columns=['Drug_Target_SIPlist'])

    candidate_drug_target_df.to_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP.xlsx", index=False)

def main():
    '''
    find enriched pathway
    :return:
    '''

    virus = "SARS-CoV"
    # virus = "SARS-CoV-2"
    sip_G = dh.load_obj(f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph")
    sip_G_nodes = sip_G.nodes()
    sip_G_nodes_df = pd.DataFrame(sip_G_nodes,columns=['ID'])
    sip_G_nodes_df.to_csv(f"../result/{virus}/network/{virus}_network_nodes.csv")
    # add_target_info_on_candidate_drugs(virus)
    make_extended_target_list(virus)


if __name__ == '__main__':
    main()
