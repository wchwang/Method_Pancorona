# Created by woochanghwang at 10/06/2020
'''
circos make subcellular
for generated file from "circos_make_input_files_key_genes"
MoA based SOM result
#Modified 07/21/2020
#Modified 07/19/2021
Hidden , Sorted by enriched pathways
Centrality

'''
import pandas as pd
import numpy as np



def make_circos_subcellular(subcellular_addr):
    key_genes_in_circos_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_mutlcategory_v4.tsv"
    genecards_subcellular_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Compartments/genecard_subcellular.txt"
    comparment_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Compartments/human_compartment_integrated_full.tsv"

    with open(genecards_subcellular_addr) as subcellular_f:
        subcellular = [x.strip() for x in subcellular_f.readlines()]

    comparments_data = pd.read_csv(comparment_addr, sep='\t', names=['Ensembl','Symbol','GO','Subcellular','Confidence'])

    print(comparments_data.head())

    comparments_data = comparments_data[comparments_data['Confidence'] >=4]
    print(comparments_data)
    comparments_data = comparments_data[comparments_data['Subcellular'].isin(subcellular)]
    print(comparments_data)
    comparments_data = comparments_data[['Symbol','Subcellular','Confidence']]
    print(comparments_data)

    key_genes_in_circos_df = pd.read_csv(key_genes_in_circos_addr, sep='\t', names=['Band','Protein'])
    print(key_genes_in_circos_df)

    key_genes_in_circos_subcelluar = pd.merge(left=key_genes_in_circos_df, right=comparments_data,how='left',left_on='Protein',right_on='Symbol')
    print(key_genes_in_circos_subcelluar)
    key_genes_in_circos_subcelluar = key_genes_in_circos_subcelluar.fillna('NA')
    key_genes_proteins = key_genes_in_circos_df['Protein'].to_list()
    key_genes_symbols = key_genes_in_circos_subcelluar['Symbol'].to_list()

    print(len(set(key_genes_proteins)))
    print(len(set(key_genes_symbols)))
    print(len(set(key_genes_proteins)&set(key_genes_symbols)))
    print(set(key_genes_proteins)-set(key_genes_symbols))

    key_genes_in_circos_subcelluar = key_genes_in_circos_subcelluar[['Band','Protein','Subcellular','Confidence']]
    print(key_genes_in_circos_subcelluar)
    key_genes_in_circos_subcelluar.to_csv(subcellular_addr,sep='\t',index=False)

def make_circos_subcellular_knowledgebased(key_genes_in_circos_addr, subcellular_addr):
    # key_genes_in_circos_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_mutlcategory_v4.tsv"
    genecards_subcellular_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Compartments/genecard_subcellular.txt"
    comparment_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Compartments/human_compartment_knowledge_full.tsv"

    with open(genecards_subcellular_addr) as subcellular_f:
        subcellular = [x.strip() for x in subcellular_f.readlines()]

    comparments_data = pd.read_csv(comparment_addr, sep='\t', names=['Ensembl','Symbol','GO','Subcellular','Source','Evidence_code','Confidence'])

    print(comparments_data.head())

    # comparments_data = comparments_data[comparments_data['Evidence_code'].isin(['ISS','IDA','HDA', ])] #Inferred from Direct Assay (IDA),

    comparments_data = comparments_data[comparments_data['Confidence'] >=4]
    print(comparments_data)
    comparments_data = comparments_data[comparments_data['Subcellular'].isin(subcellular)]  # multiple source --> mean
    print(comparments_data)

    comparments_data = comparments_data[['Symbol', 'Subcellular', 'Confidence']]
    comparments_data_groupby = comparments_data.groupby(['Symbol','Subcellular']).agg({'Confidence':'mean'})
    comparments_data_groupby = comparments_data_groupby.reset_index()
    print(comparments_data_groupby)


    key_genes_in_circos_df = pd.read_csv(key_genes_in_circos_addr, sep='\t', names=['Band','Protein'])
    print(key_genes_in_circos_df)

    key_genes_in_circos_subcelluar = pd.merge(left=key_genes_in_circos_df, right=comparments_data_groupby,how='left',left_on='Protein',right_on='Symbol')
    print(key_genes_in_circos_subcelluar)
    key_genes_in_circos_subcelluar = key_genes_in_circos_subcelluar.fillna('NA')
    key_genes_proteins = key_genes_in_circos_df['Protein'].to_list()
    key_genes_symbols = key_genes_in_circos_subcelluar['Symbol'].to_list()

    print(len(set(key_genes_proteins)))
    print(len(set(key_genes_symbols)))
    print(len(set(key_genes_proteins)&set(key_genes_symbols)))
    print(set(key_genes_proteins)-set(key_genes_symbols))

    key_genes_in_circos_subcelluar = key_genes_in_circos_subcelluar[['Band','Protein','Subcellular','Confidence']]
    print(key_genes_in_circos_subcelluar)
    key_genes_in_circos_subcelluar.to_csv(subcellular_addr,sep='\t',index=False)

def make_circos_for_a_subcellular(subcellular_addr, location, location_no_space, location_addr,covid_circos_position_dict, covid_circos_for_a_location_hist_addr):

    circos_subcellular_df = pd.read_csv(subcellular_addr, sep='\t')
    circos_for_a_location_df = circos_subcellular_df[circos_subcellular_df['Subcellular']==location]
    circos_for_a_location_df.to_csv(location_addr, sep='\t', index=False)

    genes_in_a_location = circos_for_a_location_df['Protein'].to_list()

    circos_covid_a_location = []

    # for gene, pos in covid_circos_position_dict.items():
    #     print(gene,pos)

    for gene in genes_in_a_location:

        gene_position = covid_circos_position_dict[gene]
        gene_position_sub = gene_position[:]
        gene_position_sub.append('1')
        gene_position_sub.append("id={}".format(location_no_space))
        print(gene_position)
        circos_covid_a_location.append(gene_position_sub)

    circos_covid_a_location = ['\t'.join(gene) for gene in circos_covid_a_location]

    with open(covid_circos_for_a_location_hist_addr,'w') as covid_circos_f:
        covid_circos_f.write('\n'.join(circos_covid_a_location))


def get_covid_circos_position(circos_position_addr):
    # circos_position_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_DIP_structure_DEP_HIDDEN_key_genes_in_network.txt"
    with open(circos_position_addr) as circos_position_f:
        circos_position_list = [x.strip().split('\t') for x in circos_position_f.readlines()]

    circos_position_dict = dict()
    for gene in circos_position_list:
        # circos_position_dict[gene[-1]] = '\t'.join(gene[:-1])
        circos_position_dict[gene[-1]] = gene[:-1]
    return circos_position_dict

def make_circos_for_a_centrality(circos_centrality_base_addr,centrality,circos_for_a_centrality_addr, covid_circos_position_dict, covid_circos_for_a_location_hist_addr):

    circos_centrality_df = pd.read_csv(circos_centrality_base_addr, sep='\t')
    circos_for_a_location_df = circos_centrality_df[circos_centrality_df['Centrality']==centrality]
    circos_for_a_location_df.to_csv(circos_for_a_centrality_addr, sep='\t', index=False)

    genes_in_a_location = circos_for_a_location_df['Protein'].to_list()

    circos_covid_a_location = []

    # for gene, pos in covid_circos_position_dict.items():
    #     print(gene,pos)

    for gene in genes_in_a_location:

        gene_position = covid_circos_position_dict[gene]
        gene_position_sub = gene_position[:]
        gene_position_sub.append('1')
        gene_position_sub.append("id={}".format(centrality))
        print(gene_position)
        circos_covid_a_location.append(gene_position_sub)

    circos_covid_a_location = ['\t'.join(gene) for gene in circos_covid_a_location]

    with open(covid_circos_for_a_location_hist_addr,'w') as covid_circos_f:
        covid_circos_f.write('\n'.join(circos_covid_a_location))

def make_circos_centrality_base(virus, circos_node_file_a, keyprotein_addr, circos_centrality_addr):
    key_gene_SARS_df = pd.read_csv(keyprotein_addr)
    key_gene_SARS = key_gene_SARS_df['Gene'].to_list()

    network_anaysis_df = pd.read_csv(f"../result/{virus}/network_analysis/{virus}_A549_24h_centrality_RWR_result_pvalue.csv")

    eigen_list = network_anaysis_df[network_anaysis_df['Eigen_pvalue']< 0.01]['Gene'].tolist()
    degree_list = network_anaysis_df[network_anaysis_df['Degree_pvalue']< 0.01]['Gene'].tolist()
    bw_list = network_anaysis_df[network_anaysis_df['Between_plvaue']< 0.01]['Gene'].tolist()
    rwr_list = network_anaysis_df[network_anaysis_df['RWR_pvalue']< 0.01]['Gene'].tolist()

    key_genes = list(set(key_gene_SARS))
    key_genes_eigen = list(set(key_gene_SARS) & set(eigen_list))
    key_genes_degree = list(set(key_gene_SARS) & set(degree_list) - set(key_genes_eigen))
    key_genes_bw =  list(set(key_gene_SARS) & set(bw_list) - set(key_genes_eigen)-set(key_genes_degree))
    key_genes_rwr =  list(set(key_gene_SARS) & set(rwr_list) - set(key_genes_eigen)-set(key_genes_degree)-set(key_genes_bw))

    key_genes_centrality = []
    for gene in key_genes_eigen:
        key_genes_centrality.append([gene, 'Eigen'])

    for gene in key_genes_degree:
        key_genes_centrality.append([gene, 'Degree'])

    for gene in key_genes_bw:
        key_genes_centrality.append([gene, 'BW'])

    for gene in key_genes_rwr:
        key_genes_centrality.append([gene, 'RWR'])

    key_genes_centrality_df = pd.DataFrame(key_genes_centrality,columns=['Protein','Centrality'])
    key_genes_in_circos_df = pd.read_csv(circos_node_file_a, sep='\t', names=['Band','Protein'])
    print(key_genes_in_circos_df)

    key_genes_in_circos_centrality = pd.merge(left=key_genes_in_circos_df, right=key_genes_centrality_df,how='left',left_on='Protein',right_on='Protein')
    key_genes_in_circos_centrality = key_genes_in_circos_centrality.dropna()

    key_genes_in_circos_centrality = key_genes_in_circos_centrality[['Band','Protein','Centrality']]
    print(key_genes_in_circos_centrality)
    key_genes_in_circos_centrality.to_csv(circos_centrality_addr,sep='\t',index=False)



def main():
    viruslist = ["SARS-CoV","SARS-CoV-2"]
    for virus in viruslist:
        circos_node_file_a = f"../result/{virus}/Circos/data/{virus}_backbone_node_info_key_genes_high_level_paths.tsv"
        keyprotein_addr = f"../result/{virus}/Sig.Genes/{virus}_key_protein_every.txt"
        circos_centrality_base_addr = f"../result/{virus}/Circos/data/{virus}_backbone_node_info_key_genes_high_level_paths_centrality.tsv"

        make_circos_centrality_base(virus, circos_node_file_a,keyprotein_addr,circos_centrality_base_addr)

        circos_data_file_a = f"../result/{virus}/Circos/data/{virus}_DIP_structure_DEP_HIDDEN_key_genes_high_level_paths.txt"
        covid_circos_position_dict = get_covid_circos_position(circos_data_file_a)


        centrality_list = ["Eigen", "Degree", "BW", "RWR"]
        for centrality in centrality_list:
            circos_for_a_centrality_addr = f"../result/{virus}/Circos/data/{virus}_backbone_node_info_key_genes_high_level_paths_centrality_{centrality}.tsv"
            covid_circos_for_a_location_hist_addr = f"../result/{virus}/Circos/data/{virus}_key_genes_high_level_paths_centrality_hist_{centrality}.txt"
            make_circos_for_a_centrality(circos_centrality_base_addr,centrality,circos_for_a_centrality_addr, covid_circos_position_dict, covid_circos_for_a_location_hist_addr)

        # ###############################
        # # Step 1: make circos file added subcellular
        # #################################
        # circos_data_file_a = '/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_DIP_structure_DEP_HIDDEN_key_genes_high_level_paths_v2.txt'
        # genecards_subcellular_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Compartments/genecard_subcellular.txt"
        # with open(genecards_subcellular_addr) as subcellular_f:
        #     subcellular = [x.strip() for x in subcellular_f.readlines()]
        #
        # covid_circos_position_dict = get_covid_circos_position(circos_data_file_a)
        #
        # for location in subcellular:
        #     location_no_space= location.replace(' ','_')
        #     circos_for_a_location_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_high_level_paths_subcellular_conf4_{}_v2.tsv".format(location)
        #     covid_circos_for_a_location_hist_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_key_genes_high_level_paths_subcellular_hist_{}_v2.txt".format(location_no_space)
        #     make_circos_for_a_subcellular(subcellular_addr, location, location_no_space, circos_for_a_location_addr, covid_circos_position_dict, covid_circos_for_a_location_hist_addr)


if __name__ == '__main__':
    main()