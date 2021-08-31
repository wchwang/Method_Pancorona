# Created by woochanghwang at 10/06/2020
'''
circos make subcellular
for generated file from "circos_make_input_files_key_genes"
MoA based SOM result
#Modified 07/21/2020
Hidden , Sorted by enriched pathways

'''
import pandas as pd
import numpy as np


def main_HPA():
    key_genes_in_circos_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_mutlcategory_v4.tsv"
    hpa_potential_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/HPA/protein_class_Potential.tsv"
    hpa_fda_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/HPA/protein_class_FDA.tsv"

    hpa_potential_df = pd.read_csv(hpa_potential_addr, sep='\t')
    # hpa_potential_df = hpa_potential_df['Gene', 'Subcellular main location',"Subcellular additional location" ]

    hpa_fda_df = pd.read_csv(hpa_fda_addr, sep='\t')
    hpa_druggable_df = pd.concat([hpa_fda_df,hpa_potential_df])
    print(hpa_potential_df)

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

def main():
    circos_node_file_a = '/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_high_level_paths_v2.tsv'
    subcellular_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_moa_based_som_reverse_subcellular_conf4_v3.tsv"
    # subcellular_df = pd.read_csv(subcellular_addr,sep='\t')
    # print(subcellular_df)

    # ###########################
    # #Step0: make subcellular data for circos
    # #07/14 : I used knowledgebased only form Compartment, confidece level >= 4
    # #############################
    # # make_circos_subcellular(subcellular_addr)
    make_circos_subcellular_knowledgebased(circos_node_file_a,subcellular_addr)

    ###############################
    # Step 1: make circos file added subcellular
    #################################
    circos_data_file_a = '/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_DIP_structure_DEP_HIDDEN_key_genes_high_level_paths_v2.txt'
    genecards_subcellular_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Compartments/genecard_subcellular.txt"
    with open(genecards_subcellular_addr) as subcellular_f:
        subcellular = [x.strip() for x in subcellular_f.readlines()]

    covid_circos_position_dict = get_covid_circos_position(circos_data_file_a)

    for location in subcellular:
        location_no_space= location.replace(' ','_')
        circos_for_a_location_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_high_level_paths_subcellular_conf4_{}_v2.tsv".format(location)
        covid_circos_for_a_location_hist_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_key_genes_high_level_paths_subcellular_hist_{}_v2.txt".format(location_no_space)
        make_circos_for_a_subcellular(subcellular_addr, location, location_no_space, circos_for_a_location_addr, covid_circos_position_dict, covid_circos_for_a_location_hist_addr)


if __name__ == '__main__':
    main()
    # main_HPA()