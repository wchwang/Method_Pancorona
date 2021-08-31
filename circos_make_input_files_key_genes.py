# Created by woochanghwang at 21/05/2020
# Modified by woochanghwang at 12/07/2020
# Modified at 18/07/2020
# Modified at 19/07/2021 for Method paper
'''
Make input files for CIRCUS
including selected drug targets in the network, key genes
# Modified by woochanghwang at 21/07/2020
- Hidden enrichment test (High level paths)
- Metabolism of RNA, Immune System, Cell Cycle, Other
# Modified by Woochang at 19/07/2021
- for Method paper
'''

import pandas as pd
import csv
import networkx as nx
import toolbox.data_handler as dh
import toolbox.visual_utilities as vu


def make_link_file(network_edge, circos_data_file_a, circos_link_file_a):
    '''
    - edge file : [from to]
    - circos data file : [chr   band    start(number)   to(number)  name]
    :return: circos link file [ band    start   to  band    start   to  ]
    '''

    # network_edge_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_visual/ULK1_sigGene_diNetwork.tsv"
    # circos_data_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/circos/ulk1/data/ulk1_sigGene_cluster_text.txt"
    # circos_link_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/circos/ulk1/data/ulk1_sigGene_cluster_link_v2.txt"

    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    circos_link = []
    for edge in network_edge:
        link_from = circos_band_dict.get(edge[0],'NA')
        link_to = circos_band_dict.get(edge[1],'NA')

        if link_from  == 'NA' or link_to == 'NA' :  continue

        circos_link.append('\t'.join(link_from+link_to))


    with open(circos_link_file_a,'w') as circos_link_f:
        circos_link_f.write('\n'.join(circos_link))



def make_link_with_symbol_file(network_edge_addr, circos_data_file_a, circos_link_file_a):
    '''
    - edge file : [from to]
    - circos data file : [chr   band    start(number)   to(number)  name]
    :return: circos link file [ band    start   to  band    start   to  ]
    '''

    # network_edge_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_visual/ULK1_sigGene_diNetwork.tsv"
    # circos_data_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/circos/ulk1/data/ulk1_sigGene_cluster_text.txt"
    # circos_link_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/circos/ulk1/data/ulk1_sigGene_cluster_link_v2.txt"

    with open(network_edge_addr) as network_edge_f:
        network_edge = [x.strip().split('\t') for x in network_edge_f.readlines()]


    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        # gene_pos = line[:-1]
        gene_band = line[0]
        circos_band_dict[gene] = gene_band

    circos_link = []

    # print(network_edge[:5])
    for edge in network_edge:
        link_from = circos_band_dict.get(edge[0],'NA')
        link_to = circos_band_dict.get(edge[1],'NA')

        if link_from  == 'NA' or link_to == 'NA' :  continue

        circos_link.append([edge[0],link_from,edge[1],link_to])
    print(circos_link[:5])

    circos_link_df = pd.DataFrame(circos_link, columns=["Gene A","Gene A Group","Gene B","Gene B Group"])

    circos_link_df.to_excel(circos_link_file_a,index=False)

    # with open(circos_link_file_a,'w') as circos_link_f:
    #     circos_link_f.write('\n'.join(circos_link))

def make_link_file_for_specific_groups(network_edge, groups, circos_data_file_a, circos_link_file_a):
    '''
    - edge file : [from to]
    - circos data file : [chr   band    start(number)   to(number)  name]
    :return: circos link file [ band    start   to  band    start   to  ]
    '''

    # network_edge_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_visual/ULK1_sigGene_diNetwork.tsv"
    # circos_data_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/circos/ulk1/data/ulk1_sigGene_cluster_text.txt"
    # circos_link_file_a = "/Users/woochanghwang/PycharmProjects/LifeArc/General/circos/ulk1/data/ulk1_sigGene_cluster_link_v2.txt"

    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    group_genes = []

    for circos_genes in circos_band_data:
        if circos_genes[0] in groups:
            # print(circos_genes)
            group_genes.append(circos_genes[-1])

    # print(group_genes)

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    circos_link = []
    for edge in network_edge:
        if len(set(edge)&set(group_genes)) >=1:
            link_from = circos_band_dict.get(edge[0],'NA')
            link_to = circos_band_dict.get(edge[1],'NA')

            if link_from  == 'NA' or link_to == 'NA' :  continue

            circos_link.append('\t'.join(link_from+link_to))


    with open(circos_link_file_a,'w') as circos_link_f:
        circos_link_f.write('\n'.join(circos_link))

def make_rwr_hist_file_for_heatmap (circos_data_file_a, rwr_result_file_addr, rwr_hist_file_a):
    '''
        Gene RWR for histogram
        :return: [clsuter, gene, fold change]
        '''

    with open(rwr_result_file_addr) as rwr_result_f:
        gene_rwr_result = [x.strip().split('\t') for x in rwr_result_f.readlines()]

    gene_rwr_list = [x[:2] for x in gene_rwr_result[1:]]    # header = False


    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    circos_gene_fc_data = []

    for gene_rwr in gene_rwr_list:
        gene_pos = circos_band_dict[gene_rwr[0]]
        rwr = str(gene_rwr[1])
        a_gene = [gene_pos[0], gene_pos[1], gene_pos[1], rwr, 'id=rwr']
        circos_gene_fc_data.append('\t'.join(a_gene))


    with open(rwr_hist_file_a, 'w') as circos_gene_rwr_f:
        circos_gene_rwr_f.write('\n'.join(circos_gene_fc_data))



def get_whole_node():
    core_network_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/data/string_interactions_GBM_ULK_e150_in_TCGA.tsv"
    # core_network_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/data/string_interactions_new_symbol.tsv"

    core_network_df= pd.read_table(core_network_addr,sep='\t')

    core_G = nx.from_pandas_edgelist(core_network_df,'node1','node2')

    print(len(core_G.nodes))

    key_gene_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/data/GBM_OT_TCGA_ULK.txt"

    key_gene_df = pd.read_table(key_gene_addr,sep='\t')

    key_gene = key_gene_df['Gene']

    print(len(set(key_gene)-set(core_G.nodes())))
    from_ppi_genes = list(set(core_G.nodes()) - set(key_gene))

    new_gene_list = []

    for gene in from_ppi_genes:
        new_gene_list.append(['STRING',gene])

    for string_gene in new_gene_list:
        key_gene_df = key_gene_df.append(pd.Series(string_gene,index=['Group','Gene']),ignore_index=True)

    print(key_gene_df)

    key_gene_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/ULK/data/GBM_OT_TCGA_STRING_ULK.txt",'\t',index=False,header=False)



def make_tcga_fc_file_for_heatmap(tcga_fc_addr, circos_data_file_a, circos_gene_fc_file_a):

    tcga_logfc_result = dh.load_obj(tcga_fc_addr)

    circos_data_df = pd.read_csv(circos_data_file_a,sep='\t',names=['Group','posA','posB','Gene'])
    whole_gene = list(circos_data_df['Gene'])
    print(whole_gene)

    # df_final_result = pd.DataFrame(columns=['Gene', 'FoldChange(log2)'])

    whole_fc_gene = []
    for gene in whole_gene:

        whole_fc_gene.append(tcga_logfc_result.get(gene, ''))

    gene_fc_list = list(zip(whole_gene,whole_fc_gene))
    # gene_fc_df = pd.DataFrame(gene_fc_list, columns=['Gene', 'FoldChange(log2)'])

    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    # for gene, gene_pos in circos_band_dict.items():
    #     gene_fc_df = gene_fc_df.replace(gene, '\t'.join(gene_pos))

    circos_gene_fc_data = []

    for gene_fc in gene_fc_list:
        gene_pos = circos_band_dict[gene_fc[0]]
        fc = str(gene_fc[1])
        if fc == 'inf': fc = ''
        a_gene = [gene_pos[0],gene_pos[1], gene_pos[1],fc,'id=fc']
        circos_gene_fc_data.append('\t'.join(a_gene))



    with open(circos_gene_fc_file_a,'w') as circos_gene_fc_f:
        circos_gene_fc_f.write('\n'.join(circos_gene_fc_data))



def make_drug_score_file_for_heatmap(drug_score_file, circos_data_file_a, circos_data_final_file_a):
    '''
        Gene RWR for histogram
        :return: [clsuter, gene, fold change]
        '''


    with open(drug_score_file) as drug_score_f:
        drug_score_list = [x.strip().split('\t') for x in drug_score_f.readlines()]

    drug_score_list = drug_score_list[1:]   #header = False

    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    circos_gene_drug_data = []
    for drug_score in drug_score_list:
        gene_pos = circos_band_dict.get(drug_score[0],'NA')
        if gene_pos == 'NA':    continue
        score = str(drug_score[1])
        a_gene = [gene_pos[0], gene_pos[1], gene_pos[1], score, 'id=drug']
        circos_gene_drug_data.append('\t'.join(a_gene))


    with open(circos_data_final_file_a, 'w') as circos_gene_drug_f:
        circos_gene_drug_f.write('\n'.join(circos_gene_drug_data))

def make_depmap_score_file_for_heatmap(depmap_score_file_a, circos_data_file_a, circos_depmap_score_file_a):

    with open(depmap_score_file_a) as depmap_score_f:
        drug_depmap_data = [x.strip().split('\t') for x in depmap_score_f.readlines()]

    depmap_score_list = []
    for drug in drug_depmap_data[1:]:   #Header = False
        a_depmap = [drug[0],drug[-1]]
        depmap_score_list.append(a_depmap)

    print(depmap_score_list[:10])
    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    circos_gene_depmap_data = []
    for depmap_score in depmap_score_list:
        gene_pos = circos_band_dict.get(depmap_score[0],'NA')
        if gene_pos == 'NA':    continue
        score = str(depmap_score[1])
        a_gene = [gene_pos[0], gene_pos[1], gene_pos[1], score, 'id=depmap']
        circos_gene_depmap_data.append('\t'.join(a_gene))


    with open(circos_depmap_score_file_a, 'w') as circos_gene_depmap_f:
        circos_gene_depmap_f.write('\n'.join(circos_gene_depmap_data))


def make_final_target_file_for_heatmap(gene_list_file_ulk1,gene_list_file_ulk2, circos_data_file_a, circos_final_targets_file_a):
    '''
        Gene RWR result from files
        :return: [clsuter, gene, fold change]
        '''


    result_ulk1_df = pd.read_table(gene_list_file_ulk1, sep='\t')
    result_ulk2_df = pd.read_table(gene_list_file_ulk2, sep='\t')

    print(list(result_ulk1_df))
    print(list(result_ulk2_df))

    final_col_names=['Gene','Final score normalized']
    ulk1_final_df = result_ulk1_df[final_col_names]
    ulk2_final_df = result_ulk2_df[final_col_names]

    ulk1_top20_df = ulk1_final_df.sort_values(by=['Final score normalized'],ascending=False).iloc[:20]
    ulk2_top20_df = ulk2_final_df.sort_values(by=['Final score normalized'],ascending=False).iloc[:20]


    ulk1_2_final_top20_df = pd.concat([ulk1_top20_df,ulk2_top20_df]).drop_duplicates(subset='Gene',keep='last').reset_index(drop=True)

    # #################
    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        gene_pos = line[:-1]
        circos_band_dict[gene] = gene_pos

    for gene,gene_pos in circos_band_dict.items():
        ulk1_2_final_top20_df = ulk1_2_final_top20_df.replace(gene,'\t'.join(gene_pos))

    ulk1_2_final_top20_df.insert(2,'id','id=fc')

    print(ulk1_2_final_top20_df)

    ulk1_2_final_top20_df[['cluster','pos_x','pos_y']]=ulk1_2_final_top20_df.Gene.str.split("\t",expand=True)

    ulk1_2_final_top20_df_write = ulk1_2_final_top20_df[['cluster','pos_x','pos_y','Final score normalized','id']]

    print(ulk1_2_final_top20_df_write)
    ulk1_2_final_top20_df_write.to_csv(circos_final_targets_file_a, sep='\t', header=False,
                               index=False,quoting=csv.QUOTE_NONE,quotechar="",  escapechar='\\')


def make_final_target_file_for_text(gene_list_file_ulk1,gene_list_file_ulk2, circos_data_file_a, circos_final_targets_file_a):

    result_ulk1_df = pd.read_table(gene_list_file_ulk1, sep='\t')
    result_ulk2_df = pd.read_table(gene_list_file_ulk2, sep='\t')

    print(list(result_ulk1_df))
    print(list(result_ulk2_df))

    final_col_names=['Gene','Final score normalized']
    ulk1_final_df = result_ulk1_df[final_col_names]
    ulk2_final_df = result_ulk2_df[final_col_names]

    ulk1_top20_df = ulk1_final_df.sort_values(by=['Final score normalized'],ascending=False).iloc[:20]
    ulk2_top20_df = ulk2_final_df.sort_values(by=['Final score normalized'],ascending=False).iloc[:20]


    ulk1_2_final_top20_df = pd.concat([ulk1_top20_df,ulk2_top20_df]).drop_duplicates(subset='Gene',keep='last').reset_index(drop=True)

    # #################
    circos_band_dict = dict()

    with open(circos_data_file_a) as circos_band_f:
        circos_band_data = [x.strip().split('\t') for x in circos_band_f.readlines()]

    for line in circos_band_data:
        gene = line[-1]
        # gene_pos = line[:-1]
        gene_pos = line # too keep gene name
        circos_band_dict[gene] = gene_pos

    for gene,gene_pos in circos_band_dict.items():

        # print(gene,gene_pos)
        ulk1_2_final_top20_df = ulk1_2_final_top20_df.replace(gene,'\t'.join(gene_pos))


    print(ulk1_2_final_top20_df)

    ulk1_2_final_top20_df[['cluster','pos_x','pos_y','Gene']]=ulk1_2_final_top20_df.Gene.str.split("\t",expand=True)

    ulk1_2_final_top20_df_text= ulk1_2_final_top20_df[['cluster','pos_x','pos_y','Gene']]

    print(ulk1_2_final_top20_df_text)
    ulk1_2_final_top20_df_text.to_csv(circos_final_targets_file_a, sep='\t', header=False,
                               index=False,quoting=csv.QUOTE_NONE,quotechar="",  escapechar='\\')




def make_color_file(virus, circos_data_file_a, circos_color_file_a):

    circos_data_df = pd.read_csv(circos_data_file_a, sep='\t', names=['Group','posA','posB','Gene'])

    groups = circos_data_df.Group.unique()
    print(groups)
    group_size_df = circos_data_df.groupby('Group').size()
    print(group_size_df['hs1'])
    circos_color_info = []
    ##################
    # Add Chr info
    ##################
    # grp_names = ['DIP','Virus_entry','Virus_replication','Virus/Immune','Metabolism','Anti-inflammatory','Immune_system','Ohter','DEP']
    # grp_names = ['E','M','N','nsp1','nsp2','nsp4','nsp5','nsp6','nsp7','nsp8','nsp9','nsp10','nsp12','nsp13','nsp14','nsp15','orf3a','orf6','orf7a','orf8','orf9b','orf9c','orf10',
    #              'VR_VE_MB','VE_MB','DEP','VR','AI_IR','Unknwown']
    ###v2
    # grp_names = ['E', 'M', 'N', '1', '2', '4', '5', '6', '7', '8', '9', '10', '12',
    #              '13', '14', '15', '3a', '6', '7a', '8', '9b', '9c', '10',
    #              'VRVEMB', 'DEP', 'VR', 'AIIR', 'Unknown']
    ###Localization
    if virus == "SARS-CoV":
        grp_names = ['E', 'M', 'N', '2', '6', '9', '10', '12',
                     '13', '15', '16', '3a','3b', '7a','7b', '8a','8b', '9b','S',
                     'Met.RNA','CellCycle','Immune','Met.Proteins', 'Other','DEP']
    elif virus =="SARS-CoV-2":
        grp_names = ['E', 'M', 'N', '1','4', '6', '9', '3', '7a','7b','8',
                     'Met.RNA','CellCycle','Immune','Met.Proteins', 'Other','DEP']
    # grp_names = ['T','PPT','AG','DRUG','TF','TCGA','FRANK','MIDDLE']
    chr_number = 0
    for grp,name in zip(groups,grp_names):
        chr_number += 1
        # print(grp, name)
        a_grp = ['chr','-',grp,name,'0',str(group_size_df[grp]),'chr'+str(chr_number)]
        # print(a_grp)
        circos_color_info.append('\t'.join(a_grp))

    ################
    # Add Band Info
    ################
    with open(circos_data_file_a) as circos_data_f:
        circos_dota_for_color = [x.strip().split('\t') for x in circos_data_f.readlines()]

    #band	hs1	ULK1	ULK1	0	1	grey
    for gene in circos_dota_for_color:
        a_band = ['band',gene[0],gene[-1],gene[-1],gene[1],gene[2],'grey']
        circos_color_info.append('\t'.join(a_band))

    ##############
    # make file
    ##############

    with open(circos_color_file_a,'w') as circos_color_f:
        circos_color_f.write('\n'.join(circos_color_info))


def make_circos_data_file_for_covid(virus, circos_data_file_a, circos_nodes_file_a):

    circos_covid_network_df = pd.read_csv(circos_nodes_file_a, sep='\t', names=['Group','Gene'])

    circos_covid_network_df = circos_covid_network_df.drop_duplicates(subset='Gene')

    circos_covid_network_df = circos_covid_network_df.reset_index(drop=True)


    circos_group_list = list(circos_covid_network_df['Group'])

    circos_posA_list = []
    circos_posB_list = []

    cur_group = circos_group_list[0]
    posA=0
    posB=1


    for group in circos_group_list:
        if group == cur_group:
            circos_posA_list.append(posA)
            circos_posB_list.append(posB)
            posA +=1
            posB += 1

        else:
            cur_group = group
            posA=0
            posB=1
            circos_posA_list.append(posA)
            circos_posB_list.append(posB)
            posA +=1
            posB +=1

    circos_covid_network_df['posA'] = pd.Series(circos_posA_list)
    circos_covid_network_df['posB'] = pd.Series(circos_posB_list)

    circos_covid_network_df = circos_covid_network_df[['Group','posA','posB','Gene']]
    # circos_covid_network_df = circos_covid_network_df.astype({'posA':int, 'posB':int})


    if virus == "SARS-CoV":
        ###localization
        circos_covid_network_df = circos_covid_network_df.replace({'Group':{'E':'hs1',
                                                                            'M':'hs2',
                                                                            'N':'hs3',
                                                                            'S':'hs4',
                                                                            'NSP2':'hs5',
                                                                            'NSP6':'hs6',
                                                                            'NSP9':'hs7',
                                                                            'NSP10':'hs8',
                                                                            'NSP12': 'hs9',
                                                                            'NSP13': 'hs10',
                                                                            'NSP15': 'hs11',
                                                                            'NSP16': 'hs12',
                                                                            'ORF3a':'hs13',
                                                                            'ORF3b':'hs14',
                                                                            'ORF7a':'hs15',
                                                                            'ORF7b':'hs16',
                                                                            'ORF8a':'hs17',
                                                                            'ORF8b':'hs18',
                                                                            'ORF9b':'hs19',
                                                                            'Met.RNA': 'hs20',
                                                                            'CellCycle': 'hs21',
                                                                            'Immune': 'hs22',
                                                                            'Met.Proteins': 'hs23',
                                                                            'Other': 'hs24',
                                                                            'DEP': 'hs25'
                                                                            }})

    elif virus == "SARS-CoV-2":
        circos_covid_network_df = circos_covid_network_df.replace({'Group':{'E':'hs1',
                                                                            'M':'hs2',
                                                                            'N':'hs3',
                                                                            'NSP1':'hs4',
                                                                            'NSP4':'hs5',
                                                                            'NSP6':'hs6',
                                                                            'NSP9':'hs7',
                                                                            'ORF3':'hs8',
                                                                            'ORF7a':'hs9',
                                                                            'ORF7b':'hs10',
                                                                            'ORF8':'hs11',
                                                                            'Met.RNA': 'hs12',
                                                                            'CellCycle': 'hs13',
                                                                            'Immune': 'hs14',
                                                                            'Met.Proteins': 'hs15',
                                                                            'Other': 'hs16',
                                                                            'DEP': 'hs17'
                                                                            }})
    print(circos_covid_network_df)

    circos_covid_network_df.to_csv(circos_data_file_a,sep='\t',index=False,header=False)


def get_drug_targets_in_network(network_time,candidate_drug_df):
    # print(network_time)
    candidate_drugs_in_networktime = candidate_drug_df[candidate_drug_df['Network_time']==network_time]
    candidate_drugs_target_list = candidate_drugs_in_networktime['Target Proteins in Network'].to_list()

    drug_targets  = []
    for targets in candidate_drugs_target_list:
        drug_targets += targets.split(',')
    drug_targets = list(set(drug_targets))

    return drug_targets
import re
def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)',text) ]


def get_sig_genes_by_centrality(centrality_results_6hr, centrality_bands):

    threshold = 0.01
    sig_genes_dict = dict()
    for centrality in centrality_bands:
        centrality_pvalue_col = "{}_pvalue".format(centrality)
        sig_genes = centrality_results_6hr[centrality_results_6hr[centrality_pvalue_col]<=threshold]['Gene'].to_list()
        sig_genes_dict[centrality] = sig_genes

    # common_genes = []
    # for key, value in sig_genes_dict.items():
    #     print(key, len(value),value[:5])
    #     if len(common_genes) == 0:
    #         common_genes = value
    #     else:
    #         common_genes = list(set(common_genes)&set(value))
    #
    # print("common:",len(common_genes))
    # sig_genes_without_common = dict()



    return sig_genes_dict


def sort_hidden_by_centrality(covid_network_hidden_keyGene, centrality_bands):
    centrality_results_6hr = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Centrality/6hr/round2/COVID_6hr_gene_score_by_centrality_pvalue.csv")
    centrality_results_24hr = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Centrality/24hr/round6/COVID_24hr_gene_score_by_centrality_pvalue.csv")

    centrality_genes_6hr = get_sig_genes_by_centrality(centrality_results_6hr,centrality_bands)
    centrality_genes_24hr = get_sig_genes_by_centrality(centrality_results_24hr, centrality_bands)

    ######################
    # for venn
    ######################
    centrality_genes_all_for_venn = {}
    centrality_genes_24hr_hidden_for_venn = {}
    centrality_genes_6hr_hidden_for_venn ={}

    ###################

    centrality_genes_in_hidden_dict = dict()
    common_genes = []
    for centrality in centrality_bands:
        sig_genes_6hr = centrality_genes_6hr[centrality]
        sig_genes_24hr = centrality_genes_24hr[centrality]
        centality_sig_genes= list(set(sig_genes_6hr).union(set(sig_genes_24hr)))
        centality_sig_genes_in_hidden = list(set(centality_sig_genes)&set(covid_network_hidden_keyGene))
        if len(common_genes) == 0:
            common_genes = centality_sig_genes_in_hidden
        else:
            common_genes = list(set(common_genes)&set(centality_sig_genes_in_hidden))

        centrality_genes_in_hidden_dict[centrality] = sorted(centality_sig_genes_in_hidden)
        #############################################
        centrality_genes_all_for_venn[centrality] = set(centality_sig_genes_in_hidden)
        centrality_genes_6hr_hidden_for_venn[centrality] = set(sig_genes_6hr)&set(centality_sig_genes_in_hidden)
        centrality_genes_24hr_hidden_for_venn[centrality] = set(sig_genes_24hr)&set(centality_sig_genes_in_hidden)
        #########################################

    for key, value in centrality_genes_in_hidden_dict.items():
        print(key, len(value), value[:5])
    print("common:", len(common_genes))

    import toolbox.visual_utilities as vu
    from venn import venn
    import matplotlib.pyplot as plt

    # vu.draw_venn_3group(centrality_genes_in_hidden_dict['RWR'],centrality_genes_in_hidden_dict['eigen'],centrality_genes_in_hidden_dict['degree'],group_labels=['RWR','eigen','degree'],
    #                     save_addr= "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Sig.Genes/keygenes_RWR_eigen_degree.png")

    ###########################

    venn(centrality_genes_all_for_venn)
    plt.savefig("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Sig.Genes/keygenes_RWR_eigen_degree_bw.pdf")
    plt.show()
    venn(centrality_genes_24hr_hidden_for_venn)
    plt.savefig(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Sig.Genes/24hr_keygenes_RWR_eigen_degree_bw.pdf")
    plt.show()
    venn(centrality_genes_6hr_hidden_for_venn)
    plt.savefig(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Sig.Genes/6hr_keygenes_RWR_eigen_degree_bw.pdf")
    plt.show()
    # print(centrality_results_6hr)

def sort_hidden_by_localization(covid_network_hidden_keyGene, localization_bands):
    key_gene_localization_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_moa_based_som_reverse_subcellular_conf4_v3.tsv",
                                           sep='\t')
    key_gene_localization_groupby = key_gene_localization_df.groupby('Protein')['Subcellular'].agg(list)
    print(key_gene_localization_groupby)
    subcellular = set(key_gene_localization_df['Subcellular'].to_list())
    print(subcellular-set(localization_bands))

    covid_network_hidden_keyGene.sort()
    print(covid_network_hidden_keyGene[:5])
    covid_network_hidden_keyGene_copy = covid_network_hidden_keyGene[:]

    covid_network_hidden_keyGene_localizaition_nucleus_ER = []
    covid_network_hidden_keyGene_localizaition_nucleus_no_ER = []
    covid_network_hidden_keyGene_localizaition_nucleus_with_ER = []
    covid_network_hidden_keyGene_localizaition_ER = []
    covid_network_hidden_keyGene_localizaition_else = []
    for localization in localization_bands[:1]:
        for gene in covid_network_hidden_keyGene_copy:

            if (localization in key_gene_localization_groupby[gene]) and (
                    localization_bands[1] not in key_gene_localization_groupby[gene]):
                covid_network_hidden_keyGene_localizaition_nucleus_no_ER.append(gene)
            elif (localization in key_gene_localization_groupby[gene]) and (
                    localization_bands[1] in key_gene_localization_groupby[gene]):
                covid_network_hidden_keyGene_localizaition_nucleus_with_ER.append(gene)


        covid_network_hidden_keyGene_copy = list(set(covid_network_hidden_keyGene_copy)-
                                                 set(covid_network_hidden_keyGene_localizaition_nucleus_no_ER)-
                                                 set(covid_network_hidden_keyGene_localizaition_nucleus_with_ER))
        covid_network_hidden_keyGene_copy.sort()

    for localization in localization_bands[1:2]:
        for gene in covid_network_hidden_keyGene_copy:
            if localization in key_gene_localization_groupby[gene]:
                covid_network_hidden_keyGene_localizaition_ER.append(gene)

        covid_network_hidden_keyGene_copy = list(set(covid_network_hidden_keyGene_copy)-set(covid_network_hidden_keyGene_localizaition_ER))
        covid_network_hidden_keyGene_copy.sort()

    covid_network_hidden_keyGene_localizaition_nucleus_ER += covid_network_hidden_keyGene_localizaition_nucleus_no_ER
    covid_network_hidden_keyGene_localizaition_nucleus_ER += covid_network_hidden_keyGene_localizaition_nucleus_with_ER
    covid_network_hidden_keyGene_localizaition_nucleus_ER += covid_network_hidden_keyGene_localizaition_ER

    for localization in localization_bands[2:]:
        for gene in covid_network_hidden_keyGene_copy:
            if localization in key_gene_localization_groupby[gene]:
                covid_network_hidden_keyGene_localizaition_else.append(gene)

        covid_network_hidden_keyGene_copy = list(set(covid_network_hidden_keyGene_copy)-set(covid_network_hidden_keyGene_localizaition_else))
        covid_network_hidden_keyGene_copy.sort()

    covid_network_hidden_keyGene_localizaition_else += covid_network_hidden_keyGene_copy
    covid_network_hidden_sorted = []
    covid_network_hidden_sorted.append(covid_network_hidden_keyGene_localizaition_nucleus_ER)
    covid_network_hidden_sorted.append(covid_network_hidden_keyGene_localizaition_else)
    return covid_network_hidden_sorted


def sort_hidden_by_enrichedPaths(covid_network_hidden_keyGene):
    hidden_enrichedPaths_dict = dh.load_obj("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/data/hidden_to_keyPaths_dict")

    hidden_sorted_dict = {}
    hidden_met_rna = []
    hidden_immune = []
    hidden_cellcycle = []
    hidden_others = []
    covid_network_hidden_keyGene_copied = covid_network_hidden_keyGene[:]
    for gene in covid_network_hidden_keyGene_copied:
        if gene in hidden_enrichedPaths_dict["Met.RNA"]:
            hidden_met_rna.append(gene)
        elif gene in hidden_enrichedPaths_dict["CellCycle"]:
            hidden_cellcycle.append(gene)
        elif gene in hidden_enrichedPaths_dict["Immune"]:
            hidden_immune.append(gene)
        else:
            hidden_others.append(gene)
    hidden_sorted_dict["Met.RNA"] = sorted(hidden_met_rna,reverse=True)
    hidden_sorted_dict["CellCycle"] = sorted(hidden_cellcycle,reverse=True)
    hidden_sorted_dict["Immune"] = sorted(hidden_immune, reverse=True)
    hidden_sorted_dict["Other"] = sorted(hidden_others,reverse=True)

    covid_network_hidden_sorted = []
    for label, genes in hidden_sorted_dict.items():
        print(label, genes)
        for gene in genes:
            covid_network_hidden_sorted.append([label,gene])


    return covid_network_hidden_sorted

def sort_hidden_by_centrality_rwr(gene_set, centrality_dict):

    sorted_gene_set = []

    gene_set_eigen = list(set(gene_set)&set(centrality_dict['key_eigen']))
    gene_set_degree = list(set(gene_set)&set(centrality_dict['key_degree']))
    gene_set_bw = list(set(gene_set)&set(centrality_dict['key_bw']))
    gene_set_rwr = list(set(gene_set)&set(centrality_dict['key_rwr']))
    gene_set_other = list(set(gene_set)-set(gene_set_eigen)-set(gene_set_degree)-set(gene_set_bw)-set(gene_set_rwr))

    sorted_gene_set = sorted(gene_set_eigen,reverse=True) \
                      + sorted(gene_set_degree,reverse=True) \
                      + sorted(gene_set_bw,reverse=True) \
                      + sorted(gene_set_rwr,reverse=True) \
                      + sorted(gene_set_other,reverse=True)

    return sorted_gene_set


def sort_hidden_by_enrichedPaths_centrality(covid_network_hidden_keyGene, centrality_dict):
    hidden_enrichedPaths_dict = dh.load_obj("../Data/hidden_to_keyPaths_dict")

    hidden_sorted_dict = {}
    hidden_met_rna = []
    hidden_immune = []
    hidden_cellcycle = []
    hidden_met_proteins = []
    hidden_others = []
    covid_network_hidden_keyGene_copied = covid_network_hidden_keyGene[:]
    for gene in covid_network_hidden_keyGene_copied:
        if gene in hidden_enrichedPaths_dict["Met.RNA"]:
            hidden_met_rna.append(gene)
        elif gene in hidden_enrichedPaths_dict["CellCycle"]:
            hidden_cellcycle.append(gene)
        elif gene in hidden_enrichedPaths_dict["Immune"]:
            hidden_immune.append(gene)
        elif gene in hidden_enrichedPaths_dict["Met.Proteins"]:
            hidden_met_proteins.append(gene)
        else:
            hidden_others.append(gene)

    hidden_sorted_dict["Met.RNA"] = sort_hidden_by_centrality_rwr(hidden_met_rna, centrality_dict)
    hidden_sorted_dict["CellCycle"] = sort_hidden_by_centrality_rwr(hidden_cellcycle, centrality_dict)
    hidden_sorted_dict["Immune"] = sort_hidden_by_centrality_rwr(hidden_immune, centrality_dict)
    hidden_sorted_dict["Met.Proteins"] = sort_hidden_by_centrality_rwr(hidden_met_proteins, centrality_dict)
    hidden_sorted_dict["Other"] = sort_hidden_by_centrality_rwr(hidden_others, centrality_dict)

    covid_network_hidden_sorted = []
    for label, genes in hidden_sorted_dict.items():
        print(label, genes)
        for gene in genes:
            covid_network_hidden_sorted.append([label,gene])


    return covid_network_hidden_sorted



def make_backbone_key_genes_file(virus, key_gene_SARS, circos_node_file_a):

    graph_SARS = dh.load_obj(f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph")

    covid_network_all_nodes = list(set(graph_SARS.nodes()))
    covid_network_all_nodes = list(set(covid_network_all_nodes))

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

    centrality_dict = {
        "key_genes": key_genes,
        "key_eigen" : key_genes_eigen,
        "key_degree" : key_genes_degree,
        "key_bw" : key_genes_bw,
        "key_rwr" : key_genes_rwr
    }

    print("key_target",len(key_genes), len(key_genes_eigen), len(key_genes_degree), len(key_genes_bw), len(key_genes_rwr))
    dip_df = pd.read_csv(f"../Data/DIP/{virus}_DIP_no_duple.csv")
    dip_protein = pd.read_csv(f"../Data/DIP/{virus}_DIP_no_duple.csv")['gene_name'].tolist()
    covid_network_dip = list(set(covid_network_all_nodes)&set(dip_protein))

    covid_network_dip_keyGene = list(set(covid_network_dip)&set(key_genes))
    covid_network_dip_keyGene_eigen = list(set(covid_network_dip_keyGene)&set(key_genes_eigen))
    covid_network_dip_keyGene_degree = list(set(covid_network_dip_keyGene)&set(key_genes_degree))
    covid_network_dip_keyGene_bw = list(set(covid_network_dip_keyGene) & set(key_genes_bw))
    covid_network_dip_keyGene_rwr = list(set(covid_network_dip_keyGene) & set(key_genes_rwr))

    # print(len(covid_network_dip_keyGene),len(covid_network_dip_keyGene_degree),len(covid_network_dip_keyGene_bw), len(covid_network_dip_keyGene_eigen))


    dep_SARS_protein = pd.read_csv(f"../Data/DEP/{virus}_DEP.csv")['Gene_name'].tolist()

    dep_protein = list(set(dep_SARS_protein)-set(dip_protein))  #

    covid_network_dep = list(set(covid_network_all_nodes)&set(dep_protein))
    covid_network_dep_keyGene = list(set(covid_network_dep)&set(key_genes) - set(covid_network_dip_keyGene))
    covid_network_dep_keyGene_eigen = list(set(covid_network_dep_keyGene) & set(key_genes_eigen))
    covid_network_dep_keyGene_degree = list(set(covid_network_dep_keyGene) & set(key_genes_degree))
    covid_network_dep_keyGene_bw = list(set(covid_network_dep_keyGene) & set(key_genes_bw))
    covid_network_dep_keyGene_rwr = list(set(covid_network_dep_keyGene) & set(key_genes_rwr))

    covid_network_hidden = list(set(covid_network_all_nodes)-set(dep_protein)-set(dip_protein))
    covid_network_hidden_keyGene = list(set(covid_network_hidden)&set(key_genes))
    # print(len(covid_network_all_nodes))

    # print(len(key_genes), len(key_genes))

    print(len(covid_network_dip_keyGene))
    print(len(covid_network_dep_keyGene))
    print(len(covid_network_hidden_keyGene))



    ##########
    ## DIP to sub structure
    #########

    # moa_to_pathway_groupby_df = moa_to_pathway_df.groupby('MoA_Category')['Pathways'].agg(list)
    # print(dip_df)
    dip_preygene_df = dip_df.groupby('bait_name')['gene_name'].agg(list)

    print(dip_preygene_df)
    dip_bait_gene_dict = dict()
    for key,value in dip_preygene_df.items():
        # print(key)
        # string_gene = dip_stringgene_df[key]
        # bait_gene = list(set(value).union(set(string_gene)))
        dip_bait_gene_dict[key] = value

    dip_bait_gene_dict_key_sorted = sorted(list(dip_bait_gene_dict.keys()),key=natural_keys)
    dip_bait_gene_dict_sorted = dict()
    for key in dip_bait_gene_dict_key_sorted:
        # print(key)
        dip_bait_gene_dict_sorted[key] = dip_bait_gene_dict[key]

    # print(dip_bait_gene_dict_sorted.keys())
    covid_network_dip_node_with_structure = []
    covid_network_dip_coreGene_backup = covid_network_dip_keyGene[:]
    for structure,bait_gene in dip_bait_gene_dict_sorted.items():
        key_genes_in_structrure = list(set(covid_network_dip_keyGene)&set(bait_gene))
        #######################
        # for sort
        #######################
        structure_degree_gene = list(set(covid_network_dip_keyGene_degree)&set(key_genes_in_structrure))
        structure_bw_gene = list(set(covid_network_dip_keyGene_bw)&set(key_genes_in_structrure))
        structure_eigen_gene = list(set(covid_network_dip_keyGene_eigen)&set(key_genes_in_structrure))
        structure_rwr_gene = list(set(covid_network_dip_keyGene_rwr)&set(key_genes_in_structrure))

        structure_degree_gene.sort()
        structure_bw_gene.sort()
        structure_eigen_gene.sort()
        structure_rwr_gene.sort()

        for gene in structure_eigen_gene:
            covid_network_dip_node_with_structure.append([structure,gene])
        for gene in structure_degree_gene:
            covid_network_dip_node_with_structure.append([structure, gene])
        for gene in structure_bw_gene:
            covid_network_dip_node_with_structure.append([structure, gene])
        for gene in structure_rwr_gene:
            covid_network_dip_node_with_structure.append([structure, gene])
        ########################
        # covid_network_dip_keyGene = list(set(covid_network_dip_keyGene) - set(key_genes_in_structrure))
    # print(dip_bait_gene_dict.keys())
    # print(covid_network_dip_node_with_structure)
    # print(len(covid_network_dip_node_with_structure))

    # covid_network_dep_node = [['DEP',x] for x in covid_network_dep_keyGene]
    ################
    # sort dep genes
    ################

    covid_network_dep_keyGene_degree.sort()
    covid_network_dep_keyGene_bw.sort()
    covid_network_dep_keyGene_eigen.sort()
    covid_network_dep_keyGene_rwr.sort()
    covid_network_dep_keyGene_sorted = covid_network_dep_keyGene_eigen+covid_network_dep_keyGene_degree+covid_network_dep_keyGene_bw + covid_network_dep_keyGene_rwr
    covid_network_dep_node = [['DEP', x] for x in covid_network_dep_keyGene_sorted]


    ###################
    ## HIDDEN
    ##############
    # HIDDEN enriched pathways
    ##############

    covid_network_hidden_sorted = sort_hidden_by_enrichedPaths_centrality(covid_network_hidden_keyGene, centrality_dict)

    covid_network_background_node = []

    covid_network_background_node += covid_network_dip_node_with_structure
    covid_network_background_node += covid_network_hidden_sorted
    covid_network_background_node += covid_network_dep_node


    ############################################

    covid_network_background_df = pd.DataFrame(covid_network_background_node, columns=['Mode','Gene'])
    print(covid_network_background_df)
    covid_network_background_df.to_csv( circos_node_file_a,sep='\t',index=False, header=False)
#

def make_key_edges(network_edge_list, key_gene_list):
    key_edges = []

    for edge in network_edge_list:
        if len(set(edge)&set(key_gene_list)) == 2:
            key_edges.append(edge)

    print(len(key_edges))
    return key_edges

def read_band_name(band_name_a):
    with open(band_name_a) as band_name_f:
        band_name_data = [x.strip().split('\t') for x in band_name_f.readlines()]

    band_name_dict = dict()
    for band in band_name_data:
        band_name_dict[band[1]] = band[0]

    return band_name_dict


def write_edges(write_addr , links):
    links = ['\t'.join(x) for x in links]
    with open(write_addr,'w') as link_write_f:
        link_write_f.write('\n'.join(links))



def get_selected_drug_targets(selected_drugs_df, selected_pathways):
    reactome_gmt_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Reactome/ReactomePathways_for_gProfiler.gmt"

    reactome_gene_dict = dict()

    with open(reactome_gmt_addr) as reactome_gmt_f:
        reactome_path_gene = [x.strip().split('\t') for x in reactome_gmt_f.readlines()]

    for path in reactome_path_gene:
        reactome_gene_dict[path[1]] = path[2:]

    selected_pathway_proteins = []

    for path in selected_pathways:
        path_genes =reactome_gene_dict[path]
        selected_pathway_proteins += path_genes[:]

    selected_pathway_proteins = list(set(selected_pathway_proteins))
    # print(len(selected_pathway_proteins))
    # print(selected_pathway_proteins)

    selected_drugtargets_6hr = get_drug_targets_in_network('6hr',selected_drugs_df)
    selected_drugtargets_24hr = get_drug_targets_in_network('24hr',selected_drugs_df)
    selected_drugtargets_common = get_drug_targets_in_network('Common',selected_drugs_df)

    selected_drugtargets_6hr = list(set(selected_drugtargets_6hr)&set(selected_pathway_proteins))
    selected_drugtargets_24hr = list(set(selected_drugtargets_24hr) & set(selected_pathway_proteins))
    selected_drugtargets_common = list(set(selected_drugtargets_common) & set(selected_pathway_proteins))

    selected_drugtargets_only_6hr = list(set(selected_drugtargets_6hr) - set(selected_drugtargets_common))
    selected_drugtargets_only_24hr = list(set(selected_drugtargets_24hr) - set(selected_drugtargets_common)- set(selected_drugtargets_6hr))
    # print(len(selected_drugtargets_only_6hr),len(selected_drugtargets_only_24hr),len(selected_drugtargets_common))

    return selected_drugtargets_only_6hr+selected_drugtargets_only_24hr+selected_drugtargets_common



def main():


    # virus_list = ['SARS-CoV', 'SARS-CoV-2']
    for virus in ['SARS-CoV', 'SARS-CoV-2']:

        #######################
        ## Step 0: make backbone list
        #######################
        key_gene_SARS_df = pd.read_csv(f"../result/{virus}/Sig.Genes/{virus}_key_protein_every.txt")
        key_gene_SARS = key_gene_SARS_df['Gene'].to_list()
        # print(key_gene_SARS)
        # key_gene_SARS.append('COVID19')


        circos_node_file_a = f"../result/{virus}/Circos/data/{virus}_backbone_node_info_key_genes_high_level_paths.tsv"

        # make_backbone_key_genes_file(virus, key_gene_SARS, circos_node_file_a)

        # ###################
        # # Step 1: make circos node list
        # ##################
        circos_data_file_a = f"../result/{virus}/Circos/data/{virus}_DIP_structure_DEP_HIDDEN_key_genes_high_level_paths.txt"

        # make_circos_data_file_for_covid(virus, circos_data_file_a, circos_node_file_a)

        ########################
        ## Step 1-1 : make color file
        ##########################
        circos_color_file_a = f"../result/{virus}/Circos/data/{virus}_background_color_key_genes_high_level_paths_dip_structure_in_network.txt"
        # make_color_file(virus, circos_data_file_a, circos_color_file_a)
        #
        #
        #################
        # Step 2-1: make key target edge
        #################

        graph_SARS = dh.load_obj(f"../result/{virus}/network/{virus}_All_Structure_All_Shortest_Paths_Graph")


        network_edge_SARS = list(graph_SARS.edges())

        network_node_SARS = list(graph_SARS.nodes())

        network_key_edge_SARS = make_key_edges(network_edge_SARS, key_gene_SARS)

        network_key_edge_SARS_addr = f"../result/{virus}/Circos/data/{virus}_network_key_links_high_level_paths_with_symbol.txt"

        print(network_key_edge_SARS[:3])
        write_edges(network_key_edge_SARS_addr, network_key_edge_SARS)

        circos_SARS_key_link_file_a = f"../result/{virus}/Circos/data/{virus}_network_key_links_high_level_paths_dip_structure.txt"


        make_link_file(network_key_edge_SARS, circos_data_file_a, circos_SARS_key_link_file_a)



def main_for_revision():
    '''
    For NBT Revision
    :return:
    '''

    #####################
    # FOr NBT rebuttal, make symple links
    # Oct 24
    ######################

    circos_node_file_a = '/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/network_backbone_node_info_key_genes_high_level_paths_v2.tsv'

    network_key_edge_6hr_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_6hr_network_key_links_high_level_paths_with_symbol_v2.txt"
    network_key_edge_24hr_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/COVID_24hr_network_key_links_high_level_paths_with_symbol_v2.txt"

    circos_6hr_key_link_file_a = '/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/supple_COVID_6hr_network_key_links_in_circos.xlsx'
    circos_24hr_key_link_file_a = '/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Circos/data/supple_COVID_24hr_network_key_links_in_circos.xlsx'

    make_link_with_symbol_file(network_key_edge_6hr_addr, circos_node_file_a, circos_6hr_key_link_file_a)
    make_link_with_symbol_file(network_key_edge_24hr_addr, circos_node_file_a, circos_24hr_key_link_file_a)

if __name__ == '__main__':

    main()
    # main_for_revision()