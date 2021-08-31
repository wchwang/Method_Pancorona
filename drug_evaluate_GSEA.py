# Created by woochanghwang at 11/08/2021
import pandas as pd
import gseapy as gp

def get_drugtargets(virus, cluster):
    if virus == "SARS-CoV":
        pathway_cluster_addr = "../result/SARS-CoV/Drug/SARS-CoV_candidate_drug_info_target.SIP.CT.v2.xlsx"
    elif virus == "SARS-CoV-2":
        pathway_cluster_addr = "../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_target.SIP.CT.v3.xlsx"

    drug_info_df = pd.read_excel(pathway_cluster_addr)
    drug_info_df_cluster = drug_info_df[drug_info_df['Cluster_N']==cluster]
    drug_targets = drug_info_df_cluster['Drug_Target_SIP'].tolist()

    drug_targets = ','.join(drug_targets)
    drug_targets = drug_targets.split(',')
    drug_targets = list(set(drug_targets))

    print(drug_targets)
    return drug_targets

def get_Reactome_GeneSet(virus, cluster, drug_targets):
    if virus == "SARS-CoV":
        pathway_cluster_addr = "../result/SARS-CoV/SOM/SARS-CoV_pathway_SOM_Top2.xlsx"
    elif virus == "SARS-CoV-2":
        pathway_cluster_addr = "../result/SARS-CoV-2/SOM/SARS-CoV-2_pathway_SOM_Top2_v3.xlsx"

    cov2_pathway_cluster = pd.read_excel(pathway_cluster_addr)
    print(cov2_pathway_cluster)
    cov2_pathway_cluster_c3 = cov2_pathway_cluster[cov2_pathway_cluster['Cluster']==cluster]
    c3_pathways = cov2_pathway_cluster_c3['SARS-CoV-2'].tolist()

    print(c3_pathways)
    c3_pathway_reactome = dict()

    reactome_gmt_gprofiler_addr =  "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/Reactome/ReactomePathways_for_gProfiler.gmt"
    with open(reactome_gmt_gprofiler_addr) as reactome_gmt_f:
        reactome_gmt_original = [x.strip().split('\t') for x in reactome_gmt_f.readlines()]

    reactome_gmt_gprofiler = dict()

    for path in reactome_gmt_original:
        # reactome_gmt_gprofiler.append([path[1],path[0],'\t'.join(path[2:])])
        reactome_gmt_gprofiler[path[1]] = path[2:]

    for path in c3_pathways:
        geneSet = reactome_gmt_gprofiler[path]
        c3_pathway_reactome[path] = list(set(geneSet)&set(drug_targets))

    print(c3_pathway_reactome)

    return c3_pathway_reactome

def get_FoldChange(virus):

    fc_df = pd.read_excel("../Data/Stukalov_data/Proteome_data.xlsx")
    print(fc_df)

    fc_df = fc_df[["Gene",virus]]
    fc_df = fc_df.sort_values(by=virus, ascending=True)

    return fc_df

def main():
    virus_list = ['SARS-CoV', 'SARS-CoV-2']
    virus = virus_list[1]

    cluster_drug_targets = get_drugtargets(virus,"C3")
    reactome_geneSet = get_Reactome_GeneSet(virus,"C3", cluster_drug_targets)
    print(reactome_geneSet)

    gene_FC = get_FoldChange(virus)

    corona_gsea_pre = gp.prerank(rnk=gene_FC,
                                 gene_sets=reactome_geneSet,
                                 processes=6,
                                 permutation_num=1000,
                                 min_size=5,
                                 outdir=f"../result/{virus}/GSEA",format='png',seed=6)

    corona_gsea_pre.res2d.to_excel(f"../result/{virus}/GSEA/{virus}_gsea_result.xlsx")




    # lusc_depmap_df = pd.read_csv(f"/Users/woochanghwang/PycharmProjects/MTIProject/LifeArc/Cancer_unbiased/LUSC/Depmap/broad_depmap_{cancer_type}.csv")
    # lusc_depmap_df = lusc_depmap_df[['Gene',cancer_type]]
    # # lusc_depmap_df[cancer_type] = lusc_depmap_df[cancer_type]*-1
    # lusc_depmap_df = lusc_depmap_df.sort_values(by=cancer_type,ascending=True)
    #
    # print(lusc_depmap_df)
    #
    # cc_gsea_pre = gp.prerank(rnk=lusc_depmap_df,
    #                          gene_sets=community_geneSet_dict,
    #                          processes=6,
    #                          permutation_num=1000,
    #                          min_size=10,
    #                          outdir=f"../LUSC/GSEA/{cancer_type}",format='png',seed=6)
    #
    # cc_gsea_pre.res2d.to_excel(f"../LUSC/GSEA/{cancer_type}/{cancer_type}_gsea_result_v3.xlsx")




if __name__ == '__main__':
    main()
