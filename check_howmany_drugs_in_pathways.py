# Created by woochanghwang at 16/07/2021
import pandas as pd

def main():
    # threshold = 5.0e-2
    # threshold = 1.0e-12
    # virus = "SARS-CoV-2"
    virus = "SARS-CoV"
    threshold = 1.0e-13
    # fi_matrix_df = pd.read_csv( f"../result/SARS-CoV/Drug/F1/{virus}_drug_to_low_level_pathways_{threshold}_matrix.csv", index_col=0)
    fi_matrix_df = pd.read_csv( f"../result/{virus}/Drug/F1/{virus}_extendAll_drug_to_low_level_pathways_{threshold}_matrix.csv", index_col=0)
    # fi_matrix_df = pd.read_csv(f"../result/F1/Trypanosoma_drug_to_low_level_pathways_{threshold}_matrix_900.csv", index_col=0)
    fi_matrix_df_T = fi_matrix_df.T
    # fi_matrix_df['count_path'] = fi_matrix_df.iloc[:,:1559]
    fi_matrix_df_T['count_drugs'] = (fi_matrix_df_T>0).sum(1)
    print(fi_matrix_df_T)
    drug_list = list(fi_matrix_df_T)
    len_drugs = len(drug_list)-1
    # fi_matrix_df_T.to_excel( f"../result/SARS-CoV/Drug/F1/{virus}_drug_to_low_level_pathways_{threshold}_matrix_T.csv")
    fi_matrix_df_T.to_csv(f"../result/{virus}/Drug/F1/{virus}_extendAll_drug_to_low_level_pathways_{threshold}_matrix_T.csv")
    # fi_matrix_df_T_sub = fi_matrix_df_T[(fi_matrix_df_T['count_drugs'] > 23) ]
    # fi_matrix_df_T_sub = fi_matrix_df_T[(fi_matrix_df_T['count_drugs'] > (len_drugs/4)) & (fi_matrix_df_T['count_drugs'] < 127)]    #SARS-CoV-2
    fi_matrix_df_T_sub = fi_matrix_df_T[(fi_matrix_df_T['count_drugs'] > (len_drugs*0.3)) & (fi_matrix_df_T['count_drugs'] < (len_drugs * 0.6))]    #SARS-CoV
    # fi_matrix_df_T_sub = fi_matrix_df_T[(fi_matrix_df_T['count_drugs'] > 50) & (fi_matrix_df_T['count_drugs'] < 127)]
    fi_matrix_df_sub = fi_matrix_df_T_sub.T
    fi_matrix_df_sub['count_path'] = (fi_matrix_df_sub>0).sum(1)
    print(fi_matrix_df_sub[fi_matrix_df_sub['count_path']==0])
    print(fi_matrix_df_sub[fi_matrix_df_sub['count_path']==0].index)
    count_path = fi_matrix_df_sub['count_path'].tolist()
    count_path.sort()
    print(count_path)
    # f1_matrix_df = fi_matrix_df_sub.drop(columns = ['count_path'])
    # f1_matrix_df = f1_matrix_df.drop(axis=1, index = 'count_drugs')
    # print(f1_matrix_df)
    #
    # candidate_drug_df = pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP.xlsx")
    # print(list(candidate_drug_df))
    # candidate_drug_df = candidate_drug_df[['source','name']]
    # f1_matrix_df = pd.merge(left = candidate_drug_df,
    #                         right = f1_matrix_df,
    #                         how='left',
    #                         left_on='source',
    #                         right_index=True)
    # f1_matrix_df = f1_matrix_df.drop(columns=['source'])
    # f1_pathway = list(f1_matrix_df)
    # f1_matrix_df.to_csv(f"../result/{virus}/Drug/F1/{virus}_extendAll_drug_to_low_level_{len(f1_pathway[1:])}pathways_matrix.csv",index=False)
    #
    # f1_pathway = list(f1_matrix_df)
    # print(f1_pathway)
    # f1_pathway_df = pd.DataFrame(f1_pathway[1:],columns=['pathway'])
    # f1_pathway_df.to_csv(f"../result/{virus}/Drug/F1/{virus}_f1_pahtways.csv", index=False)




if __name__ == '__main__':
    main()
