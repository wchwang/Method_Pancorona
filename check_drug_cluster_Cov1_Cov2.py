# Created by woochanghwang at 21/07/2021
import pandas as pd

def main():
    cov1_df = pd.read_excel("../result/SARS-CoV/Drug/SARS-CoV_candidate_drug_info_target.SIP.xlsx")
    cov2_df = pd.read_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_target.SIP.CT.xlsx")

    cov1_drug_cluster = cov1_df[['source','name','Cluster']]
    cov2_drug_cluster = cov2_df[['source','name','Cluster']]

    cov1_drug_cluster = cov1_drug_cluster.rename(columns={'source': 'source_cov1',
                                                          'name': 'name_cov1',
                                                          'Cluster': 'Cluster_cov1'})

    cov1_cov2_cluster = pd.merge(left=cov1_drug_cluster,
                                 right = cov2_drug_cluster,
                                 how="outer",
                                 left_on = 'source_cov1',
                                 right_on='source')

    print(cov2_drug_cluster)
    print(cov1_cov2_cluster)

    cov1_cov2_cluster.to_excel("../result/Cov1_Cov2_cluster_compare_v2.xlsx",index=False)


if __name__ == '__main__':
    main()
