# Created by woochanghwang at 23/07/2021
import pandas as pd

def main():
    virus_list = ["SARS-CoV", "SARS-CoV-2"]
    for virus in virus_list:
        if virus == "SARS-CoV": pvalue = "1e-13"
        elif virus == "SARS-CoV-2": pvalue = "1e-12"
        virus_top3_df = pd.read_excel(f"../result/{virus}/Drug/{virus}_extendAll_Reactome_low_level_{pvalue}_top3.xlsx")
        som_pathway_df = pd.read_excel(f"../result/{virus}/SOM/{virus}_pathway_SOM.xlsx")
        som_pathway_df = som_pathway_df[[virus,'SOM','Cluster']]
        virus_top3_df = virus_top3_df[['Top3','Top2']]
        virus_top3_df = virus_top3_df.drop_duplicates()
        print(virus_top3_df)
        print(som_pathway_df)

        som_pathway_top2_df = pd.merge(left = som_pathway_df,
                                       right = virus_top3_df,
                                       how = 'left',
                                       left_on= virus,
                                       right_on= 'Top3')

        print(som_pathway_top2_df)

        som_pathway_top2_df.to_excel(f"../result/{virus}/SOM/{virus}_pathway_SOM_Top2.xlsx")




if __name__ == '__main__':
    main()
