# Created by woochanghwang at 12/07/2021

import toolbox.drugbank_handler_postgres as drug_pg
import pandas as pd


def get_DDI_of_a_drug(drug):

    # drug_a = row['DrugA_ID']
    # drug_b = row['DrugB_ID']
    sql = "SELECT * from structured_drug_interactions where subject_drug_drugbank_id = \'{}\' " \
          "or affected_drug_drugbank_id = \'{}\'".format(drug, drug)
    print(sql)
    sql_result= drug_pg.query_postgres_get_DF(sql,**drug_pg.DRUGBANK_DB_INFO)

    return sql_result

def get_CT_of_a_condition(condition):

    # drug_a = row['DrugA_ID']
    # drug_b = row['DrugB_ID']
    sql = f"SELECT * from drug_clinicaltrial_2 where condition_id = \'{condition}\' "

    print(sql)
    sql_result= drug_pg.query_postgres_get_DF(sql,**drug_pg.DRUGBANK_DB_INFO)

    return sql_result

def add_CT_on_simulation_result():
    candidate_drugs_CT_df = pd.read_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_CT_v2.xlsx",sheet_name="CT" )
    print(candidate_drugs_CT_df)
    CT_drugs = candidate_drugs_CT_df['drugbank_id'].tolist()
    # drug_info_df = pd.read_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_target.SIP.xlsx")
    drug_info_df = pd.read_excel("../result/SARS-CoV/Drug/SARS-CoV_candidate_drug_info_target.SIP.xlsx")

    # drug_info_df['Clincical_Trial'] = drug_info_df
    drug_info_df['Clincical_Trial'] = drug_info_df['source'].isin(CT_drugs)
    print(drug_info_df)
    # drug_info_df.to_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_target.SIP.CT.v2.xlsx", index=False)
    drug_info_df.to_excel("../result/SARS-CoV/Drug/SARS-CoV_candidate_drug_info_target.SIP.CT.v2.xlsx", index=False)

def main():
    # condition = "DBCOND0129755"
    #
    # CT_result = get_CT_of_a_condition(condition)
    # # print(CT_result)
    # # CT_result.to_excel("../Data/SARS-CoV-2_CT_drugs_v2.xlsx", index=False)
    #
    # candidate_drugs_df = pd.read_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_drugbank.xlsx")
    # # candidate_drugs_df = pd.read_excel("../result/SARS-CoV-2/Drug/")
    # print(candidate_drugs_df)
    #
    # candidate_drugs_CT_df = pd.merge(left=candidate_drugs_df,
    #                                  right = CT_result,
    #                                  how='left',
    #                                  left_on='source',
    #                                  right_on='drugbank_id')
    # print(candidate_drugs_CT_df)
    # # candidate_drugs_CT_df.to_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_CT.xlsx", index=False)
    # candidate_drugs_CT_df.to_excel("../result/SARS-CoV-2/Drug/SARS-CoV-2_candidate_drug_info_CT_v2.xlsx", index=False)


    add_CT_on_simulation_result()
if __name__ == '__main__':
    main()
