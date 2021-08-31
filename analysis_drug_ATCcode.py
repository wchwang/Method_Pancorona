# Created by woochanghwang at 28/06/2020
# Modified by woochang at 21/07/2021
'''
Count ATC code
Draw Bar graph
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def count_ATC_of_candidate_drugs():
    virus = "SARS-CoV-2"
    # candidate_drug_atc_df = pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP.CT.xlsx")
    candidate_drug_atc_df = pd.read_excel(f"../result/{virus}/Drug/{virus}_candidate_drug_info_target.SIP.xlsx")
    candidate_drug_atc_df = candidate_drug_atc_df[candidate_drug_atc_df['atc_codes']!='None']

    print(candidate_drug_atc_df.head())
    drug_atc_codes = candidate_drug_atc_df['atc_codes'].to_list()

    all_atc_codes = []
    for drug in drug_atc_codes:
        drug = str(drug).strip().split('|')
        all_atc_codes+= drug

    print(len(all_atc_codes))

    all_atc_codes = list(set(all_atc_codes)-set(['many']))

    all_atc_code_4 = [x[:4] for x in all_atc_codes]
    all_atc_code_3 = [x[:3] for x in all_atc_codes]

    print(all_atc_code_4[:3])
    print(all_atc_code_3[:3])

    all_atc_code_4_df = pd.DataFrame(all_atc_code_4, columns=['ATCcode4'])
    print(all_atc_code_4_df)
    all_atc_code_4_groupby_df = all_atc_code_4_df.groupby(['ATCcode4']).size().reset_index(name='counts')
    all_atc_code_4_groupby_df = all_atc_code_4_groupby_df.sort_values('counts', ascending=False)

    print(all_atc_code_4_groupby_df)

    all_atc_code_3_df = pd.DataFrame(all_atc_code_3, columns=['ATCcode3'])
    print(all_atc_code_4_df)
    all_atc_code_3_groupby_df = all_atc_code_3_df.groupby(['ATCcode3']).size().reset_index(name='counts')
    all_atc_code_3_groupby_df = all_atc_code_3_groupby_df.sort_values('counts', ascending=False)

    # all_atc_code_3_groupby_df = all_atc_code_3_groupby_df[all_atc_code_3_groupby_df['counts']>1]
    print(all_atc_code_3_groupby_df)
    # all_atc_code_3_groupby_df = all_atc_code_3_groupby_df.set_index('ATCcode3')

    ax = all_atc_code_3_groupby_df[['counts']].plot(kind='bar', title='ATC code', figsize = (15,10), legend=True, fontsize = 12)
    ax.set_xlabel("ATC Code", fontsize=14)
    ax.set_ylabel("Counts", fontsize=14)
    # plt.savefig("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/ATC/ATC_code_3.pdf")
    # plt.show()

    atc_code_desc = pd.read_excel("/Users/woochanghwang/PycharmProjects/MTIProject/COVID-19/result/Drug/COVID19_candidiate_v4_ATC.xlsx",sheet_name="ATC_second",names=['ATCcode3','Desc'])
    print(atc_code_desc)
    all_atc_code_3_groupby_df['ATCcode3'] = all_atc_code_3_groupby_df['ATCcode3'].str.strip()
    atc_code_desc['ATCcode3'] = atc_code_desc['ATCcode3'].str.strip()
    all_atc_code_3_groupby_desc_df = all_atc_code_3_groupby_df.merge(atc_code_desc,how='left')
    print(all_atc_code_3_groupby_desc_df)
    all_atc_code_3_groupby_desc_df = all_atc_code_3_groupby_desc_df.set_index('Desc')
    ax =  all_atc_code_3_groupby_desc_df[['counts']].plot(kind='bar', figsize = (40,25), legend=True, fontsize = 30, color = "#ff7f0e")
    ax.set_xlabel("ATC Code", fontsize=30)
    ax.set_ylabel("Counts", fontsize=30)
    plt.rcParams['font.family'] = 'Helvetica'
    plt.tight_layout()
    plt.savefig(f"../result/{virus}/Drug/ATC/{virus}_ATC_code_3_desc_candidate_drug_v2.pdf")
    plt.show()
    all_atc_code_3_groupby_desc_df.to_excel(f"../result/{virus}/Drug/ATC/{virus}_candidiate_ATC_grouby_v2.xlsx")


def check_ATC_drugs():
    candidate_drugs = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/ATC/DrugID.txt")['DrugID'].to_list()

    drugbank_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/v5.1.6/drugbank.tsv",sep='\t')

    candidate_drugbank_df = drugbank_df[drugbank_df['drugbank_id'].isin(candidate_drugs)]
    print(candidate_drugbank_df)

    candidate_drugbank_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/ATC/candidate_drugs_atc_all_info.xlsx")

def make_candidate_drugs_ATC_code():
    candidate_drug_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.v5.multimoa.category.subcategory.xlsx")
    # print(candidate_drug_df.columns)
    candidate_drugs = candidate_drug_df['Drug_ID'].to_list()
    drugbank_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/v5.1.6/drugbank.tsv",
                              sep='\t')
    candidate_drugbank_df = drugbank_df[drugbank_df['drugbank_id'].isin(candidate_drugs)]
    candidate_drugbank_df = candidate_drugbank_df[candidate_drugbank_df.atc_codes.notnull()]
    # candidate_drugbank_df =candidate_drugbank_df[candidate_drugbank_df['atc_codes']!= np.nan]

    # print(set(candidate_drugbank_df['atc_codes'].to_list()))
    candidate_drugs_in_drugbank = candidate_drugbank_df['drugbank_id'].to_list()
    candidate_not_drugbank_df = candidate_drug_df[~candidate_drug_df['Drug_ID'].isin(candidate_drugs_in_drugbank)]

    print(candidate_drugbank_df.columns)
    print(candidate_not_drugbank_df.columns)

    candidate_drugbank_simple_df = candidate_drugbank_df[['drugbank_id','name','atc_codes']]
    candidate_not_drugbank_simple_df = candidate_not_drugbank_df[['Drug_ID','Name','atc_codes']]

    candidate_drugbank_simple_df = candidate_drugbank_simple_df.rename(columns={'drugbank_id': 'Drug_ID',
                                                                                'name': 'Name'})

    candidate_drug_atc_df = pd.concat([candidate_drugbank_simple_df, candidate_not_drugbank_simple_df])
    print(candidate_drug_atc_df)

    candidate_drug_atc_df = candidate_drug_atc_df.sort_values(by='Drug_ID')
    candidate_drug_atc_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidiate_v5_ATC.xlsx", index=False)
def count_number_of_drugs_in_aATC(aATC):
    candidate_drug_atc_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidiate_v4_ATC.xlsx",
        sheet_name='v2')
    print(candidate_drug_atc_df)
    candidate_drug_aATC_df = candidate_drug_atc_df[candidate_drug_atc_df['atc_codes'].str.startWith(aATC)]
    print(len(candidate_drug_aATC_df))

def count_drugs_top_ATCs():
    ATCs_in_candidate_drugs_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/ATC/203drugs/COVID19_candidiate_v5_ATC_grouby_v1.xlsx")

    print(ATCs_in_candidate_drugs_df)
    candidate_drug_atc_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.v5.multimoa.category.subcategory.xlsx")
    candidate_drug_atc_df = candidate_drug_atc_df[candidate_drug_atc_df['atc_codes']!='None']
    atc_first_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidiate_v4_ATC.xlsx",sheet_name='ATC_first')
    print(atc_first_df)
    atc_first_list = atc_first_df['ATC'].to_list()

    atc_first_drug_dict = dict()
    number_of_drugs = len(candidate_drug_atc_df)
    for atc in atc_first_list:
        # candidate_drug_aATC_df = candidate_drug_atc_df[candidate_drug_atc_df['atc_codes'].str.contains(atc)]

        candidate_drug_aATC_df = candidate_drug_atc_df[(candidate_drug_atc_df['atc_codes'].str.startswith(atc)) | (candidate_drug_atc_df['atc_codes'].str.contains('\|'+atc))]
        drug_count = len(candidate_drug_aATC_df)
        percent = drug_count/number_of_drugs*100
        atc_first_drug_dict[atc] = [drug_count, percent]
        print(atc, drug_count, percent)

    # print(atc_first_drug_dict)
    atc_first_drug_df = pd.DataFrame.from_dict(atc_first_drug_dict,orient='index',columns=['Drugs','Percentage'])
    atc_first_drug_df = atc_first_drug_df.reset_index()
    atc_first_drug_df = atc_first_drug_df.merge(atc_first_df,left_on='index',right_on='ATC')
    atc_first_drug_df = atc_first_drug_df[['ATC','Desc','Drugs','Percentage']]
    atc_first_drug_df = atc_first_drug_df.sort_values('Drugs', ascending=False)
    atc_first_drug_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/ATC/203drugs/COVID19_candidate_v5_count_ATC_first.xlsx")

    atc_first_drug_df = atc_first_drug_df.set_index('Desc')

    ax = atc_first_drug_df[['Percentage']].plot(kind='bar', figsize=(30, 30), legend=True,
                                                         fontsize=30)
    ax.set_xlabel("ATC First level", fontsize=30)
    ax.set_ylabel("Drug Percentage", fontsize=30)
    plt.rcParams['font.family'] = 'Helvetica'
    plt.tight_layout()
    plt.savefig("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/ATC/240drugs/ATC_code_first_desc_drug_percentage.pdf")
    plt.show()




if __name__ == '__main__':
    count_ATC_of_candidate_drugs()
    # check_ATC_drugs()
    # # count_number_of_drugs_in_aATC('L01')
    # count_drugs_top_ATCs()
    #
    # # make_candidate_drugs_ATC_code()   # to upgrade drubank5.1.5 --> 5.1.6