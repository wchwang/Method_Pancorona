# Created by woochanghwang at 12/07/2021

# Created by woochanghwang at 14/03/2021

from src_drug.p02_add_drug_info_on_result import add_drugbank_info, add_drug_target, groupby_drug_targets

def main():

    # ############################
    # ## Add drug name in result
    # ############################

    disease = "SARS-CoV"
    drug_proximity_f  = f"../result/{disease}/Drug/{disease}_drug_network_proximity.csv"
    result_file_prefix = f"../result/{disease}/Drug/{disease}_drug_network_proximity"
    # result_f = "{}.tsv".format(result_file_prefix)

    drugbank_data_file = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/Drugbank/v5.1.6/drugbank.tsv"
    drugbank_result_f = "{}_drugbank.tsv".format(result_file_prefix)
    add_drugbank_info(drug_proximity_f, drugbank_data_file, drugbank_result_f)


    # #############################
    # ## Add drug target in result
    # ############################
    drug_target_file = "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/STITCH/2020.12/9606.protein_chemical.links.v5.0.drugbank.v5.1.6.target_symbol.s900.onlyTarget.tsv"
    drug_target_result_f = "{}_target.tsv".format(result_file_prefix)
    add_drug_target(drugbank_result_f,drug_target_file, drug_target_result_f , type="symbol")


    # ########################
    # ## groupby drug targets
    # #######################

    drug_targets_addr = "{}_drug_to_target.csv".format(result_file_prefix)
    drug_groupby_target_addr = "{}_drug_groupby_target.tsv".format(result_file_prefix)

    groupby_drug_targets(drug_target_result_f, drug_targets_addr, drug_groupby_target_addr)




if __name__ == '__main__':
    main()
