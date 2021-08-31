# Created by woochanghwang at 12/07/2021

# Created by woochanghwang at 01/05/2020

import scipy.stats as stats

def main():
    '''
                    Candidate drugs Not_candidate
    COVID-19        44              209
    Non COVID-19    152             1512
    :return:
                    Candidate drugs Not_candidate
    COVID-19        45              257
    Non COVID-19    151             1464
    :return:

    '''

    '''
                Candidate drugs Not_candidate
    SARS        21              232
    Non SARS    81             1583
    :return:
    '''


# oddratio, pvalue = stats.fisher_exact([[21,118],[225,5605]])
    # oddratio, pvalue = stats.fisher_exact([[26, 92], [220, 3163]])  #approved + trial
    # oddratio, pvalue = stats.fisher_exact([[29, 121], [217, 2931]])  # approved + trial
    # oddratio, pvalue = stats.fisher_exact([[24, 88], [222, 1755]])  # approved

    ## Clinical trial Only 0520
    # oddratio, pvalue = stats.fisher_exact([[28, 121], [172, 2797]])  # approved + trial 06/26
    # print(pvalue)
    ## Clinical trial Only- 0723
    # oddratio, pvalue = stats.fisher_exact([[40, 232], [160, 2581]])  # approved + trial 06/26

    # # Clinical trial + Literatures
    # oddratio, pvalue = stats.fisher_exact([[72, 78], [168, 2840]])  # approved + trial 06/26
    # print(pvalue)

    ## Approved only
    # oddratio, pvalue = stats.fisher_exact([[44, 209], [152, 1514]])  # approved + trial 06/26


    oddratio, pvalue = stats.fisher_exact([[44, 152], [209, 1512]])  # covid-19

    # oddratio, pvalue = stats.fisher_exact([[21, 81], [232, 1583]])  # SARS

    # oddratio, pvalue = stats.fisher_exact([[45, 151], [257, 1464]])  # covid-19_v2


    print(pvalue)
if __name__ == '__main__':
    main()