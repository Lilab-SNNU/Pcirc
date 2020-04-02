import bisect
import pandas as pd

def get_orf(filename):
    file = open(filename, 'r')
    info_all = file.read().split('//')[:-1]
    orf_length_dict = {}
    gene_names, gene_name_list = [], []
    counts_num = 0

    for i in info_all:
        info_list = i.split('\n')[:-1]
        orf_length_list = []
        counts_num = counts_num + 1

        for info in info_list:

            if 'LOCUS' in info:
                mid_str = '/'.join(info.split())
                gene_name, gene_length = '>' + mid_str.split('/')[1], int(mid_str.split('/')[2])/2
                gene_names.append(gene_name)

            elif '/dna_len=' in info:
                orf_length = int(info.split('=')[-1])
                orf_length_list.append(orf_length)

            n, l = len(orf_length_list), []
            index = bisect.bisect(orf_length_list, gene_length)

            if n > 0 and index > 0:
                orf_max = orf_length_list[index - 1]
                orf_coverage = orf_max / gene_length
                l.append(orf_max)
                l.append(orf_coverage)
                orf_length_dict[gene_name] = l

            else:
                orf_max = 0
                orf_coverage = orf_max / gene_length
                l.append(orf_max)
                l.append(orf_coverage)
                orf_length_dict[gene_name] = l

    df_orf = pd.DataFrame(orf_length_dict, index=['orf_max', 'orf_coverage'])
    df_orf = df_orf.T
    df_orf.sort_index(inplace=True)

    return df_orf