import itertools as it
import pandas as pd

def combination(n):
    list_n, list_all = list(range(1, n + 1)), []
    for i in list_n:
        mylist = list(it.product(['A', 'C', 'G', 'T'], repeat=i))
        for j in mylist:
            string = ''.join(j)
            list_all.append(string)
    return list_all

def extractNuc(filename, n=4):
    if '.fa' in filename:
        outputfile = filename.replace('.fa', str(n) + '.csv')
    elif '.fasta' in filename:
        outputfile = filename.orflace('.fasta', str(n) + '.csv')
    else:
        outputfile = filename + '.csv'

    infile, outfile = open(filename, 'r'), open(outputfile, 'w')
    feature_list = combination(n)
    feature_N = {}
    features = []
    header_list = []

    for line in infile:
        if line[0] == '>':
            header = line.rstrip('\n')
            header_list.append(header)

        elif line[0] != '>':
            seq = line.rstrip('\n')
            length = len(seq)

            for feature in feature_list:
                N_len = len(feature)
                seq = seq + seq[0:N_len - 1]
                num = seq.count(feature)
                fre = num * N_len / length
                features.append(fre)

        feature_N[header] = features
        features = []

    df = pd.DataFrame(feature_N, index=feature_list)
    df = df.T
    df['gc'] = df.apply(lambda x: x['G'] + x['C'], axis=1)
    df.sort_index(inplace=True)

    return df