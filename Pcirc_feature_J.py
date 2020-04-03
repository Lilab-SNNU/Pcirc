import pandas as pd

def transform(seq):
    transform_info_dict = {'A': 1, 'C': 2, 'G': -2, 'T': -1, 'N': 0}
    seq_info = []
    feature_list = []

    for i in seq:
        seq_info.append(transform_info_dict[i])

    df = pd.DataFrame(seq_info)
    df = df.T
    array = df.values[:, :]

    for i in array:
        feature_list.append(list(i))

    return feature_list


def feature_J(file_name):
    '''if '.txt' in file_name:
        res_name = file_name.replace('.txt', '.csv')
    else:
        res_name = file_name + '.csv'
    '''

    file = open(file_name, 'r')
    # res_file = open(res_name, 'w')
    file_info = file.read().split('>')[1:]
    file_info_dict = {}
    transform_info_dict = {'A': 1, 'C': 2, 'G': -2, 'T': -1, 'N': 0}
    index_list = []

    for i in range(200):
        if i < 100:
            index = 'up_site' + str(i)
            index_list.append(index)
        elif i >= 100:
            index = 'down_site' + str(i - 100)
            index_list.append(index)

    n = 0
    header_list = []
    for line in file_info:
        n = n + 1
        file_info_list = []
        header, seq = '>' + line.split('\n')[0], line.split('\n')[1]

        for i in range(len(seq)):
            feature = transform_info_dict[seq[i]]
            file_info_list.append(feature)

        header_list.append(header)
        file_info_dict[header] = file_info_list

    df = pd.DataFrame(file_info_dict, index=index_list)
    df = df.T
    df.sort_index(inplace=True)

    return df