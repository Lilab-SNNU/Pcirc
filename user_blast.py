import pandas as pd
def exactReads(filename):
    df = pd.read_table(filename, sep='\t', header=None)
    df = df[df.duplicated(subset=0, keep=False)]  
    df = df.drop_duplicates(subset=[0, 1, 3, 6, 7], keep=False)  
    return df

def cutDatafram(df):  # cut df based on chromnome
    set_chr, df_list = set(df[1].values), []  # get all chromnome
    for chr_num in set_chr:
        index_list = df[df[1] == chr_num].index.tolist()
        df_cut = df.loc[index_list, :9]
        df_cut = df_cut[df_cut.duplicated(subset=0, keep=False)]
        df_list.append(df_cut)
    return df_list

def getCircInformation(df, file):
    res = open(file, 'a')
    set_name = set(df[0].values)
    circ_start_end = []

    for name in set_name:
        index_list = df[df[0] == name].index.tolist()
        little_df = df.loc[index_list, :9]
        little_df.sort_values(6, inplace=True)
        for i in range(len(index_list) - 1):
            res_list = []
            if little_df.iloc[i, 1] == little_df.iloc[i+1, 1]:
                chr_name, querry_name = little_df.iloc[i, 1], name
                sum1 = little_df.iloc[i, 9] - little_df.iloc[i, 8]
                sum2 = little_df.iloc[i+1, 9] - little_df.iloc[i+1, 8]
                spacing = little_df.iloc[i + 1, 6]-little_df.iloc[i, 7]
                if (sum1 > 0 and sum2 > 0) and (little_df.iloc[i, 8] >= little_df.iloc[i+1, 9]) and (spacing in range(-3, 5)):
                    res_list.append(querry_name)
                    res_list.append(chr_name)
                    res_list.append('+')
                    res_list.append(little_df.iloc[i+1, 8])
                    res_list.append(little_df.iloc[i, 9])
                    circ_start_end.append(res_list)
                    res_list = []
                elif (sum1 < 0 and sum2 < 0) and (little_df.iloc[i, 8] <= little_df.iloc[i+1, 9]) and (spacing in range(-3, 5)):
                    res_list.append(querry_name)
                    res_list.append(chr_name)
                    res_list.append('-')
                    res_list.append(little_df.iloc[i, 9])
                    res_list.append(little_df.iloc[i+1, 8])
                    circ_start_end.append(res_list)
                    res_list = []
            else:pass

    res_df = pd.DataFrame(circ_start_end)
    res_df.to_csv(res, header=None, index=None)
    return res_df
