import argparse
import multiprocessing
import os
import pickle
import pandas as pd
import Pcirc_deal_blast_result
import Pcirc_circ_seq_extract
import Pcirc_feature_N
import Pcirc_feature_O
import Pcirc_feature_J
import Pcirc_double_circRNA

path = os.path.realpath(__file__)  # Get the install path of Pcirc
abs_dir = path[:path.rfind('/')]
parser = argparse.ArgumentParser(description='Predict plant circRNAs')

parser.add_argument("-i",
                    help='path to the result of blast, and the file format please refer to the format of '
                         'blast output format 6, you can refer to the test.blast')
parser.add_argument("-g",
                    help='path to genome, and the file format is fasta')
parser.add_argument("-m",
                    default=abs_dir +
                            "/model.pkl",
                    help='path to model, if you do not change the path of model, please do not set it')

args = parser.parse_args()  # Get all args
info_file = args.i  # Get the -i
genome_file = args.g  # Get the -g
model_filename = args.m  # Get the -m

# Pcirc running
alignment_res = '/'.join(info_file.split('/')[:-1]) + '/alignment_res'
df_alis = Pcirc_deal_blast_result.exactReads(info_file)
df_list = Pcirc_deal_blast_result.cutDatafram(df_alis)
pool = multiprocessing.Pool(len(df_list))
print('Alignment reads begining')

for dataframe in df_list:
    pool.apply_async(Pcirc_deal_blast_result.getCircInformation,
                     (dataframe, alignment_res, ))
pool.close()
pool.join()
print('Alignment reads ending')

junction_filename, circ_filename = Pcirc_circ_seq_extract.get_seq(alignment_res, genome_file)
res_filename = '/'.join(circ_filename.split('/')[:-1]) + '/circ'
circ_file = open(circ_filename, 'r')
res_file = open(res_filename, 'w')

# Ready for seq info
print('Ready for seq info')
info_circ = circ_file.read().split('>')[1:]
circ_dict = {}
for info in info_circ:
    header = '>' + info.split('\n')[0]
    seq = info.split('\n')[1]
    circ_dict[header] = seq

orf_p = Pcirc_double_circRNA.double(circ_filename)
orf_f = orf_p + '.gb'
os.system("ugene find-orfs --in=%s --out=%s" % (orf_p, orf_f))  # invoke the ugene

df_N, df_J, df_O = Pcirc_feature_N.extractNuc(circ_filename), \
                   Pcirc_feature_J.feature_J(junction_filename), \
                   Pcirc_feature_O.get_orf(orf_f)

df = pd.concat([df_N, df_J, df_O], axis=1)  # Merge the datafram
replace_value = 0.0

print('Ready for seq info ending')
print('Ready for prediction')

model_file = open(model_filename, 'rb')
model = pickle.load(model_file)
headers = list(df.index)
pre = model.predict(df.values[:, :])

for (index, value) in enumerate(list(pre)):
    if value == 1:
        name = headers[index]
        res_file.write(name + '\n' + circ_dict[name] + '\n')
    else:
        pass

res_file.close()
print('Predict ending')
