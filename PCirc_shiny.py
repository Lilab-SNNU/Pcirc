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

'''
PCirc is a pipeline to predict plant circular RNA (CircRNA) which based on Python3. 
It can identify circRNA from a given RNA-seq data by high-throughput. 
You can upload files in this GUI.
You can also get more information in github.
(https://github.com/Lilab-SNNU/Pcirc) 
'''

path = os.path.realpath(__file__)  # Get the install path of Pcirc
abs_dir = path[:path.rfind('/')]

parser = argparse.ArgumentParser(description='Predict plant circRNAs')

parser.add_argument("-g",
                    "--genome",
                    help='reference genome files are Fasta')
parser.add_argument("-G",
                    "--gtf",
                    help='reference gtf files')
parser.add_argument("-m",
                    default=abs_dir + "/model.pkl",
                    help='path to model, if you do not change the path of model, please do not set it')
parser.add_argument("-p",
                    "--num_threads",
                    default=4)
parser.add_argument("-o",
                    "--output_dir",
                    help='the output folder which save the results ')
parser.add_argument("RNA_seq",
                    nargs='*',
                    help='data of RNA_seq, the input format is <reads1[,reads2,...]>')

args = parser.parse_args()
inputs = args.RNA_seq
threads = int(args.num_threads)
output_dir = args.output_dir
gtf_filename = args.gtf
genome_filename = args.genome
model_filename = args.m

genome_path = os.path.realpath(genome_filename)
genome_index = genome_path[:genome_path.rfind('/')] + '/genome_index'

print('Alignment reads begining')
print('='*20 + 'Bowtie begin' + '='*20)
os.system("bowtie2-build %s %s" %
          (genome_filename, genome_index))
print('='*20 + 'Bowtie end' + '='*20)

print('='*20 + 'Tophat begin' + '='*20)
num_fq = len(inputs)
if num_fq == 1:
    os.system("tophat2 -p %d -o %s -G %s %s %s" %
              (threads, output_dir, gtf_filename, genome_index, inputs[0]))

elif num_fq == 2:
    os.system("tophat2 -p %d -o %s -G %s %s %s %s" %
              (threads, output_dir, gtf_filename, genome_index, inputs[0], inputs[1]))
print('='*20 + 'Tophat end' + '='*20)

os.system(r'''samtools view %s/unmapped.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' - > %s/unmapped.fasta''' %
          (output_dir, output_dir))

genome_database = genome_path[:genome_path.rfind('/')] + '/genome_database'
os.system("makeblastdb -in %s -dbtype nucl -out %s" %
          (genome_filename, genome_database))

os.system("blastn -db %s -query %s/unmapped.fasta -out %s/unmapped.blast -outfmt 6 -evalue 1e-5 -num_thread %d" %
          (genome_database, output_dir, output_dir, threads))

info_file = '%s/unmapped.blast' % output_dir


# Pcirc running
alignment_res = '/'.join(info_file.split('/')[:-1]) + '/alignment_res'
df_alis = Pcirc_deal_blast_result.exactReads(info_file)
df_list = Pcirc_deal_blast_result.cutDatafram(df_alis)
pool = multiprocessing.Pool(len(df_list))

for dataframe in df_list:
    pool.apply_async(Pcirc_deal_blast_result.getCircInformation,
                     (dataframe, alignment_res, ))
pool.close()
pool.join()
print('Alignment reads ending')

junction_filename, circ_filename = Pcirc_circ_seq_extract.get_seq(alignment_res, genome_filename)
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
