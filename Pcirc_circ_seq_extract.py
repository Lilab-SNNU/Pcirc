import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def reverse(seq):
    my_seq = Seq(seq, IUPAC.unambiguous_dna)
    re_seq = str(my_seq.reverse_complement())
    return re_seq

def get_seq(info_file, genome_file):
    print('start sequence finding')
    df = pd.read_csv(info_file, header=None)
    df = df.drop_duplicates(subset=[1, 2, 3])
    genome = open(genome_file, 'r')
    genome_list = genome.read().split('>')[1:]
    genome_dict = {}
    up_filename = info_file + '_circ_up'
    down_filename = info_file + '_circ_down'
    junc_filename = info_file + '_junc_seq'
    circ_filename = info_file + '_aligned_seq'
    file_up, file_down, file_circ, file_junc = open(up_filename, 'w'), \
                                               open(down_filename, 'w'), \
                                               open(circ_filename, 'w'), \
                                               open(junc_filename, 'w')

    for genome_info in genome_list:
        genome_chr = '>' + genome_info.split('\n')[0].split(' ')[0]
        genome_seq = ''.join(genome_info.split('\n')[1:-1])
        genome_dict[genome_chr] = genome_seq

    for index, row in df.iterrows():
        circ_name, chr_name, strand, start_site, end_site = row[0], '>' + row[1], row[2], int(row[3]), int(row[4])

        if start_site >= 50 and len(genome_dict[chr_name]) - end_site >= 50 and strand == '+':
            circ_seq = genome_dict[chr_name][start_site - 1:end_site]
            circ_up_seq = genome_dict[chr_name][start_site - 50:start_site + 50]
            circ_down_seq = genome_dict[chr_name][end_site - 50:end_site + 50]

        elif start_site < 50 and len(genome_dict[chr_name]) - end_site >= 50 and strand == '+':
            circ_seq = genome_dict[chr_name][start_site - 1:end_site]
            circ_up_seq = genome_dict[chr_name][0:100]
            circ_down_seq = genome_dict[chr_name][end_site - 50:end_site + 50]

        elif start_site >= 50 and len(genome_dict[chr_name]) - end_site < 50 and strand == '+':
            circ_seq = genome_dict[chr_name][start_site - 1:end_site]
            circ_up_seq = genome_dict[chr_name][start_site - 50:start_site + 50]
            circ_down_seq = genome_dict[chr_name][-100:]

        elif start_site >= 50 and len(genome_dict[chr_name]) - end_site >= 50 and strand == '-':
            circ_seq = genome_dict[chr_name][-end_site:-(start_site - 1)]
            circ_up_seq = genome_dict[chr_name][-start_site - 50:-start_site + 50]
            circ_down_seq = genome_dict[chr_name][-end_site - 50:-end_site + 50]
            circ_seq, circ_up_seq, circ_down_seq = reverse(circ_seq), reverse(circ_up_seq), reverse(circ_down_seq)

        elif start_site < 50 and len(genome_dict[chr_name]) - end_site >= 50 and strand == '-':
            circ_seq = genome_dict[chr_name][-end_site:-(start_site - 1)]
            circ_up_seq = genome_dict[chr_name][-100:]
            circ_down_seq = genome_dict[chr_name][-end_site - 50:-end_site + 50]
            circ_seq, circ_up_seq, circ_down_seq = reverse(circ_seq), reverse(circ_up_seq), reverse(circ_down_seq)

        elif start_site >= 50 and len(genome_dict[chr_name]) - end_site < 50 and strand == '-':
            circ_seq = genome_dict[chr_name][-end_site:-(start_site - 1)]
            circ_up_seq = genome_dict[chr_name][-start_site - 50:-start_site + 50]
            circ_down_seq = genome_dict[chr_name][0:100]
            circ_seq, circ_up_seq, circ_down_seq = reverse(circ_seq), reverse(circ_up_seq), reverse(circ_down_seq)


        if len(circ_seq) < 10000:
            file_circ.write('>' + circ_name + '_%s' % chr_name[1:] + '_%s' % strand + '_%s|%s' % (start_site, end_site) + '\n' + circ_seq + '\n')
            file_up.write('>' + circ_name + '_%s' % chr_name[1:] + '_%s' % strand + '_%s|%s' % (start_site, end_site) + '\n' + circ_up_seq + '\n')
            file_down.write('>' + circ_name + '_%s' % chr_name[1:] + '_%s' % strand + '_%s|%s' % (start_site, end_site) + '\n' + circ_down_seq + '\n')
            file_junc.write('>' + circ_name + '_%s' % chr_name[1:] + '_%s' % strand + '_%s|%s' % (start_site, end_site) + '\n' + circ_up_seq + circ_down_seq)

    print('sequence finding　ending')

    return junc_filename, circ_filename
