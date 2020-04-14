import argparse
import os

parser = argparse.ArgumentParser(description='Sequence Aligment \n '
                                             'usage: \n '
                                             '[options] <reads1[,reads2,...]>')

parser.add_argument("-g",
                    "--genome",
                    help='reference genome files are Fasta')
parser.add_argument("-G",
                    "--gtf",
                    help='reference gff or gtf files')
parser.add_argument("-p",
                    "--num_threads",
                    default=4)
parser.add_argument("-o",
                    "--output_dir")
parser.add_argument("RNA_seq",
                    nargs='*',
                    help='data of RNA_seq, the input format is <reads1[,reads2,...]>')

args = parser.parse_args()

inputs = args.RNA_seq
threads = args.num_threads
output_dir = args.output_dir
gtf_filename = args.gtf
genome_filename = args.genome

genome_path = os.path.realpath(genome_filename)
genome_index = genome_path[:genome_path.rfind('/')] + '/genome_index'
os.system("bowtie2-build %s %s" %
          (genome_filename, genome_index))

num_fq = len(inputs)
if num_fq == 1:
    os.system("tophat2 -p %d -o %s -G %s %s" %
              (threads, output_dir, gtf_filename, inputs[0]))

elif num_fq == 2:
    os.system("tophat2 -p %d -o %s -G %s %s %s" %
              (threads, output_dir, gtf_filename, inputs[0], inputs[1]))

os.system("samtools view %s/unmapped.bam | "
          "awk \'{OFS=\"\t\"; print \">\"$1\"\n\"$10}\' - > "
          "%s/unmapped.fasta" %
          (output_dir, output_dir))