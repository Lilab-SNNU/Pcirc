import argparse
import os

parser = argparse.ArgumentParser(description='Unmapped reads alignment')

parser.add_argument("-g",
                    "--genome",
                    help='reference genome files are Fasta')
parser.add_argument("-q",
                    "--query",
                    help='reference query files are Fasta')
parser.add_argument("-o",
                    "--output_file",
                    help='Output_file saved reasult')
parser.add_argument("-e",
                    "--evalue",
                    default=1e-5,
                    help='Expectation value (E) threshold for saving hits')
parser.add_argument("-p",
                    "--num_threads",
                    default=4)


args = parser.parse_args()

evalue = args.evalue
num_threads = args.num_threads
output_file = args.output_file
query_filename = args.query
genome_filename = args.genome

genome_path = os.path.realpath(genome_filename)
genome_database = genome_path[:genome_path.rfind('/')] + '/genome_database'

os.system("makeblastdb -in %s -dbtype nucl -out %s" %
          (genome_filename, genome_database))

os.system("blastn -db %s -query %s -out %s -outfmt 6 -evalue %f -num_thread %d" %
          (genome_database, query_filename, output_file, evalue, num_thread))
