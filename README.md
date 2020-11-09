# Pcirc

## Introduction

PCirc is a pipeline to predict plant circular RNA (CircRNA) which based on ***Python3***.  It can identify circRNA from  a given ***RNA-seq*** data by high-throughput. The pipeline can alignment RNA_seq and extract features ,then based on Random Forest classification, the final circRNA will be predicted. The user only needs to ***run the python scripts*** step by step to get the predicted circular RNA.

![PCirc](https://github.com/ChasingWind/Pcirc_v1.0/blob/master/Figure/Pcirc.jpg)

---

## Requirement
### Data:

- Genome sequence (fasta format)
- Genome annotation file (gtf/gff format)
- RNA_seq data (fastq format)
### Software:

#### Alignment

- NCBI-blast(v.2.9.0+):(ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- bowtie2 (v.2.2.6+): (https://github.com/BenLangmead/bowtie2)
- tophat2 (v.2.1.1+): (http://ccb.jhu.edu/software/tophat)
- samtools (v.0.1.19+): (https://github.com/samtools/samtools)
- Python3 (v.3.6.5+): (https://www.python.org/)

#### Predicting

- UGENE(1.30.0+):(http://ugene.net/)

### python3 package:

- Biopython (v.1.72+): (https://pypi.org/project/biopython/)
- Pandas (v.0.23.3+): (https://pypi.org/project/pandas/)
- Pyyaml (v.5.3.1): (https://pypi.org/project/pyyaml/)
- Scikit-learn (v.0.21.2): (https://scikit-learn.org/stable/)

### Supported operating systems：
![](https://s2.ax1x.com/2019/08/30/mX31PK.th.jpg)

---
## Download
  Open the terminal and input:
  ```bash
  git clone https://github.com/Lilab-SNNU/Pcirc.git
  ```
---
## Usage

### - You can run Pcirc step by step use command line

  1. Alignment

     In this step need given file:

     - Genome sequence (fasta format)
     - Genome annotation file (gtf/gff format)
     - RNA_seq data (fastq format)

​	Then, you will get a file named ***unmapped.fasta*** in ***output_dir***, it content the information of unmapped reads and will be used to the ***step blast*** 

  ```bash
   python3 Pcirc_aligment.py -g genome.fa -G genome.gtf -o output_dir <reads_1[,reads_2]>
  ```

  2. Blast

     Blast alignment of unmapped reads sequence and genomic sequences, in this step need given file:

     - Genome sequence (fasta format)
     - Unmapped reads from ***step Alignment*** (fasta format)

     Then, you will get a file named unmapped.blast in the dir which unmapped.fasta in, it content the blast information  and will be used to the ***step Circular RNA predicting***

  ```bash
   python3 Pcirc_unmapped_reads_blast.py -g genome.fa -q unmapped.fasta -o output_dir
  ```
  3. Circular RNA predicting

     Final this step, you will get the circular RNA information saved in a file name circ. 

  ```bash
   python3 Pcirc_predict.py -i unmapped.blast -g genome.fa
  ```

### - You can also run Pcirc in graphical user interface(GUI) 
   If you use the GUI, please ensure the R, Rstudio and Rshiny in your PC.

---
## Run example

You can download the required sra file from [NCBI-SRA](https://www.ncbi.nlm.nih.gov/sra/SRR3495992) or fasta file from [EMBL-EBI](https://www.ebi.ac.uk/), and run ***Pcirc*** or use the *test.blast* ***run step Circular RNA predicting***. 

---
## Contact us

If you encounter any problems while using Pcirc, please send an email (glli@snnu.edu.cn) or submit the issues on GitHub (https://github.com/Lilab-SNNU/Pcirc/issues) and we will resolve it as soon as possible.
