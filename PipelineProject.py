import Bio
import datasets
import bowtie2
import samtools
import spades
import blastn
import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Wrapper script for HCMV pipeline.")
parser.add_argument("--first_name", required=True)
parser.add_argument("--last_name", required=True)
parser.add_argument("--sra1_r1", required=True)
parser.add_argument("--sra1_r2", required=True)
parser.add_argument("--sra2_r1", required=True)
parser.add_argument("--sra2_r2", required=True)
args = parser.parse_args()

project_dir = "/home/2025/sali32/PipelineProject_Sara_Ali"

#Problem 2: when testing sample data, in this code you would change SRR5660030_1.fastq to sampleSRR5660030_1.fastq. This is where I stored my sample data
os.system("datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report") #downloaded the dataset of HCMV
os.system("unzip ncbi_dataset.zip") #unzipped the dataset to aquire more files

os.system("bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV") #built bowtie2 index for mapping

#running the Bowtie2 mapping below and change the fastq files to sam files to store the read alignments
os.system("bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S mapped_2dpi.sam") #changed to sample here twice
os.system("bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S mapped_6dpi.sam") #changed here as well

#this counts the number of mapped reads for both files
os.system("samtools view -c -F 4 mapped_2dpi.sam")
os.system("samtools view -c -F 4 mapped_6dpi.sam")

#below, both files were change to samples
#logs the number of read pairs before and after filtering into PipelineProject.log
os.system("echo \"Donor 1 (2dpi) had $(expr $(wc -l < SRR5660030_1.fastq) / 4) read pairs before Bowtie2 filtering and $(samtools view -c -F 4 mapped_2dpi.sam) read pairs after.\" >> PipelineProject.log")
os.system("echo \"Donor 1 (6dpi) had $(expr $(wc -l < SRR5660033_1.fastq) / 4) read pairs before Bowtie2 filtering and $(samtools view -c -F 4 mapped_6dpi.sam) read pairs after.\" >> PipelineProject.log")



#Problem 3:

os.system("samtools view -b -F 4 mapped_2dpi.sam > mapped_2dpi.bam") #converts SAM file to BAM format for processing
os.system("samtools fastq -1 mapped_2dpi.1.fastq -2 mapped_2dpi.2.fastq -s mapped_2dpi.unpaired.fastq mapped_2dpi.bam") #extracts paired and unpaired reads from BAM file
os.system("samtools view -b -F 4 mapped_6dpi.sam > mapped_6dpi.bam") #same conversion for 6dpi sample
os.system("samtools fastq -1 mapped_6dpi.1.fastq -2 mapped_6dpi.2.fastq -s mapped_6dpi.unpaired.fastq mapped_6dpi.bam") #extracts reads for 6dpi sample
#assembles the contigs using SPAdes assembler
os.system("spades.py -k 99 -t 1 --isolate -1 mapped_2dpi.1.fastq -2 mapped_2dpi.2.fastq -1 mapped_6dpi.1.fastq -2 mapped_6dpi.2.fastq -o HCMV_combined_assembly/")
os.system("echo \"Size of contigs.fasta: $(du -h HCMV_combined_assembly/contigs.fasta | cut -f1)\" >> PipelineProject.log") #logs the size of the assembled contigs
print("Finished.")




#Problem 4:
import os
contigs_file = "HCMV_combined_assembly/contigs.fasta"
num_contigs = 0
total_bp = 0
with open(contigs_file, "r") as file: #opens the contigs file to count large contigs >1000 bp in the assembly
    for line in file:
        if line.startswith(">"): #identifies headers in FASTA format
            length = int(line.split("_")[3]) #extracts length from contig header 
            if length > 1000: #checks if contig length is greater than 1000 bp
                num_contigs += 1
                total_bp += length #adds length of large contig to total base pairs

log_file = "PipelineProject.log"
with open(log_file, "a") as log:
    log.write(f"There are {num_contigs} contigs > 1000 bp in the assembly.\n")
    log.write(f"There are {total_bp} bp in the assembly.\n")
print("Finished.")
print(f"There are {num_contigs} contigs > 1000 bp in the assembly.")
print(f"There are {total_bp} bp in the assembly.")




#Problem 5:

#performs BLAST search of longest contig against Betaherpesvirinae database
os.system("datasets download virus genome taxon Betaherpesvirinae --include genome")
os.system("unzip ncbi_dataset.zip") #this uses the same method as 2 to obtain the files
os.system("mv ncbi_dataset/data/genomic.fna Betaherpesvirinae.fna")  #moves and renames the genome file so that its easier to use
os.system("makeblastdb -in Betaherpesvirinae.fna -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl") #creates BLAST nucleotide database

contigs_file = "HCMV_combined_assembly/contigs.fasta"
longest_contig = None
max_length = 0

for record in SeqIO.parse(contigs_file, "fasta"): #goes through assembled contigs to find the longest one
    if len(record.seq) > max_length: #compares lengths of contigs
        max_length = len(record.seq)
        longest_contig = record #stores longest contig

longest_contig_file = "problem5_contig_file.fasta" 
with open(longest_contig_file, "w") as output:
    SeqIO.write(longest_contig, output, "fasta") #writes longest contig to a new FASTA file

print(f"Longest contig saved to {longest_contig_file} with length {max_length} bp.")

#this will perform blast search of longest contig against Betaherpesvirinae database
os.system(
    'blastn -query problem5_contig_file.fasta -db Betaherpesvirinae '
    '-out myresults.tsv -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" '
    '-max_target_seqs 10 -max_hsps 1'
)
blast_file = "myresults.tsv"
log_file = "PipelineProject.log"

#appends the results into log
with open(blast_file, "r") as blast:
    blast_results = blast.readlines()
with open(log_file, "a") as log:
    log.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
    log.writelines(blast_results)

print("BLAST results successfully appended to PipelineProject.log.")
