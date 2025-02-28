PipelineProject_Sara_Ali
For the first problem, both SRA files were dowloaded using wget. This was found from the links that were provided in the homework. The links I chose in the SRA Normalized row. 
I used the codes: wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030, and wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

I made a python script that produced a text tile with output named pipelineproject.log in a directory named PipelineProject_Sara_Ali. This was done  in terminal using the commands: mkdir -p PipelineProject_Sara_Ali, then cd PipelineProject_Sara_Ali, touch PipelineProject.log 

Then I searched the accession number for HCMV in genome assembly'. Then I clicked ‘datasets’ and found the code. This dataset was copied into the command and same method was used for problem 5.
datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report
datasets download virus genome taxon Betaherpesvirinae --include genome

Make sure to include the dependencies in the code:
import Bio
import datasets
import bowtie2
import samtools
import spades
import blastn
import os
import argparse
from Bio import SeqIO

This was the documentation that includes what dependencies are needed to be installed to run your code.
