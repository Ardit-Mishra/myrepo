
# PipelineProject_Ardit_Mishra

# Install BioPython
#pip install biopython

# Install Bowtie2
#sudo apt-get install bowtie2

# Install SPAdes
#sudo apt-get install spades

# Install the necessary libraries
import os 
import Bio
from Bio import SeqIO
from Bio import Entrez

# #1
# Create a new directory for the project
os.system("mkdir PipelineProject_Ardit_Mishra")
os.chdir("PipelineProject_Ardit_Mishra")

# Generate the log file
log_file = open("PipelineProject.log", "w")
log_file.write("Log file generated.")
log_file.close()

# Download the SRA files
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045")

# Convert SRA files to paired-end fastq files

#fastq-dump uncompresses the SRA data
#-I : appends _1 and _2 after read files they are paired-end
#--split-files: splits the paired-end reads into two files
os.system("fastq-dump -I --split-files SRR5660030")
os.system("fastq-dump -I --split-files SRR5660033")
os.system("fastq-dump -I --split-files SRR5660044")
os.system("fastq-dump -I --split-files SRR5660045")

# #2
# Assembling these transciptome reads
#Before assembly, need to map reads to the HCMV genome

# Download the HCMV reference genome from NCBI
Entrez.email = 'arditmishra@gmail.com'  # Put your email address here
handle = Entrez.efetch(db='nucleotide', id='NC_006273.2', rettype='fasta', retmode='text')
fasta = handle.read()
handle.close()

# Save the fasta file to disk
with open('HCMV.fasta', 'w') as f:
    f.write(fasta)

# Before mapping, we need to create "index" to map
# Use Bowtie2 to create an index
os.system('bowtie2-build HCMV.fasta HCMV_index')

# map the transcriptomes to the HCMV genome using Bowtie2
os.system("bowtie2 -x HCMV_index -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S SRR5660030.sam")
os.system("bowtie2 -x HCMV_index -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S SRR5660033.sam")
os.system("bowtie2 -x HCMV_index -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S SRR5660044.sam")
os.system("bowtie2 -x HCMV_index -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S SRR5660045.sam")

#Count the number of reads before and after mapping and save to log file
#Define a function to count the number of reads in a fastq file

def count_reads(fastq_file):
    count = 0
    with open(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
    return count

# Count the number of reads before mapping
log_file = open("PipelineProject.log", "a")
log_file.write("Number of reads before mapping:\n")
log_file.write("SRR5660030: " + str(count_reads("SRR5660030_1.fastq")) + "\n")
log_file.write("SRR5660033: " + str(count_reads("SRR5660033_1.fastq")) + "\n")
log_file.write("SRR5660044: " + str(count_reads("SRR5660044_1.fastq")) + "\n")
log_file.write("SRR5660045: " + str(count_reads("SRR5660045_1.fastq")) + "\n")
log_file.close()

# Count the number of reads after mapping
log_file = open("PipelineProject.log", "a")
log_file.write("Number of reads after mapping:\n")
log_file.write("SRR5660030: " + str(count_reads("SRR5660030.sam")) + "\n")
log_file.write("SRR5660033: " + str(count_reads("SRR5660033.sam")) + "\n")
log_file.write("SRR5660044: " + str(count_reads("SRR5660044.sam")) + "\n")
log_file.write("SRR5660045: " + str(count_reads("SRR5660045.sam")) + "\n")
log_file.close()

# #3
# Run SPAdes to assemble all four transcriptomes together
os.system("spades.py -k 77,99,127 -t 2 --only-assembler -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -o assembly")

#Write the SPAdes command used to the log file
log_file = open("PipelineProject.log", "a")
log_file.write("\nSPAdes command used: spades.py -k 77,99,127 -t 2 --only-assembler -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -o assembly")
log_file.close()

# #4
# Count the number of contigs with length > 1000 and calculate the total length of these contigs
#Specify where is the assembly file
assembly_file = 'assembly/contigs.fasta'

# Counters
contig_count = 0
total_length = 0

# Parse the assembly file in fasta format
records = SeqIO.parse(assembly_file, "fasta")

# Iterating through contigs in the assembly file
for record in records:
    if len(record.seq) > 1000:
        contig_count += 1
        total_length += len(record.seq)

print(f"There are {contig_count} contigs > 1000 bp in the assembly.")
print(f"There are {total_length} bp in the assembly.")

# Write the results to the log file
with open("PipelineProject.log", "a") as log_file:
    log_file.write(f"There are {contig_count} contigs > 1000 bp in the assembly.\n")
    log_file.write(f"There are {total_length} bp in the assembly.\n")

# #5

# Define the taxid
taxid = "txid10359"

# Retrieve the longest contig from the SPAdes assembly
assembly_file = "assembly/contigs.fasta"
contigs = SeqIO.parse(assembly_file, "fasta")
longest_contig = max(contigs, key=lambda x: len(x.seq))
print(f"The longest contig is {longest_contig.id} with length {len(longest_contig.seq)}")


# Create a local database of just sequences from the Betaherpesvirinae subfamily
Entrez.email = "arditmishr@gmail.com" # Put your email address here
handle = Entrez.esearch(db="nucleotide", term="Betaherpesvirinae[Organism]", idtype="acc")
record = Entrez.read(handle)
handle.close()
id_list = record["IdList"]
handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta")
with open("Betaherpesvirinae.fna", "w") as f:
    f.write(handle.read())
handle.close()
os.system("makeblastdb -in Betaherpesvirinae.fna -dbtype nucl")

# Use blastn to query the nr nucleotide database with the longest contig
query_file = open("query.fasta", "w")
query_file.write(">" + longest_contig.id + "\n" + str(longest_contig.seq) + "\n")
query_file.close()

blast_cmd = "blastn -db nr -query query.fasta -outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle' -out blast_output.txt -max_hsps 1 -max_target_seqs 10 -remote -entrez_query 'Betaherpesvirinae[Organism]'"
os.system(blast_cmd)

# Write the top 10 hits to the log file
log_file = open("PipelineProject.log", "a")
log_file.write("Subject accession\tPercent identity\tAlignment length\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject\tBit score\tE-value\tSubject Title\n")
with open("blast_output.txt", "r") as f:
    for line in f:
        log_file.write(line)
log_file.close()

print("blast_output.txt")
