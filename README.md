# HaemoSeq
HaemoSeq is a bioinformatic pipeline for whole genome sequencing based analysis of human-related Haemophilus isolates including H. influenzae (typeable and NTHi), H. haemolyticus (subsp. intermedius), H. parahaemolyticus, H. paraphrohaemolyticus, H. parainfluenzae, H. sputorum, H. pittmaniae and H. ducreyi.
It runs on linux and accepts FastQ (illumina) and FastA files as input.
It includes:

1.	Quality control of raw sequence reads
[Input: FastQ files; Required tools: FastQC and multiQC]
2.	Contamination detection 
[Input: FastQ files or FastA files; Required tools: kraken2]
3.	Preprocessing of raw sequence reads (e.g. adapter removal) 
[Required tools: fastp]
4.	Haemophilus (sub)species and serotype prediction using a custom marker database 
[Input: FastQ files; Required tools: srst2 and GNU parallel; Required database: Class_Haemophilus_and_Serotyping-v2.0.fasta]
5.	Multi-locus sequence typing (MLST)
[Input: FastQ files; Required tools: srst2; Required database: pubMLST]
6.	Assembly of sequence reads
[Input: FastQ files; Required tools: Shovill; Output: FastA files]
7.	Fast phylogenetic analysis using Mashtree
[Input: FastQ files or FastA files; Required tools: Mashtree]
8.	Detection of known plasmids
[Input: FastQ files; Required tools: srst2 and seqkit; Required database: PLSDB or custom]
9.	De novo prediction of plasmid contigs
[Input: FastQ files; Required tools: platon and/or plasmidspades]
10.	Resistance and virulence gene prediction
[Input: FastA files; Required tools: AMRfinder+]


# Getting Started
For complete installation instructions, description and usage examples please send a mail to mdiricks@fz-borstel.de.

# Citation
Diricks et al., 2022: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01017-x

