# happypenguin
Find crossovers in DNA long reads from F2 populations.

*happypenguin* uses variation between the genotypes to determine their genotype. Reads that exhibit a switch in single nucleotide polymorphisms (SNPs) between the two genomes indicate a crossover point during F1 meiosis.

## Usage
> ./happypenguin alignment_1.sam alignment_2.sam syri.out name

### alignment_1.sam
All reads aligned to reference genome 1.

For example, to align using minimap2 (https://github.com/lh3/minimap2):
> minimap2 --eqx -a genomeA.fasta reads.fasta

### alignment_2.sam
All reads aligned to reference genome 2.

### syri.out
SyRI (https://schneebergerlab.github.io/syri/fileformat.html) is used to call variation between the two parental genomes. 
First, align the parental genomes together:
> minimap2 --eqx -ax asm10 genomeA.fa genomeB.fa > aln.sam
> syri -c aln.sam -r genomeA.fa -q genomeB.fa -k -F S --cigar
Further filtering of SNPs for homozygosity is recommended.

## Output
Each read is sorted into one of the following categories:
- **Genotype 1**
- **Genotype 2**
- **Low Quality Crossover**
- **Medium Quality Crossover**
- **High Quality Crossover**
- **Chomp** reads can be chomped (inconclusive) because:
1. Bad Alignment: Read mapped to different chromosomes
2. Bad Genotype: Consensus SNP Matrix switches gt more than once
3. No SNPs: Zero SNPs in Consensus SNP Matrix
4. Multimap Rec: Read appears to a crossover but has multiple mappings

## Installation
**happypenguin** is compiled using:
> g++ happypenguin.cc -std=c++11 -O3 -o happypenguin

## Algorithm
**happypenguin** compares the SNPs of both alignments. In the case that both alignments contain the SNP or neither contain the SNP, the SNP is ignored in categorization. This allows *happypenguin* to be relatively robust to sequencing errors, heterozygosity in variation calling, and inconsistency in alignments.

## Troubleshooting
**happypenguin** only uses variations that occur in syntenic regions. If synteny is low, many reads will be lost to **ChompNoSNPs**. In this case, different genomic alignment methods and/or parameters should be tested.
