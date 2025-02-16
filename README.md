# happypenguin
A tool to find crossovers from Nanopore DNA reads from F2 populations.
happypenguin uses variation between the genomes to determine their genotype. Reads that exhibit a switch in single nucleotide polymorphisms (SNPs) between the two genomes indicate a crossover point during F1 meiosis.

## Usage
> ./happypenguin alignment_A.sam alignment_B.sam syri.out [output_name]

### alignment_A.sam
All reads aligned to reference genome A, with match/mismatch information in the CIGAR string.

For example, to align using the minimap2 aligner (https://github.com/lh3/minimap2):
> minimap2 --eqx -a genomeA.fasta reads.fasta

### alignment_A.sam
All reads aligned to reference genome B, like above.

### syri.out
SYRI (https://schneebergerlab.github.io/syri/fileformat.html) is used to call variation between the two parental genomes. 
First, align the parental genomes together:
> minimap2 --eqx -ax asm10 genomeA.fa genomeB.fa > aln.sam 

Now, run syri:
> syri -c aln.sam -r genomeA.fa -q genomeB.fa -k -F S --cigar

## Output
Each read is sorted into one of the following categories:
- **GenotypeA**
- **GenotypeB**
- **LowQualityCrossover**
- **MediumQualityCrossover**
- **HighQualityCrossover**
- **Chomp** reads can be chomped (inconclusive) because:
1. **BadGt** Aligned to different chromosomes
2. **BadGT** Consensus SNP Matrix switches gt more than once
3. **NoSNPs** Less than two SNPs in Consensus SNP Matrix
4. **MultimapRec** Read appears to a crossover but has multiple mappings

## Installation
**happypenguin** is compiled using:
> g++ happypenguin.cc -std=c++11 -O3 -o happypenguin

## Algorithm
**happypenguin** compares the SNPs of both alignments. In the case that both alignments contain the SNP or neither contain the SNP, the SNP is ignored in categorization. This allows **happypenguin** to be robust to sequencing errors, heterozygosity in variation calling, and inconsistency in alignments.

## Troubleshooting
**happypenguin** only uses variations that occur in syntenic regions. If synteny is low, many reads will be lost to **ChompNoSNPs**. In this case, different genomic alignment methods and/or parameters should be tested.
