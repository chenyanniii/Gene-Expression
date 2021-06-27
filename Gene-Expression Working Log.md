# Gene Expression Working Log
***
Ideal Procedure
1. Generating a Trinity de novo RNA-Seq assembly
1. Evaluating the quality of the assembly
1. Quantifying transcript expression levels
1. Identifying differentially expressed (DE) transcripts
1. Functionally annotating transcripts using Trinotate and predicting coding regions using TransDecoder
1. Examining functional enrichments for DE transcripts using GOseq
1. Interactively Exploring annotations and expression data via TrinotateWeb
***


## Previous Records (till 6/26/21)
### *De novo* Genome Assembly Using Trinity

#### Trinity
Trinity assembles transcript sequences from Illumina RNA-Seq data. More details could be find on [Trinity Wiki](https://https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#running-trinity-using-singularity). (See DNA reads for genome assembly and their differences.)

Trinity could be run after installation (Code 1), or run through singularity image (Code 2). Certain things need to pay attention:
    * Trinity needs use RAM in the node, quick transfer lots of written small file. Would prefer hight RAM small computer over super computer with low RAM.
    * Trinity is not able to handle too many extra threads or memory.

```
## Code 1
/opt/apps/Software/trinityrnaseq-Trinity-v2.8.5/Trinity \
 --seqType fq \
 --max_memory 50G \
 --left left_1.gz \
 --right right_2.gz \
 --no_normalize_reads
```
```
## Code 2
singularity exec -e Trinity.simg  Trinity \
          --seqType fq \
          --left `pwd`/reads_1.fq.gz  \
          --right `pwd`/reads_2.fq.gz \
          --max_memory 1G --CPU 4 \
          --output `pwd`/trinity_out_dir
```

#### Post Assembly Evaluation
1. Count the mumber of contigs of Trinity.fasta
```
wc -l <file name> 
## results divide four = number of contigs
```
2. use TrinityStats.pl to estimate the quality
```
## installed trinity
/opt/apps/Software/trinityrnaseq-Trinity-v2.8.5/util/TrinityStats.pl Trinity_BOGR.fasta (sphagnum)

## singularity image
singularity exec -e ~/trinityrnaseq.v2.12.0.simg sh -c '$TRINITY_HOME/util/TrinityStats.pl Trinity.fasta' > Trinity_assembly.metrics (nocona)
```
3. Map reads to the assemblies
```
##### Build Bowtie2 Index ##### 
singularity exec -e ~/trinityrnaseq.v2.12.0.simg bowtie2-build ../trinity/trinity_out_dir2/Trinity.fasta Trinity_BOGR (nocona)

##### map reads and calculate alignment statics ##### 
sbatch scripts/map_to_bowtie2_index_BOGR.sh bowtie2/Trinity_BOGR \ fastp_reads/BOGRSW1_R1.fastq.gz,fastp_reads/BOGRSW2_R1.fastq.gz,fastp_reads/BOGRSW3_R1.fastq.gz,fastp_reads/BOGRSW4_R1.fastq.gz,fastp_reads/BOGRW1_R1.fastq.gz,fastp_reads/BOGRW2_R1.fastq.gz,fastp_reads/BOGRW3_R1.fastq.gz,fastp_reads/BOGRW4_R1.fastq.gz \ fastp_reads/BOGRSW1_R2.fastq.gz,fastp_reads/BOGRSW2_R2.fastq.gz,fastp_reads/BOGRSW3_R2.fastq.gz,fastp_reads/BOGRSW4_R2.fastq.gz,fastp_reads/BOGRW1_R2.fastq.gz,fastp_reads/BOGRW2_R2.fastq.gz,fastp_reads/BOGRW3_R2.fastq.gz,fastp_reads/BOGRW4_R2.fastq.gz

##### scripts/map_to_bowtie2_index_BOGR.sh 
#!/bin/bash
#SBATCH -J map_to_bowtie2_BOGR
#SBATCH -p nocona
#SBATCH -o array-%A_%a.out
#SBATCH -e array-%A_%a.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chenyanniii@gmail.com

readonly SINGULARITY_EXEC='singularity exec -e /home/yannchen/trinityrnaseq.v2.12.0.simg'

#1 $1 name of your assembly (without the .fasta suffix)
#2 comma separated list of left read file names
#3 comma separated list of right read file names

${SINGULARITY_EXEC} bowtie2 -p 64 -q --no-unal -k 20 -x $1 -1 $2 -2 $3  2>align_stats_unfixrm_BOGR.txt| ${SINGULARITY_EXEC} samtools view -@64 -Sb -o bowtie2_unfixrm_BOGR.bam
```
4. BUSCO Score
```
##### Runing BUSCO in interactive session
busco -i trinity/trinity_out_dir2/Trinity.fasta -o BUSCO_BOGR -m transcriptome --auto-lineage \
--out_path BUSCO -f


##### Runing BUSCO through job submission
#!/bin/bash
#SBATCH -J BUSCO_SACO
#SBATCH -p nocona
#SBATCH -o array-%A_%a.out
#SBATCH -e array-%A_%a.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=249G
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chenyanniii@gmail.com

# $1 input fasta file (your assembly, e.g. Trinity.fasta)
# $2 output directory (BUSCO will prepend run_)
# $3 output path

busco -i $1 -o $2 -m transcriptome --auto-lineage -f --cpu 64
```

4. Evaluate assembled transcript by comparing the known proteins.

```
# singularity exec ~/trinityrnaseq_latest.sif Trinity --seqType fq --left SACO_R1.gz --right SACO_R2.gz --max_memory 499G --SS_lib_type RF --output trinity_out_dir_SACO
```

#### Quantifying transcript expression levels
##### Identify KAI2 Orthologues
HMMER: remote homologues search based on protein family domain. 
```
## Build KAI2 homologue search protein database in hmm format, using KAI2 genes identified within 1 KP project of 31 species which have whole genome.

hmmbuild hmmbuild KAI2_31genome.hmm aligned_KAI2_31genomespecies.sto
```
Transdecoder identifies candidate coding regions within transcript sequences. 
```
## Build the open reading frame 
TransDecoder.LongOrfs -t trinity/trinity_out_dir2/Trinity.fasta 

## Search KAI2 homolog within the open frame
hmmsearch ~/KAI2/KAI2_31genome.hmm longest_orfs_BOGR.pep > faa_BOGR.out

## TransDecoder.Predict
## TransDecoder.Predict needs to run after TransDecoder.LongOrfs; plus the output_dir need to be the same in order to let ".Predit" use ".LongOrfs"

TransDecoder.Predict -t ../trinity/trinity_out_dir2/Trinity.fasta --output_dir TransDecoder_BOGR.__checkpoints
```














#### Reference Workshops/Tutorials
1. Trinity for de novo assembly (general follow with QC)
https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
2. Using Trinity from Sigularity container
https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-without-a-genome.html#gsc.tab=0
3. Post-assembly analysis (Cornell)
https://biohpc.cornell.edu/doc/RNA-Seq-2019-exercise2-2.html
4. Informatics for RNA-Seq Analysis (Canadian Bioinformatics Workshops)
https://bioinformaticsdotca.github.io/rnaseq_2017_tutorial6
5. Gene-level RNA-Seq Data Analysis (Gitbook)
https://ycl6.gitbook.io/guide-to-rna-seq-analysis/
6. HMMER 
http://www.csb.yale.edu/userguides/seq/hmmer/docs/node5.html
7. Transdecoder
https://github.com/TransDecoder/TransDecoder/wiki

