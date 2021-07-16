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
* Stage 1: Pilot Analysis, Standard Procedure
* Stage 2: Bias & Contaminatioon Control
*** extra tips:
The code is only use as a working record and reference. The files name and location may not be consistant since the analysis were conducted in two different sever at different times.


## Searching Sequence Homologs (7/12)
Use HMMER package to build sequence databases and search for sequence homologs.

### Build hmmdatabase
```
## align sequences before use
clustalo -i KAI2_31genomespecies.faa -o aligned_KAI2_31genomespecies.sto --outfmt=st --threads=10

hmmbuild hmmbuild KAI2_31genome.hmm aligned_KAI2_31genomespecies.sto
```

### Transdecoder identifies candidate coding regions within transcript sequences
```
## TransDecoder.Predict needs to run after TransDecoder.LongOrfs; plus the output_dir need to be the same in order to let ".Predit" use ".LongOrfs"

TransDecoder.LongOrfs -t trinity/trinity_out_dir2/Trinity.fasta 
```

### hmmsearch for the homologs
```
hmmsearch ~/KAI2/KAI2_31genome.hmm longest_orfs_BOGR.pep > faa_BOGR.out
```

## Gene Ontology Enrichment (7/10 - 7/11)
### Extra GO assignment per gene
```bash
## Conda version of Trinotate missing obo/go-basic.obo.gz in conda environment
## find the script online to put in the environment, showed in the error
## cd /home/yannchen/conda/lib/site_perl/5.26.2/x86_64-linux-thread-multi
## mkdir obo
## cd obo
## wget https://github.com/Trinotate/Trinotate/raw/master/PerlLib/obo/go-basic.obo.gz

extract_GO_assignments_from_Trinotate_xls.pl \
                         --Trinotate_xls trinotate_blastxp_pfam.xls \
                         -G --include_ancestral_terms \
                         > go_annotations.txt       
```
### Run GOseq
Prepare the essential documents for GOseq
```bash
## prepare factor_labeling.txt, if there are specific factor wants to test.
## nano factor_labeling.txt
## Inside the file ${gene_name} (tab) ${gene_length}
## Otherwise could use '--genes_single_factor', a list of genes will be extract from input.

## prepare gene.length.txt
## the orginal perl script in the container include the header line ("shebang") as #!/usr/local/bin perl, which is not creat for the signularity environment.
## download the orginal fasta_seq_length.pl, modified the header line ("shebang") to #!$HOME/conda/bin perl

perl fasta_seq_length.pl  /lustre/scratch/yannchen/NovaSeq_2021/BOGR/trinity/trinity_BOGR.fasta > Trinity.fasta.seq_lens


## the orginal python script in the container include the header line ("shebang") as #!/usr/local/bin perl, which is not creat for the signularity environment.
## download the orginal TPM_weighted_gene_length-2.py, modified the header line ("shebang") to #!$HOME/conda/bin python
## Tips: Python version; '--TPM_matrix use transcript expression (isoform)'

python TPM_weighted_gene_length-2.py \
         --gene_trans_map ../trinity/trinity_postR_BOGR.Trinity.fasta.gene_trans_map \
         --trans_lengths Trinity.fasta.seq_lens \
         --TPM_matrix ../GeneExpression/TMM.isoform.TMM.EXPR.matrix > Trinity.gene_lengths.txt
         
## the orginal perl script in the container include the header line ("shebang") as #!/usr/local/bin perl, which is not creat for the signularity environment.
## download the orginal run_GOseq.pl, modified the header line ("shebang") to #!$HOME/conda/bin perl

## Install "goseq" through bioconductor
## if (!requireNamespace("BiocManager", quietly = TRUE))
##    install.packages("BiocManager")
## BiocManager::install("goseq")

## For my purpose, using the DESeq2 results for '--gene_single_factor' and '--background'

## BOGRSW-UP gene group
perl run_GOseq.pl --genes_single_factor \
../GeneExpression/DESeq2.44719_5.dir/TMM.gene.counts.matrix.BOGRSW_vs_BOGRW.DESeq2.DE_results.P1e-3_C2.BOGRSW-UP.subset \
--GO_assignments go_annotations.txt \
--lengths Trinity.gene_lengths.txt \
--background ../GeneExpression/DESeq2.44719_5.dir/TMM.gene.counts.matrix.BOGRSW_vs_BOGRW.DESeq2.DE_results 


## BOGRW-UP gene group
perl run_GOseq.pl --genes_single_factor \
../GeneExpression/DESeq2.44719_5.dir/TMM.gene.counts.matrix.BOGRSW_vs_BOGRW.DESeq2.DE_results.P1e-3_C2.BOGRW-UP.subset \
--GO_assignments go_annotations.txt \
--lengths Trinity.gene_lengths.txt \
--background ../GeneExpression/DESeq2.44719_5.dir/TMM.gene.counts.matrix.BOGRSW_vs_BOGRW.DESeq2.DE_results 

```








## Functional Annotation of Transcripts (7/10/21)
### Add annotations to expression matrix
Generate a map of feature identifier to an annotated feature identifier
```bash
Trinotate_get_feature_name_encoding_attributes.pl \
                  Trinotate_blastxp_pfam.xls > annot_feature_map.txt
```
Integrate functional annotations
```bash
## the orginal perl script in the container include the header line ("shebang") as #!/usr/local/bin perl, which is not creat for the signularity environment.
## download the orginal rename_matrix_feature_identifiers.pl, modified the header line ("shebang") to #!$HOME/conda/bin perl

perl rename_matrix_feature_identifiers.pl ../GeneExpression/TMM.gene.counts.matrix annot_feature_map.txt > TMM_gene.counts.wAnnot.matrix
```




## Gene Annotation (7/9/2021 - 7/10/2021)
### Identification of likely protein-coding region using TransDecoder
```bash
TransDecoder.LongOrfs -t ../trinity/trinity_BOGR.fasta

TransDecoder.Predict -t ../trinity/trinity_BOGR.fasta
```
### Sequence homology searches
blastx against SWISSPROT database to identify likely full-length transcript. (## Using conda install to install trinotate and prepare related databases.)
```bash
## conda install trinotate
## Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
## makeblastdb -in uniprot_sprot.pep -dbtype prot
## gunzip Pfam-A.hmm.gz
## hmmpress Pfam-A.hmm

blastx -db uniprot_sprot.pep \
            -query ../trinity/trinity_BOGR.fasta \
            -num_threads 8 \
            -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
            > swissprot.blastx.outfmt6
```

Look for sequence homologies by just searching our predicted protein sequences rather than using the entire transcript.
```bash
blastp -query trinity_BOGR.fasta.transdecoder.pep \
             -db uniprot_sprot.pep -num_threads 8 \
             -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
              > swissprot.blastp.outfmt6
```
Identify conserved domains that might be indicative or suggestive of function, running a HMMER search against the Pfam database.
```bash
hmmscan --cpu 128 --domtblout trinotatePFAM128.out \
              Pfam-A.hmm \
              trinity_BOGR.fasta.transdecoder.pep
```
### Generate Trinotate Annotation Report
Load the Trinity transcripts and predicted protein sequences
```bash
## If input wrong file into the Trinotate.squlite, use the following comment to rebuild a new one.
## Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

Trinotate Trinotate.sqlite init \
--gene_trans_map ../trinity/trinity_postR_BOGR.Trinity.fasta.gene_trans_map \
--transcript_fasta ../trinity/trinity_BOGR.fasta \
--transdecoder_pep trinity_BOGR.fasta.transdecoder.pep
```
Load various outputs generate in search (blasx, blastp, pfam, signal.out)
```bash
Trinotate Trinotate.sqlite \
           LOAD_swissprot_blastx swissprot.blastx1282.outfmt6

Trinotate Trinotate.sqlite \
           LOAD_swissprot_blastp swissprot.blastp1282.outfmt6

Trinotate Trinotate.sqlite LOAD_pfam trinotatePFAM1282.out
```
### Generate the Trinotate Annotation Report
```bash
Trinotate Trinotate.sqlite report > Trinotate_blastxp_pfam.xls
```


## Differential Expression Analysis (7/8/2021)
### Running Differential Expression Analysis 
Outputs of list of differential analysis, and volcano plot
```bash
## The "--method" and "--samples" could be adjusted based on needs
## demo with samples of 2 vs 3, using methods DESeq2

singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix TMM.gene.counts.matrix \
--method DESeq2 \
--samples fastp_reads/samples_BOGR_replicates_5.txt'
```
### Extracting and clustering differentially expressed transcripts
```bash
cd DESeq2.44719_5.dir

singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../TMM.gene.TMM.EXPR.matrix -P 1e-3 -C 2 \
--samples ../fastp_reads/samples_BOGR_replicates_5.txt'
```
### Automatically Partitioning
```bash
singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
-R  /scratch/yannchen/Transcriptome/NovaSeq_2021/DESeq2.44719_5.dir/diffExpr.P1e-3_C2.matrix.RData \
--Ptree 60'
```





## QC Samples and Biological Replicates
### Compare gene counts between biological samples (7/8/2021)
```bash
singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/Analysis/DifferentialExpression/PtR --matrix TMM.gene.counts.matrix \
                  --samples fastp_reads/samples_BOGR_replicates.txt --log2 --CPM \
                  --min_rowSums 10 \
                  --compare_replicates'
```
### Compare gene counts across samples
```bash
singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/Analysis/DifferentialExpression/PtR --matrix TMM.gene.counts.matrix \
                  --samples fastp_reads/samples_BOGR_replicates.txt --log2 --CPM \
                  --min_rowSums 10 \
                  --sample_cor_matrix'
```
### PCA analysis amount samples
```bash
singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/Analysis/DifferentialExpression/PtR --matrix TMM.gene.counts.matrix \
                  --samples fastp_reads/samples_BOGR_replicates.txt --log2 --CPM \
                  --min_rowSums 10 \
                  --prin_comp 2'
```


## Transcript expression quantitation (7/7/2021)
### Transcript expression quantitation using Salmon

```bash
## the sequence files and the sample.txt need to be in one folder
singularity exec -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/util/align_and_estimate_abundance.pl \
--seqType fq \
--SS_lib_type RF \
--samples_file /scratch/yannchen/Transcriptome/NovaSeq_2021/fastp_reads/samples_BOGR.txt \
--transcripts /scratch/yannchen/Transcriptome/NovaSeq_2021/trinity/trinity_BOGR.fasta \
--est_method salmon \
--trinity_mode \
--thread_count 32 \
--prep_reference'
```

Creat a list of "quant.sf" files
```bash
find BOGR* -name 'quant.sf' | tee quant_files.list
```

Using the "quant_files.list" to generation the count and expression matrices (isoforms and genes)

Tips: Delete cetain genes in trinity_BOGR.fasta, saved as fixedtry2.fasta, since salmon process delete some of the file.
>TRINITY_DN0_c0_g2_i1
>TRINITY_DN0_c0_g1_i2
>TRINITY_DN0_c0_g5_i1
```bash
singularity exec -H /scratch/yannchen -e ~/trinityrnaseq.v2.12.0.simg sh -c \
'$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
--out_prefix Trinity \
--name_sample_by_basedir \
--quant_files quant_files.list \
--gene_trans_map /scratch/yannchen/Transcriptome/NovaSeq_2021/trinity/fixedtry2.fasta'
```








## Assembly Evaluation (NovoSeq)
### How many transcriptome reads in the assembly?
```bash
grep '>' trinity_BOGR.fasta | wc -l
```
### Assembly stats
Using the scripts inside the trinity to collect some basic statistical information about the assembly.
```bash
singularity exec ~/trinityrnaseq.v2.12.0.simg sh -c $TRINITY_HOME/util/TrinityStats.pl trinity_BOGR.fasta' > BOGR_assembly.metrics
```
### Quantify how well the reads could support for the assembly
Using bowtie2-build to build the bowtie2 reference of the assembly, then map the reads back to the reference.If the overall map rate is 70% or higher, the assembly quality is decent.
```bash
singularity exec -e ~/trinityrnaseq.v2.12.0.simg bowtie2-build trinity/trinity_BOGR.fasta trinity/trinity_BOGR.fasta

singularity exec -e ~/trinityrnaseq.v2.12.0.simg bowtie2 --local -p 64 -q --no-unal -k 20 -x trinity/trinity_BOGR.fasta -1 fastp_reads/BOGRSW1_R1.fastq.gz,fastp_reads/BOGRSW2_R1.fastq.gz,fastp_reads/BOGRSW3_R1.fastq.gz,fastp_reads/BOGRSW4_R1.fastq.gz,fastp_reads/BOGRW1_R1.fastq.gz,fastp_reads/BOGRW2_R1.fastq.gz,fastp_reads/BOGRW3_R1.fastq.gz,fastp_reads/BOGRW4_R1.fastq.gz,fastp_reads/BOGR_R1.fastq.gz -2 fastp_reads/BOGRSW1_R2.fastq.gz,fastp_reads/BOGRSW2_R2.fastq.gz,fastp_reads/BOGRSW3_R2.fastq.gz,fastp_reads/BOGRSW4_R2.fastq.gz,fastp_reads/BOGRW1_R2.fastq.gz,fastp_reads/BOGRW2_R2.fastq.gz,fastp_reads/BOGRW3_R2.fastq.gz,fastp_reads/BOGRW4_R2.fastq.gz,fastp_reads/BOGR_R2.fastq.gz  2>align_stats_postR_BOGR.txt |samtools view -@64 -b -o bowtie2_BOGR.bam
```
### BUSCO
```bash
## busco is in a separate conda environment on sphagnum, need to use conda activate and deativate to triger or kill the software. 
busco -i ../trinity/trinity_BOGR.fasta -o BOGR_BUSCO -m transcriptome --auto-lineage -f --cpu 64
```





## Previous Records (till 6/26/21; MiSeq)
### *De novo* Genome Assembly Using Trinity



#### Trinity
Trinity assembles transcript sequences from Illumina RNA-Seq data. More details could be find on [Trinity Wiki](https://https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#running-trinity-using-singularity). (See DNA reads for genome assembly and their differences.)

Trinity could be run after installation (Code 1), or run through singularity image (Code 2). Certain things need to pay attention:
    * Trinity needs use RAM in the node, quick transfer lots of written small file. Would prefer hight RAM small computer over super computer with low RAM.
    * Trinity is not able to handle too many extra threads or memory.

```bash
## Code 1
/opt/apps/Software/trinityrnaseq-Trinity-v2.8.5/Trinity \
 --seqType fq \
 --max_memory 50G \
 --left left_1.gz \
 --right right_2.gz \
 --no_normalize_reads
```
```bash
## Code 2
singularity exec -e Trinity.simg  Trinity \
          --seqType fq \
          --left `pwd`/reads_1.fq.gz  \
          --right `pwd`/reads_2.fq.gz \
          --max_memory 1G --CPU 4 \
          --output `pwd`/trinity_out_dir
```
Before Running Trinity, there are certain steps should be precessed to ensure the input quality and avoid overwhelming data set. For more detail could refer to [Harvard FAS Informatics](https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html).
* Quality Control: Fastqc
* Trim adapter and low quality bases (PolyA, PolyG): fastp
* Remove erroneous k-mers: rCorrector 
* Remove unfixable reads: FilterUncorrectabledPEfastq.py
* Remove unwanted (rRNA reads): SILVA + bowtie2 (Optional)
* Quality Control (second run): Fastqc
* Remove over-represented sequences (Optional)

#### Post Assembly Evaluation
1. Count the mumber of contigs of Trinity.fasta
```bash
wc -l <file name> 
## results divide four = number of contigs
```
2. use TrinityStats.pl to estimate the quality
```bash
## installed trinity
/opt/apps/Software/trinityrnaseq-Trinity-v2.8.5/util/TrinityStats.pl Trinity_BOGR.fasta (sphagnum)

## singularity image
singularity exec -e ~/trinityrnaseq.v2.12.0.simg sh -c '$TRINITY_HOME/util/TrinityStats.pl trinity_BOGR.fasta' > Trinity_assembly.metrics (nocona)
```
3. Map reads to the assemblies
```bash
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
```bash
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

```bash
# singularity exec ~/trinityrnaseq_latest.sif Trinity --seqType fq --left SACO_R1.gz --right SACO_R2.gz --max_memory 499G --SS_lib_type RF --output trinity_out_dir_SACO
```

#### Quantifying transcript expression levels
##### Identify KAI2 Orthologues
HMMER: remote homologues search based on protein family domain. 
```bash
## Build KAI2 homologue search protein database in hmm format, using KAI2 genes identified within 1 KP project of 31 species which have whole genome.

hmmbuild hmmbuild KAI2_31genome.hmm aligned_KAI2_31genomespecies.sto
```
Transdecoder identifies candidate coding regions within transcript sequences. 
```bash
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

