# Reference Information

This workshop is developed based on the [RNA-seq Bioinformatics Course](https://rnabio.org/) by the [Griffith Lab](https://griffithlab.org/) at Washington University, with extensive help from [ChatGPT](https://chatgpt.com/).

# What is RNA-Seq Analysis?

RNA-Seq (RNA sequencing) is a powerful technique for studying the transcriptome, the complete set of RNA scripts produced by a genome at a given time. It allows researchers to identify and quantify RNA molecules, providing insights into gene expression, alternative splicing, non-coding RNAs, and other molecular phenomena. 

__Attention!! - RNA-seq is not working on RNA directly but cDNA, here's why:__

* __RNA isolation:__ The first step in RNA-seq is extracting RNA from the sample, which consists of various types of RNA, such as mRNA, rRNA, and non-coding RNAs. 
* __cDNA synthesis:__ Since sequencing technologies generally work with DNA, the extracted RNA is reverse-transcribed into cDNA using the enzyme reverse transcriptase. This is because the sequencing platforms (like Illumina, for example) read DNA sequences, not RNA. 
* __Sequencing:__ The cDNA, which mirrors the RNA sequence (except for the replacement of uracil with thymine), is then fragmented, and libraries are prepared for sequencing. 

# Before the analyses

Should write about how the experimental design, library prepare, and sequencing technology used can influnce our downstream analyses, such as the use of different pipelines and software, quality control etc.... 

There are many RNA-seq library construction strategies. 

# Types of RNA-seq Analysis  

RNA-seq is a versatile technology with diverse applications, and the specifics of experimental design, library preparation, sequencing technology, and downstream analysis vary significantly based on the research goal.

__Summary of applications and which one we are focusing on?__

* Differential gene expression (DGE)

## Differential Gene Expression (DEG) Analysis

* __Goal__: Compare gene expression levels between conditions (e.g., treated vs. untreated, diseased vs. healthy).
* __Design__: Biological replicates are critical. 
* __Library Prep__: Typically uses poly(A) enrichment or rRNA depletion for mRNA-focused studies. 
* __Sequencing__: Moderate depth (20-50M reads per sample) is sufficient. 
* __Downstream Analysis__: 
    * Alignment tools: STAR, HISAT2.
    * Differential gene expression tools: DESeq2, edgeR.

## Transcript Isoform Analysis

* __Goal__: Detect alternative splicing, isoform usage, or novel transcripts. 
* __Design__: High sequencing depth for robust isoform detection. 
* __Library Prep__: Use protocols that preserve strand information (e.g., stranded RNA-seq).
* __Sequencing__: Paired-end sequencing with ~150bp reads.
* __Downstream Analysis__:
    * Tools: StringTie, Cufflinks, or IsoQuant.
    * Focus: Identify isoforms, splicing junctions, and exon usage. 

## Single-cell RNA-seq

* __Goal__: Investigate gene expression at the resolution of individual cells. 
* __Design__: Optimise for cell capture efficiency and number of cells. 
* __Library Prep__: Droplet-based methods (e.g., 10x Genomics) or plate-based methods (e.g., Smart-seq2).
* __Sequencing__: Shallow sequencing per cell (~50K-100K reads per cell) but many cells (>10,000).
* __Downstream Analysis__: 
    * Tools: Seurat, Scanpy, or Monocle.
    * Focus: Clustering, trajectory analysis, cell type identification. 

## Small RNA Analysis

* __Goal__: Study non-coding RNAs (e.g., miRNAs, siRNAs).
* __Design__: Small RNA-specific library prep kits.
* __Library Prep__: Capture RNAs with sizes between 15-50bp.
* __Sequencing__: Use short-read sequencing platforms.
* __Downstream Analysis__:
    * Tools: miRDeep, miRBase.
    * Focus: miRNA quantification and novel miRNA discovery. 

## Other Types of RNA-seq Analyses 

* De novo transcriptome assembly 
* RNA editing and modifications
* Fusion gene detection 
* Metatranscriptomics

# General computational steps of RNA-seq workflows

Each type of RNA-seq analysis has distinct requirements and challenges but also a common theme.

1. Quality control
2. Trimming
3. Read alignment
4. Read quantification 
5. Downstream analysis
    * Differential gene expression analysis
    * Functional annotation and enrichment
    * Alternative splicing 
    * Isoform quantification 
    * and more...

### Quality Control

* The first step after obtaining raw data. 
* To check read quality, adapter contamination, GC content, etc.
* Tools: `FastQC`, `MultiQC`

__What is adapter contamination?__

Adapter contamination refers to the presence of adapter sequences in RNA-seq reads that were not fully removed during the library preparation process. Adapters are short, synthetic DNA sequences that are ligated to RNA or cDNA fragments during library preparation to facilitate sequencing. They are essential for binding the fragments to the sequencing platform and allow for amplification and sequencing. However, sometimes fragments of RNA or cDNA are too short, and the sequencing process reads into the adapter sequence. This results in "adapter contamination" in the sequencing data. 

### Trimming 

* To remove low-quality bases and adapter sequences. 
* Tools: `Trimmomatic`, `Cutadapt`, `fastp`

### Read Alignment 

* Align reads to reference genome/transcriptome to identify their origins. 
* Tools: `STAR`, `HISAT2`
* Pseudo-alignment: `Kallisto`, `Salmon` for faster processing.
* Outputs: SAM/BAM files (aligned reads). 

__When to align to genome and when to align to transcriptome?__

Deciding whether to align RNA-seq reads to a genome or transcriptome depends on your research goals and the type of analysis you are performing. 

When to align to the genome:

* __Novel transcript discovery__: If you want to identify new transcripts, novel splice junctions, or previously unannotated regions of transcription, aligning to the genome is essential. Examples such as identifying non-coding RNAs, novel isoforms, or novel genes. 
* __Alternative Splicing Analysis__: For detailed analysis of splicing events (e.g., exon skipping, intron retention), genome alignment is crucial since splice junctions need to be inferred. 
* __Poor Transcriptome Annotation__: If the transcriptome of your organism is imcomplete or poorly annotated, genome alignment ensures that all reads are mapped to their correct locations. 
* __Non-Model Organisms__: For non-model organisms, where a high-quality transcriptome may not exist, aligning to the genome is necessary to extract transcript-level information. 
* __Other Specialised Applications__: Studying intronic/intergenic regions. Detecting unspliced pre-mRNA. 

When to align to the transcriptome:

* __Gene/Transcript Quantification__: If your primary goal is to quantify gene or transcript expression levels (e.g., for differential expression analysis), aligning to the transcriptome is typically faster and sufficient. 
* __Well-Annotated Transcriptomes__: When working with well-studied organisms that have high-quality transcriptome annotations, transcriptome alignment can simplify the analysis.
* __Pseudoalignment Methods__: For fast, lightweight quantification (e.g., tools like Salmon or Kallisto) that use pseudoalignment, transcriptome alignment is standard. 
* __Single-Cell RNA-seq__: Single-cell RNA-seq analysis often uses transcriptome alignment or pseudoalignment, as the focus is typically on transcript quantification rather than novel discovery. 

__What is pseudo-alignmend and why Kallisto and Salmon are faster?__

Pseudoalignment is a computational method used in RNA-seq analysis to assign reads directly to transcripts or genes without performing traditional base-by-base alignment to a reference genome or transcriptome. Instead of aligning reads precisely to a sequence, pseusoalignment focused on identifying the compatibility of a read with one or more transcripts. 

How does pseudoalignment work?

* __Indexing the transcriptome__: A database of transcript sequences is preprocessed to create an index, which maps k-mers (short subsequences of a fixed length) to transcripts. Tools like Kallisto or Salmon use this approach. This step builds an efficient data structure, such as a k-mer hash table or a de Bruijn graph.
* __Read querying__: For each RNA-seq read, its k-mers are extracted and compared to the indexed transcriptome. Instread of finding the exact position where the read aligns, pseudoalignment determines which transcripts contain the k-mers in the read. 
* __Transcript compatibility__: The method determines a set of compatible transcripts for each read - i.e., the transcripts where the read could have originated. This avoids computing an exact alignment but still provides information about transcript-level expression. 
* __Quantification__: Using a statistical model, pseudoalignment tools estimate the abundance of each transcript based on the number of reads compatible with it. 

### Read Quantification

* To count reads that mapped to genes or transcripts to generate a count matrix. 
* Tools: `HTSeq`, `featureCounts`

__When do I need to do read quantification? Is it every time for every type of analyses?__

Read quantification is the process of counting how many sequencing reads map to particular gene, transcript, or genomic feature. These counts represent the expression levels of genes or transcripts and are essential for comparing expression across samples.

When do you need read quantification?

* Differential expression analysis (DEA)
* Gene/transcript expression profiling
* Single-cell RNA-seq analysis
* Gene set enrichment analysis (GSEA)

When is read quantification not necessary?

* Novel transcript discovery
* Alternative splicing analysis
* Variant calling
* Fusion transcript detection

### Normalisation

* Normalise counts to correct for sequencing depth and library size.
* Tools: Built into `DESeq2`, `edgeR`, or `limma`.

### Differential Gene Expression (DGE) Analysis

* To identify genes that are differentially expressed between experimental conditions. 
* Tools: `DESeq2`, `edgeR`, or `limma`.

### Functional Annotation and Enrichment

* Perform Gene Ontology (GO) or KEGG pathway analysis to interpret biological significance. 
* Tools: `DAVID`, `clusterProfiler`, `GOseq`.

__What is gene ontology and KEGG pathway? What are their application scenarios?__



### Visualisation

* Visulaise results with PCA, heatmaps, volcano plots, and more. 
* Tools: `ggplot2`, `pheatmap`, `Seurat`.

### Other Advanced Analyses

* A


...

# title ????

...............

## What is Sequencing Depth?

...

## How long are transcripts?

...

# Data 

