# Before class

This workshop is developed based on the [RNA-seq Bioinformatics Course](https://rnabio.org/) by the [Griffith Lab](https://griffithlab.org/) at Washington University, with the help from [ChatGPT](https://chatgpt.com/).

# Introduction to RNA-seq (Overview)

## What is RNA-Seq Analysis?

RNA-Seq (RNA sequencing) is a powerful technique for studying the transcriptome, the complete set of RNA scripts produced by a genome at a given time. It allows researchers to identify and quantify RNA molecules, providing insights into gene expression, alternative splicing, non-coding RNAs, and other molecular phenomena. 

__Attention!! - RNA-seq is not working on RNA directly but cDNA, here's why:__

* __RNA isolation:__ The first step in RNA-seq is extracting RNA from the sample, which consists of various types of RNA, such as mRNA, rRNA, and non-coding RNAs. 
* __cDNA synthesis:__ Since sequencing technologies generally work with DNA, the extracted RNA is reverse-transcribed into cDNA using the enzyme reverse transcriptase. This is because the sequencing platforms (like Illumina, for example) read DNA sequences, not RNA. 
* __Sequencing:__ The cDNA, which mirrors the RNA sequence (except for the replacement of uracil with thymine), is then fragmented, and libraries are prepared for sequencing. 

## Why is it used in biological research?

...Not discussing it here, but it's something you need to know.

...More like what types of biological questions you can solve by using RNA-seq.

...Including gene expression profiling, discovery of new transcripts, identification of alternative splicing events, comparative transcriptomics, single-cell transcriptomics, functional genomics, comprehensive and unbiased data, studying non-coding RNAs, and understanding disease mechanisms etc. 

...Not discussing why and how RNA-seq can solve these questions. 

## Differences between bulk RNA-seq and single-cell RNA-seq

### Bulk RNA-seq

* RNA is extracted from __a population of cells__, often from a tissue or homogenised sample. The RNA represens the collective gene expression of all the cells in the sample.
* The data generated is a __composite view__ of gene expression, averaging across the entire sample. 
* It provides an overall view of gene expression in the entire sample, but it __lacks information__ about gene expression in individual cells. This means that any heterogeneity (variations between individual cells) is averaged out. 
* It is ideal when you want to compare expression levels across different conditions or groups (e.g., healthy vs. diseased tissues).
* The data is more straightforward because it reflects the total gene expression profile of a population. The analysis is less complex as it doesn't require resolving individual cell-level variation. 
* Statistical analyses generally focus on __differential gene expression__ between conditions, identifying genes that are upregulated or downregulated. 
* Typically, fewer sequencing reads are required, and the data is easier to generate because you're working with a pool of cells. 
* The cost and complexity of sequencing are generally __lower__ compared to single-cell RNA-seq. 
* It provides an average expression profile for the entire sample, so it __misses rare cell types__ or subpopulations that may only make up a small fraction of the tissue. 
* Rare cell types may contribute to the overall expression profile, but they are not distinctly identified. 
* Used for more general, large-scale gene expression analysis. Common applications include:
    * Comparing gene expression between different conditions (e.g., cancer vs. normal tissues).
    * Identifying differentially expressed genes (DEGs) across experimental groups. 
    * Transcriptomic profiling in diseases, developmental stages, and responses to treatments. 
* The interpretation is more straightforward, focusing on differences in average gene expression between experimental groups. 

### Single-cell RNA-seq

* Single-cell RNA-seq focuses on analysing gene expression at the __individual cell level__. Each cell is isolated and its RNA is captured and sequenced separately. 
* The method provides a __granular view__ of gene expression, revealing cell-to-cell variability within a population.
* Single-cell RNA-seq allows for the detection of __gene expression differences__ between individual cells, revealing cellular __heterogeneity__ within a population. This is particularly useful for studying complex tissues where cells may have distinct roles or behaviours, such as tumors, immune cells, or neuronal populations. 
* The data is much more complex due to the higher dimensionality of the dataset. Since it provides expression profiles for each cell, there's a need for specialised analysis techniques to deal with __noise__, sparsity (many genes may not be detected in each cell), and __cell clustering__.
* It often involves techniques like __dimension reduction__ (PCA, t-SNE, UMAP), clustering (to identify subpopulations), and __trajectory analysis__ (to study developmental pathways or differentiation). 
* Each cell's RNA must be captured and sequenced separately, which significantly increases the cost and complexity of the experiment. 
* It also requires special methods for isolating and capturing individual cells (e.g., droplet-based systems, microfluidics).
* While sequencing depth per cell is lower due to the need to capture many cells, the overall __sequencing effort__ is higher. 
* Single-cell RNA-seq is particularly useful for identifying __rare cell populations__ that may be important in diseases, development, or responses to treatment. 
* It can resolve subtle differences in gene expression that may be lost in bulk analysis. 
* More specialised, often used for:
    * Studying __cellular heterogeneity__ within a population (e.g., tumours, immune cell subsets, neuronal diversity).
    * __Identifying rare cell types__ or states, such as stem cells or differentiating cells. 
    * __Developmental biology__, where understanding the transitions between cell states over time is critical. 
    * Investigating __gene expression dynamics__ at the single-cell level during processes like differentiation, immune responses, or disease progression. 
* Data interpretation is more challenging and involves clustering cells based on expression profiles, examining the relationships between cells, and understanding the underlying biological processes at a much higher resolution. 

| Feature                | **Bulk RNA-seq**                            | **Single-cell RNA-seq**                      |
|------------------------|---------------------------------------------|---------------------------------------------|
| **Sample Composition**  | Population of cells (averaged data)        | Individual cells (high-resolution data)     |
| **Gene Expression**     | Average expression across the sample        | Cell-to-cell variability, reveals heterogeneity |
| **Data Complexity**     | Less complex (bulk analysis)                | High complexity (individual cell analysis)  |
| **Rare Cell Types**     | May miss rare cell populations              | Can identify rare cell populations          |
| **Applications**        | General gene expression comparison          | Identifying cell subtypes, dynamics, and heterogeneity |

# Before the analyses

Should write about how the experimental design, library prepare, and sequencing technology used can influnce our downstream analyses, such as the use of different pipelines and software, quality control etc.... 

There are many RNA-seq library construction strategies. 

# General Steps of RNA-seq Analysis

1. 
2. 
3. 
4. 
5. 

# Types of RNA-seq Analysis  

RNA-seq is a versatile technology with diverse applications, and the specifics of experimental design, library preparation, sequencing technology, and downstream analysis vary significantly based on the research goal.

We will focus on the easiest one in this workshop, differential gene expression (DEG) analysis.

## 1. Differential Gene Expression (DEG) Analysis

* __Goal__: Compare gene expression levels between conditions (e.g., treated vs. untreated, diseased vs. healthy).
* __Design__: Biological replicates are critical. 
* __Library Prep__: Typically uses poly(A) enrichment or rRNA depletion for mRNA-focused studies. 
* __Sequencing__: Moderate depth (20-50M reads per sample) is sufficient. 
* __Downstream Analysis__: 
    * Alignment tools: STAR, HISAT2.
    * Differential gene expression tools: DESeq2, edgeR.

## 2. Transcript Isoform Analysis

* __Goal__: Detect alternative splicing, isoform usage, or novel transcripts. 
* __Design__: High sequencing depth for robust isoform detection. 
* __Library Prep__: Use protocols that preserve strand information (e.g., stranded RNA-seq).
* __Sequencing__: Paired-end sequencing with ~150bp reads.
* __Downstream Analysis__:
    * Tools: StringTie, Cufflinks, or IsoQuant.
    * Focus: Identify isoforms, splicing junctions, and exon usage. 

## 3. Single-cell RNA-seq

* __Goal__: Investigate gene expression at the resolution of individual cells. 
* __Design__: Optimise for cell capture efficiency and number of cells. 
* __Library Prep__: Droplet-based methods (e.g., 10x Genomics) or plate-based methods (e.g., Smart-seq2).
* __Sequencing__: Shallow sequencing per cell (~50K-100K reads per cell) but many cells (>10,000).
* __Downstream Analysis__: 
    * Tools: Seurat, Scanpy, or Monocle.
    * Focus: Clustering, trajectory analysis, cell type identification. 

## 4. Small RNA Analysis:

* __Goal__: Study non-coding RNAs (e.g., miRNAs, siRNAs).
* __Design__: Small RNA-specific library prep kits.
* __Library Prep__: 




# ---------------------------------


### 1. Experimental Design

* __Biological Replicates:__ Ensure enough replicates for robust statistical analysis. 
* __Treatment Conditions:__ Plan comparisons (e.g., control vs. treated).
* __Sequencing Depth:__ Decide the number of reads needed based on the complexity of the transcriptome. 

### 2. RNA Extraction and Library Preparation

* __RNA Extraction:__ Isolate high-quality RNA from cells or tissues. 
* __Quality Control:__ Assess RNA integrity (e.g., using a Bioanalyser).
* __Library Preparation:__ Convert RNA to cDNA, fragment it, and add sequencing adapters. 

### 3. Sequencing

* Use high-throughput sequencing platforms (e.g., Illumina) to generate millions of short reads. 


......

# Common Analysis Goals 

* Gene expression and differential expression
* Alternative expression analysis
* Transcript discovery and annotation
* Allele specific expression (relating to SNPs or mutations)
* Mutation discovery
* Fusion detection
* RNA editing

# General themes of RNA-seq workflows

* Each type of RNA-seq analysis has distinct requirements and challenges but also a common theme:
    1. Obtain raw data
    2. Align/asemble reads
    3. Process alignment with a tool specific to the goal (e.g., `cufflinks` for expression analysis, `defuse` for fusion detection, etc.)
    4. Post process (import into downstream software, e.g., R and Matlab)
    5. Summarise and visualise (create gene lists, prioritise candidates for validation etc.)

# --------------------------------------------------------

# 2. Bulk RNA-seq Analysis Workflow

* Preprocessing: Quality control (fastqc, trimming).
* Alignment: Mapping reads to a reference genome (STAR, HISAT2).
* Quantification: Counting reads mapping to genes (featureCounts, HTSeq).
* Differential Expression Analysis: Tools like DESeq2, edgeR.
* Visualisations: MA plots, heatmaps, PCA.

# 3. Single-cell RNA-seq Analysis Workflow

* Special considerations for single-cell data (e.g., dropout events, sparsity).
* Preprocessing: Filtering cells, normalisation.
* Dimensionality Reduction: PCA, t-SNE, UMAP.
* Clustering and cell type identification.
* Differential Expression in single cells: Tools like Seurat or Scanpy.
* Visualisation: Cluster heatmaps, UMAP plots, trajectory analysis. 

# 4. Hands-On Section 

* Provide a small dataset for analysis (either bulk or single-cell).
* Walk the participants through each step of the analysis (preprocessing, alignment, quantification, visualisation).
* Encourage participants to explore different tools and methods. 

# 5. Challenges and Best Practises

* Common pitfalls in RNA-seq analysis.
* Tips for optimising experiments and analyses. 

# 6. Q&A and Discussion 

* Open the floor for questions. 
* Discuss real-world examples and ongoing research. 

# ...???

## Steps in RNA-Seq Analysis

?

### 4. Quality Control of Raw Data 

* __Tools:__ Use tools like FastQC to check read quality, GC content, and adapter contamination. 
* __Trimming:__ Remove low-quality reads and adapter sequences (e.g., using Trimmomatic or Cutadapt). 

### 5. Alignment to Reference Genome 

* Align reads to a reference genome or transcriptome using tools like STAR, HISAT2, or Bowtie2.
* Generate alignment files in formats like BAM or SAM. 

### 6. Transcript Assembly and Quantification 

* __Transcript Assembly:__ Use tools like StringTie or Cufflinks to reconstruct transcripts. 
* __Quantification:__ Measure gene/transcript abundance with tools like Salmon or Kallisto (quantify without alignment, if desired).

### 7. Differential Expression Analysis 

* Identidy differentially expressed genes (DEGs) between conditions using statistical tools such as DESeq2, EdgeR, and Limma-voom. 

### 8. Functional Analysis 

* Perform Gene Ontology (GO) and pathway enrichment analyses to interpret the biological significance of DEGs. 
* Use tools like DAVID, GSEA, or Enrichr. 

### 9. Data Visualisation 

* Common plots include:
    * Heatmaps (gene expression patterns).
    * Volcano plots (DEGs).
    * PCA (sample clustering and variability).

### 10. Validation 

* Confirm key findings with independent methods, such as qPCR or Western blotting. 

## Applications of RNA-seq

* __Gene Expression Profiling:__ Study differences in gene expression between conditions. 
* __Alternative Splicing:__ Identify isoform usage and alternative splicing events. 
* __Non-Coding RNAs:__ Detect and quantidy long non-coding RNAs (IncRNAs) and microRNAs. 
* __Mutation Detection:__ Discover RNA editing and single nucleotide variants (SNVs).
* __De Novo Transcriptome Assembly:__ Assemble transcriptomes without a reference genome (useful for non-model organisms). 

## Common Challenges 

* __Data Volume:__ Managing and analysing large datasets can be computationally intensive. 
* __Batch Effects:__ Variability introduced by technical or biological factors must be controlled. 
* __Annotation Accuracy:__ Results depend on the quality of the reference genome and annotation. 

## What is Sequencing Depth?

...

## How long are transcripts?

...

# Introduction to Linux Commands

...

# Data 

