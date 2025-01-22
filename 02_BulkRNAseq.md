# Steps in Bulk RNA-seq analysis

## 1. Experimental Design

* Define the biological question (e.g., comparing gene expression between conditions).
* Plan replicates to ensure statistical power (at least 3 biological replicates per group is recommended). 
* Select an appropriate sequencing depth (typically 20-50 million reads per sample).

## 2. Sample Preparation and Sequencing 

* Extract high-quality RNA from your samples.
* Assess RNA quality (e.g., using an Agilent Bioanalyser for RNA Integrity Number (RIN)).
* Prepare cDNA libraries for sequencing.
* Perform sequencing on a platform (e.g., Illumina) to generate raw reads.

## 3. Quality Control of Raw Reads

* Inspect raw sequencing data using tools like: 
    * __FastQC__: Provides an overview of quality metrics (base quality, GC content, adapter contamination). 
* Trim low-quality bases and remove adapters using tools like:
    * __Trimmomatic__ or __Cutadapt__. 

## 4. Alignment to Reference Genome

* Map the cleaned reads to a reference genome (or transcriptome) using aligners like:
    * __STAR__: Fast and widely used for RNA-seq.
    * __HISAT2__: Efficient for spliced alignments.
* Output: Aligned reads in a BAM file. 

__???What is a transcriptome? How does it different from a genome?__ 

## 5. Post-Alignment Quality Control 

* Evaluate alignment results:
    * Use __samtools flagstat__ to check the percentage of mapped reads.
    * Use __RSeQC__ to access read distribution across genomic features. __WHY???__
* Check for biases (e.g., 3' bias due to degraded RNA). __WHY???__

## 6. Quantification of Gene Expression 

* Count the number of reads mapped to each gene using tools like:
    * __HTSeq__ or __featureCounts__: Count reads based on gene annotation files (e.g., GTF/GFF).
* Output: A __count matrix__, where rows are genes and columns are samples. 

## 7. Normalisation 

* Normalise the count data to account for differences in sequencing depth and gene length. __HOW???__
* Common normalisation methods:
    * __TPM (Transcripts Per Million)__: For comparing gene expression within a sample.
    * __RPKM/FPKM__: Length-normalised, but less commonly used now. 
    * __DESeq2 or edgeR normalisation__: Scales raw counts for differential expression analysis. 

## 8. Differential Gene Expression Analysis 

* Identify genes with significant expression differences between experimental groups. 
* Common tools:
    * __DESeq2__ (R-based): Handles raw counts directly.
    * __edgeR__ (R-based): Suitable for small sample sizes. 
* Output:
    * List of __differentially expressed genes (DEGs)__ with log2 fold changes and p-values. 

## 9. Functional Annotation and Pathway Analysis

* Interpret the biological relevance of DEGs by performing:
    * __Gene Ontology (GO) Enrichment Analysis__: Identify enriched biological processes, molecular functions, or cellular components.
    * __Pathway Analysis__: Map DEGs to pathways using tools like __KEGG__, __Reactome__, or __GSEA__.

## 10. Data Visualisation 

* __Quality Control__:
    * PCA plot: Visualise sample clustering.
    * Heatmaps: Show clustering of samples/genes.
* __Differential Expression__:
    * MA plot: Log fold change vs. mean expression.
    * Volcano plot: Significant vs. log fold change. 
* __Pathway Analysis__:
    * Enrichment bar plots or network diagrams. 

## 11. Validation 

* Validate key findings using an independent method:
    * __qRT-PCR__: Validate differential expression for a subset of genes. 
* Cross-reference with existing datasets or prior research. 

## 12. Reporting 

* Compile findings into a report or publication:
    * Document methods, quality control steps, and statistical analyses.
    * Share data and scripts for reproducibility (e.g., GitHub or public repositories). 

