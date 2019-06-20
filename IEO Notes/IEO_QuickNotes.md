# Information Extraction from OMICS Technologies

* [Quality Assessment and Normalization of RNA\-seq Data](#quality-assessment-and-normalization-of-rna-seq-data)
* [Experimental Design and Batch Effect Identification](#experimental-design-and-batch-effect-identification)
* [Differential Expression Analysis (I)](#differential-expression-analysis-i)
* [Differential Expression Analysis (II)](#differential-expression-analysis-ii)
* [Differential Expression Analysis (III)](#differential-expression-analysis-iii)
* [Pulling Functional Annotations with R/BioC](#pulling-functional-annotations-with-rbioc)
* [Reproducible Research](#reproducible-research)
* [Functional Enrichment](#functional-enrichment)
* [GSEA and GSVA](#gsea-and-gsva)
* [QA, Mapping & Summarization of HTS Reads](#qa-mapping--summarization-of-hts-reads)
* [Analysis of Metabolomics Data](#analysis-of-metabolomics-data)

## Quality Assessment and Normalization of RNA-seq Data

* RNA-seq: relative concentration of RNA molecules.
* Preprocessign:
	* Align to genome, transcriptome or assemble directly; summarized at gene or transcript level.
* Result: matrix of counts. Rows: genes (entrez ID), columns: samples.
* Gene-level counts can be either directely obtained from sequence reads mapped to the genome or by combining transcript-level counts for the corresponding gene.
* Diagnostics: library size: coverage.
* Normalization needed because of different sequencing depths:
	* Within-sample: compare across features in a sample (CPM).
	* Between-sample: compare features across samples.
		* TMM: sample-specific normalization factors. Scaling factor for each, to adjust its different RNA composition.
		* Quantile normalization.
* Diagnostics:
	* Distribution of expression across samples to see samples with distinctive RNA composition.
	* Distribution of expression levels across genes: identify lowly expressed genes.
		* Remove by simple cutoff or cutoff+gene in a minimum number of samples.
* Assess normalization by MA plots:
	* Compare two samples or two groups of samples.
	* We can see if more genes are up or down (unexpected). This is due to different sample proportions for each gene. Correct by normalization procedure that takes into account sample-specific effects (TMM).
	* We also see high fold changes at low expression levels: fold changes from high expression are more reliable.
* More on between-sample normalization:
	* Sample preparation induces differences between samples (different amounts of RNA, sequencing conditions).
	* Quantile normalization: map every value in one sample to a corresponding quantile of a reference distribution.
* MDS plots: help to identify samples with distinctive features from the rest.

## Experimental Design and Batch Effect Identification

* A good experimental design should follow:
	* Replication: repeat measurements (to minimize variation).
	* Randomization: subjects are allocated to groups randomly to control confounding effects.
	* Blocking: balance the groups.
* Bad experimental designs have confounding factors.
* Batch effect: sub-groups that are differentiated by other reason than the biological variables.
* Normalization can remove part of batch effect.
Block design allows to identify and correct batch.
* Worst case: batch correlated to the outcome of interest -> wrong conclusions.
* Surrogates of batch effect:
	* Microarrays: timestamp.
	* RNA-seq: sample information.
* Identification: hierarchical clustering.
	* Start by prior count large values (to moderate extreme fold-changes derived from low counts).
	* Calculare the distances between samples with non-parametric approaches (i.e. Spearman correlation).
	* Perform the hierarchical clustering.
	* Assessment: dendogram or MDS.
* Batch effect can be:
	* Adjusted: the choice if statistical testing will be applied (inference). Then, it can be included in a linear model.
		* SVA: captures sources of unknown heterogeneity. Gives “surrogate variables” with continuous values, and requires:
			* The matrix of counts.
			* Full model matrix: outcome of interest and adjustment terms.
			* Null model matrix: adjustment variables only.
		* Performs F-tests to compare genes changing between null and full models. Not really addequate because of relationships between mean and variance. P-values should have a uniform distribution.
	* Removed: appropriate for exploratory analysis (clustering).
		* It can remove part of the biological variability.
		* Can lead to spurious differences between groups of the outcome of interest if the distribution between batches is not balanced.
		* it can be used for inferential analysis if adjusting doesn’t work, but this is the last resource.
		* ComBat: empirical Bayes method robust to outliers in small sample sizes.
		* QR decomposition: requires a batch indicator.
		* SVD: very agressive, decomposes the gene expression matrix and removes the largest variability (assumes it is the batch).

## Differential Expression Analysis (I)

* Differential expression: comparing levels of gene expression between two conditions.
* Start from normalized data.
* We will have an outcome of interest (the conditions we want to compare, but experimental design variables might be a blocking factor. We should see if those are balanced in our outcome of interest).
* Large count values have more variability than smaller ones. Base 2 logarithms start to stabilize the variance.
* Basic DE comparison: gene by gene between the conditions.
	* For each sample in a group, take the mean of all the counts of each gene.
	* Means between groups are compared, obtaining the fold-change (base 2 logarithmic scale). It is a ratio of the mean expression values of the two conditions.
	* Then, we can get a list of the most differentially expressed genes by ranking them by the absolute value of log2 fold change.
* A problem comes: maybe some of those genes are not DE. We need to verify them. There are two options:
	* Use an independent assay (costs money).
	* Replicate: have different replicates per condition and use statistical inference.
* Statistical significance of fold changes derives from variability of expression values within each condition.
* Null hypothesis: no difference in expression for the population means.
* Test the hypotheses with a t-statistic, in which we compare the differences between group means with the differences within the group.
* The null hypothesis can be wrongly accepted or rejected. When performing multiple tests, the probability of rejecting the null hypothesis (type I error) increases, even if it is true. So multiple testing corrections are applied:
	* Bonferroni: really conservative, 
	* FDR
	* Non-specific filtering: discarding genes using a criterion that is independent of the test statistic (genes not likely to be DE, or genes with no biological interest). This allows to reduce the number of multiple tests. Examples of criteria are absence of functional annotation or low expression levels.

## Differential Expression Analysis (II)

**Introduction**

* Linear regression models capture linear relationships between a response variable and a predictor variable. Its notation is:
	* $y$: response/dependent variable
	* $x$: predictor/independent/explanatory variable.
		* If it is discrete, it is called *factor*, and its possible values *levels*.
	* $beta_0$: intercept.
	* $beta_1$: coefficient, effect if $x$ is discrete.
* What we do is building a model, which has systematic structure and random variation. When the variation is not purely random, then more $p$ dimensions are needed.
* Linear regression models can be used for DE analysis:
	* $y_{ij}$: expression values for gene $i$ in sample $j$.
	* $beta_0$: intercept, representing an unknown baseline expression level.
	* $beta_j$: unknown coefficient/parameter representing the effect of experiment $j$.
	* $epsilon_{ij}$: error or unexplained variation of the gene in that sample.
* In this models, a matrix is built, called **design matrix**. This is done for algebraic reasons, and software builds it from the definition of the model. The matrix contains rows = samples and columns = coefficients to estimate.
	* Confounding variables add columns.

**DE analysis with limma**

* The basic approach has 4 steps:
	1. Build design matrix: `model.matrix()`. This is a base R function for linear modeling.
	2. Fit linear model: `lmFit() `, which estimates the unknowns. Its inputs are the design matrix and the data.
	3. Calculate moderated t-statistics: `eBayes()`. This decides what to test (in complex models, more stuff is tested).
	4. Output results: `topTable()`.

* How to use limma with batch adjustment is done in slide 12 with a toy example.
* In the step 4, what we want is the second column (extracted with the argument `coef = 2`. Those are the p-values.
* Then, we check that they are well calibrated under the null hypothesis: we expect an equal distribution of p-values (the expected when there is no signal).
* We also observe that batch effect removal results in more false positives.
* The step 1 uses a reference. By default, it is the first factor it finds, but this can be changed with the `ref` argument.
* The function `decideTests()` can be used to obtain a summary of the results (slide 20). It uses 5% FDR by default, but this can be changed.

**Some analysis of results**

* The step 4 gives entrex identifiers. This might not be very explanative, so slide 22 explains how to add more gene data from the SE object (such as the chromosome or the gene symbol).
* We can also do an analysis of chromosome distribution of DE genes (slide 23).
* Slides 23 and 24 explain how to assess the accuracy of the analysis: use the entrez identifiers of genes from sex chromosomes. We should get a list of lung cancer for this (as it is our outcome of interest).
* We can also study the distribution of p-values and moderated t-statistics. We expect a flat distribution (and I think a diagonal for the t).

**Adjusting for known covariates**

* We just add them to the model.

**Adjusting for unknown covariates**

* We can do this with SVA (slide 32). First, we estimate them and then we add them to the model. This results in a bigger design matrix.

**Other diagnostics**

* If we consider using more than one model, we can compare the performance graphically (slide 38).
* We also perform volcano and MA plots.
* p-value distribution can reveal the presence of confounding factors and heterogeneity
* Quantile-quantile (Q-Q) plots and volcano plots give a graphical overview of the amount of DE in our data.
* Gold-standard subset of genes with documented differential expression can help to tune the analysis. This subset, however, is most of the times not available

## Differential Expression Analysis (III)

**Other issues to address**

* Mean-variance relationship in log-transformed counts from RNA-seq data (voom).
* Adjust for repeated measurements.
* Factorial and paired designs.

**Mean-variance relationship**

* Identical CPM values for different count sizes happen because:
	* RNA-seq counts vary from sample to sample for the same gene.
	* Different samples can be sequenced to different depths.
* This relationship is incorporated to the model:
	* `limma-trend`
	* `limma-voom`
* If the library size is similar, their precision and recall is similar. When it is similar, `voom` is superior to `trend`.

**Adjusting for repeated measurements**

* Repeated measurements lead to overly low p-values and false positive calls.
* Adjusting for this keeps the statistical power.
* This would be easier to adjust if all individuals were sequenced twice.
* Existence of duplicates can be checked from the cell-line phenotypic information.
* Strategies:
	* Averaging replicates.
	* Including a replicated indicator as a main effect.
	* Using a linear mixed-effects model.

**Paired measurements**

* Add the individual identifier to the design matrix.
* Samples without pair can be discarded or apply the repeated measurements strategy.

**Factorial designs**

* The outcome of interest has two or more factors.
* A factor from the combination of the factors can be used.

## Pulling Functional Annotations with R/BioC

* Functional annotations: descriptions of biological activity of structure elements of the genome and features profiled with HT platoforms.
* In this slides: what is reported on databases on the role of the sequence.
* They go from molecular mechanisms to high-level phenotypes.
* Varying coverage among organisms. Dynamic nature.

**Gene-centric packages**

* Categorized into groups (chip annotations, organism annotations or genome centric annotations).

**Genome-centric**

* Functional annotations mapped to genomic coordinates.
* Types:
	* Exon-level or transcript level.
	* SNP locations
* Usage cases:
	* Extracting genome-wide coordinates of the protein-coding sequence (CDS) of genes and their corresponding nucleotide sequences.
	* Retrieving promoter regions of genes.
	* Select genes overlapping regions of interest such as those mapped by short-reads.

**Other access**

* Packages discussed before do not contain all possible information.
* Other packages allow access to complete databases (UCSC Genome Browser, Biomart).

## Reproducible Research

* While funding in genomics increases, it doesn’t mean that results are proportional.
* Ten years ago, reproducibility started to be a concern.
* Lack of reproducibility can come from unclear methods or unavailable software.
* Methodological quality of scientific experiments does not increase with increasing rank of the journal.
* Retractions:
	* Public statement about an earlier statement which amends or cancels the original one.
	* Retractions are still small, but they increase exponentially, correlating positively with the impact factor.
	* Most retractions are due to misconduct.
* Known cases have been spotted after some “unrelated” information came, but they could have been spotted earlier if the code and exact steps were released with the publication.
* **Reproducibility:** start from the same samples/data, use the same methods, get the same results.
	* Always feasible with computationally-oriented research.
* **Replicability:** using different data, get confirmatory results.
	* It is a sum of replicability and conducting the experiment again.
	* It might be challenging in epidemiology (recruit again a cohort) or molecular biology (complex cell manipulation).
		* For this reason, reproducibility is the minimum standard.
* Reproducibility is not only about the data, but also about having software and specifically those particular versions) available. That is why cloud computing and Docker are important.

## Functional Enrichment

* Functional enrichment is used with DE genes to study the molecular processes they are involved.
* Biological pathway: series of actions among molecules in a cell that leads to certain product.
* Genes act collectively under finely tuned molecular regulatory programs.
* DE pathways can be searched within those that are defined (databases such as Reactome or KEGG).
* Pathways can be defined as gene sets, which can be obtained from Gene Ontology or MSigDB.
* Some pathway definitions don't include the ones of interest, which requires manual curation or inference.
* Pathways can be already seen visually (volcano plots), but it is possible to obtain them statistically.
* DE pathways can be found by functional enrichment:
	* Search DE genes.
	* For each set, see how many of the DE genes are there.
		* The set is enriched if it contains more DE genes than the expected by pure chance.
		* A hypergeometric test (Fisher's exact test) can be applied.
* Fisher's exact test:
	* Null hypothesis: a gene set does not include more DE genes than the expected by chance.
	* Gene universe: all considered genes.
	* The number of DE genes follows an hypergeometric distribution under the null hypothesis.
	* This test is used in the analysis of contingency tables that result from cross-classifying objects by two factors (in this case, membership to a set and DE status).
	* The lower the p-value the less consistent with chance observing those k genes.
* Gene Ontology (GO) analysis: apply the one-tailed Fisher's exact test to every gene set in the GO database.
	* All terms belong to three ontologies.
	* Terms are related to each other by parent-child relationships (from general to specific). 
	* GO terms are annotated to genes
* GO analysis with `GOstats`:
	* Parameter object with information specifying the gene universe, the set of DE genes, the annotation package to use.
	* Run the functional enrichment analysis.
	* Store and visualize results.
* GO analysis can have redundancy because it is hierarchical and some genes can be in the parent and child terms.
	* Conditional testing solves this: the child term is considered more relevant.
* Finally, terms with few or a lot of genes are not informative, so it is better to filter them out.

## GSEA and GSVA

* Classical FE based on Fisher's exact tests require a list of DE genes.
* In some contexts, this list is too small or no significant DE genes exist.
* This doesn't mean that there aren't molecular mechanisms changed: maybe genes change little but coordinately.
* GSEA adresses this: it looks for differential expression at the gene level.
	* There are multiple methods based on this idea.
* Main concepts:
	* An enrichment score is calculated for each gene set as function of the changes in gene expression of the genes forming the set.
* Data is prepared by loading molecular data and gene sets to one object.
* Simple GSEA is parametric. It has a distribution of t-statistics that is approximately normal. This allows to calculate a z-score.
* Higher Z-score: higher enrichment.
* As the z-scores follow a normal distribution, a z-test can be done by calculating p-values. Multiple tests adjustmen is needed because of gene set overlaps.
* We can also calculate the gene set overlap.
* With z-scores, sometimes the shift inmean expression can be cancelled out when half the genes are up and half are down. Chi2 is adequate to detect this (it is called change in scale).
	* For gene sets > 20 genes the chi2 follows a normal distribution.

**GSVA**

* There are different types of GSEA approaches:
	* Tested null hypothesis:
		* Competitive: involves genes in and out the set.
		* Self-contained: only genes in the set.
	* ES calculated using phenotypes or not:
		* Supervised: uses phenotypes.
		* Unsupervised: doesn't use phenotypes.
	* Units quantifies by ES:
		* Entire smple population: ES tells something about a gene set and the entire dataset.
		* Single sample: ES tells something about a gene set and an individual sample.
* GSVA: unsupervised and single sample. Non parametric, with a competitive-like approach.
* Allows pathway-centric analyses of gene expression

**Conclusions**

* GSEA detects small but consistent changes in expression within gene sets and facilitate the interpretability of results.
* Simple GSEA: allows one to detect shifts in mean expression and changes in scale. It assumes, however, that genes are expressed independently.
* Non-specific filtering methods for gene sets are pretty much limited to constraining their minimum and maximum size. Other strategies that could remove irrelevant gene set definitions upfront would likely decrease false positive calls and increase statistical power.

## QA, Mapping & Summarization of HTS Reads

* HTS instrument produce sequencing reads.
	* Those are raw data.
	* Discrete nature (opposed to continuous fluorescence from microarrays).
	* Prone to technical variability and biases.
		* Some a re unknown, as the technology evolves.
	* Stored in FASTQ, with quality encoded in ASCII characters (platform-specific meanings).
	* Quality value: integer mapping of the probability of the base call being incorrect. Current standard: Phred.
	* QA is needed to filter out low-quality reads or reads that are not useful.
* Possible QA diagnostics:
	* Nucleotide distribution throuhg all reads independently at each cycle.
		* Obtain the numbers of each nucleotide for each cicle and divide by the total nucleotides of that cycle.
	* Mean quality per cycle.
		* Coherce quality scores to a matrix of integers.
	* Average quality distribution per read (using a random sample of reads).
	* Filter out low quality reads. Filter by distance to some specific reads by obtaining the distribution of distances to decide a cutoff.
* Read mapping:
	* Align them to a reference genome. Specific tools exist for short reads.
	* BWA: used to call variants from aligned reads of DNA-seq experiments.
	* STAR: to obtain expression profiles from aligned reads of RNA-seq experiments.
	* Formats: SAM is the *de facto* standard. BAM is its binary representation.
* Summarization into features of interest (exons, transcripts, genes):
	* Approaches: 
		* Count reads falling into knownfeatures.
		* Assemble reads to identify novel features and count reads per feature again.
	* `GenomicFeatures` and `GenomicAlignments` packages. 
	* Starting point: BAM files.
	* Annotations are needed to summarize by counting.
	* Overlaps are summarized and a `SummarizedExperiment` object is returned.
	* There are different strategies for the summarization if a read does not map an entire feature: union, intersection strict and intersection non-empty.



## Analysis of Metabolomics Data

**Introduction**

* Genomics talks tells what could happen.
* Transcriptomics tells what appears to happen.
* Metabolomics tells what makes it to happen.
* Metabolomics tells what has happened.

**Definition**

* Metabolites are low weight chemical products used by organisms to sustain life.
* Metabolomics studies them. Two types:
	* Targeted, for known metabolites.
	* Untargeted: metabolite profiling.
* The metabolome is the phenotype of a cell. It can be used:
	* To find biomarkers.
	* To understand diseases.
	* To predict the outcome of a treatment.

**Procedure**

* An HPLC system uses solvents to load the sample, which goes through a chromatography column to separate them.
	* Time taken to travel trhough the column: retention time.
	* Peak intensities are relative to the abundance of each compound.
* Once separated, a mass spectometer ionizes the compounds to obtain its mass/charge ratio.
	*  Peak intensities relative to the abundance of each compound.
* Final signal:
	* Combines intensity, mass/charge and retention time for each observed peak.
	* Sparse matrix.
	* Anisotropic behaviour.

**Analysis workflow**

* For untargeted studies:
	1. Peak Detection: highly dependent on the parameters of the wavelet.
	2. Peak Annotation.
	3. Peak Alignment:
		* Samples show temporal drifts.
	4. Data Normalisation.
		* There are intensity drift effects.
		* Instrumental factors can change the output (batch effects, column change).
		* Chromatographic columns get old.
	5. Statistical Analysis.
		* Unsupervised methods: PCA and hierarchical clustering.
		* Supervised methods: PLSDA, SVM, Random Forest.
		* Univariate statistical tests +multiple test correction are not effective: statistical power decreases because there are many more variables than samples.
	6. Metabolite Identification.
	7. Functional Interpretation and Pathway Analysis

**Functional enrichment**

* Relate the experimental data with biological pathways to provide biological context.
* Database annotations can be automatically generated from literature mining, manually curated or a mixed approach.
* Gene ontology: three ontologies and a hierarchical structure reflecting the specificity of the entry.
* UniProt: protein sequence and annotation data.
* STRING: protein-protein interaction database; automatic and curated annotations
* KEGG: comprehensive and curated pathway database.
* Reactome, BioCyc, WikiPathways, SMPDB.
* Enrichment is needed for functional analysis.
	* Over representation analysis:
		* DE genes are computed.
		* Assess if the number of affected genes in a particular pathway is higher than expected by chance.
* GSEA: takes advantage on quantitative data about gene expression.
* Network-based methods take into account the topological features of biological networks in order to test pathways.
* Metabolomics is younger, so enrichment  concepts come from genomics and transcriptomics.
	* Background for statistical test is not obvious. Most of the metabolome is unknown.
	* Experimental devices can have low sensitivity and noisy signals.
	* No technology to measure the full spectrum of compounds. Large part of the background is not observable.

**Subnetwork analysis**

* Extracting subnetworks (modules) from biological structures can give insights of the affected mechanisms in the condition under inspection.
* Two approaches:
	* Signigicant Area Search:
		* Score nodes, assess an agregate score of subnetworks and build the optimal ones.
	* Difussion processes: Nodes known to be affected introduce flow to the network. Then, maximum flow subnetworks are found.

**More concepts**

* Pathway enrichment techinques bridge between metabolic pathways and LC/MS data.
* This also allows further integration with other omics.
* Relevant parts of graphs are found by applying heat to certain nodes. Then, hottest subgraphs can be analyzed.
	* Temperatures cannot be compared between nodes, and they are not explicative by themselves.
	* Null models are used to analyse the distribution of temperatures for random cases to be able to compare the temperatures.
		* The technology does not allow to measure all compounds from a database, so not all nodes can me measured. A portion of nodes introduces uncertainty.  
		* The model is splitted in two:
			* Observable pool: nodes that can be measured.
			* Latent pool: nodes that cannot be measured.    
		* We want to see how the presence of the latent pool modifies our statistical tests.
		* By introducing the latent pool all the variances grow. The magnitude of this enlargement will assess how sensitive a particular node is to experimental uncertainty.


	