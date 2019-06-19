# Information Extraction from OMICS Technologies

* [Quality Assessment and Normalization of RNA\-seq Data](#quality-assessment-and-normalization-of-rna-seq-data)
* [Experimental Design and Batch Effect Identification](#experimental-design-and-batch-effect-identification)
* [Differential Expression Analysis (I)](#differential-expression-analysis-i)
* [Differential Expression Analysis (II)](#differential-expression-analysis-ii)
* [Differential Expression Analysis (III)](#differential-expression-analysis-iii)

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

