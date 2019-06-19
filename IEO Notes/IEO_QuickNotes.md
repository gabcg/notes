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

* Differential expression is performed by comparing levels of gene expression between two conditions (male-female, healthy-disease, etc).
* Start from normalized data.
* We will have an outcome of interest (the conditions we want to compare, but experimental design variables might be a blocking factor. We should see if those are balanced in our outcome of interest).
* Stabilize the variance with base 2 logarithms. It can be further stabilized.
* Basic DE comparison: gene by gene between the conditions.
	* For each sample in a group, take the mean of all the counts of each gene.
	* Means between groups are compared, obtaining the fold-change (base 2 logarithmic scale).
	* Then, we can get a list of the most differentially expressed genes by ranking them by the absolute value of log2 fold changes.
* A problem comes: maybe some of those genes are not DE. We need to verify them. There are two options:
	* Use an independent assay (costs money).
	* Replicate: have different replicates per condition and use statistical inference.

* Null hypothesis: no difference in expression for the population means.
* Test the hypotheses with a t-statistic, in which we compare the differences between group means with the differences within the group.
* The null hypothesis can be wrongly accepted or rejected. When performing multiple tests, the probability of rejecting the null hypothesis increases, even if it is true. So multiple testing corrections are applied:
	* Bonferroni: really conservative, 
	* FDR
	* Non-specific filtering: discarding genes using a criterion that is independent of the test statistic (genes not likely to be DE, or genes with no biological interest). This allows to reduce the number of multiple tests. Examples of criteria are absence of functional annotation or low expression levels.

## Differential Expression Analysis (II)

## Differential Expression Analysis (III)

