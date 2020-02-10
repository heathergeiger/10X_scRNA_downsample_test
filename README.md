# Introduction

Single-cell RNA-sequencing is a rapidly growing field that allows one to see which genes are expressed in subsets of a tissue or particular cell types rather than just in aggregate.

10X Genomics has a very useful protocol for performing single-cell RNA-seq at high throughput. This protocol may be combined with cell hashing, which uses an oligonucleotide to label each cell with its sample of origin, to sequence cells from multiple different samples all at once (Stoeckius et al., 2018).

With single-cell RNA-seq plus cell hashing, the number of total cells is relatively fixed at 20,000 total cells. Therefore rather than trying to trade off the number of cells vs. number of reads per cell given a fixed sequencing cost, the question becomes how much we can minimize sequencing cost while still determining the same biology from those 20,000 cells.

Another consideration is maximimizing the use of the flow cell. An Illumina S4 flow cell for the NovaSeq sequencer is estimated to produce 2-2.5 billion fragments per lane [(Illumina)](https://www.illumina.com/systems/sequencing-platforms/novaseq/specifications.html). However, multiplexing allows these reads to be divided among multiple pools of 20,000 cells. Given a large number of cells to process, it would be beneficial to maximize the number of pools that can be run on one lane.

Combining 4 pools per lane would give an average number of reads per cell of at least 25,000 for a cell hashing experiment (2 billion/(20,000x4)). This is slightly higher than 10X Genomics minimum recommendation of 20,000 reads per cell [(10X Genomics)](https://kb.10xgenomics.com/hc/en-us/articles/115002022743-What-is-the-recommended-sequencing-depth-for-Single-Cell-3-and-5-Gene-Expression-libraries-). We were curious if we could potentially increase the number of pools per lane even further, perhaps to say 10 pools per lane (for an average number of reads per cell of 10,000).

In this post, we seek to use published data to compare cell type identification from a very deeply sequenced pool of cells (over 75,000 average reads per cell) to various levels of more shallow sequencing. For this analysis, we assume that random downsampling of the raw reads (FASTQs) provides an adequate simulation of lower sequencing depth. Then, we compare clustering results from Seurat to quantify how similar or distinct cell type identification is between different levels of sequencing.

A major research area of the New York Genome Center is neuropsychiatric and neurodegenerative disease, so we chose a dataset of neurons to best match the single-cell or single-neuron samples that we might encounter as part of this research in the future. The dataset we used contains over 5,000 cells from a combined cortex, hippocampus, and subventricular zone of an E18 mouse. Libraries were prepared using 10X Genomics latest "Next GEM" kit. This data is available for download in the "datasets" section of the 10X Genomics website [(5k mouse neurons, Next GEM kit)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_neuron_v3_nextgem?src=search&lss=google&cnm=sem-goog-2020-website-page-ra_g-p_brand-amr&cid=7011P000000oWx1).

# Analysis

## Downsampling the data



# References

Stoeckius, M., Zheng, S., Houck-Loomis, B. et al. Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome Biol 19, 224 (2018). https://doi.org/10.1186/s13059-018-1603-1
