# Introduction

Single-cell RNA-sequencing is a rapidly growing field that allows one to see which genes are expressed in subsets of a tissue or particular cell types rather than just in aggregate.

10X Genomics has a very useful protocol for performing single-cell RNA-seq at high throughput. This protocol may be combined with cell hashing, which uses an oligonucleotide to label each cell with its sample of origin, to sequence cells from multiple different samples all at once (Stoeckius et al., 2018).

With single-cell RNA-seq plus cell hashing, the number of total cells is relatively fixed at 20,000 total cells. Therefore rather than trying to trade off the number of cells vs. number of reads per cell given a fixed sequencing cost, the question becomes how much we can minimize sequencing cost while still determining the same biology from those 20,000 cells.

Another consideration is maximimizing the use of the flow cell. An Illumina S4 flow cell for the NovaSeq sequencer is estimated to produce 2-2.5 billion fragments per lane [(Illumina)](https://www.illumina.com/systems/sequencing-platforms/novaseq/specifications.html). However, multiplexing allows these reads to be divided among multiple pools of 20,000 cells. Given a large number of cells to process, it would be beneficial to maximize the number of pools that can be run on one lane.

Combining 4 pools per lane would give an average number of reads per cell of at least 25,000 for a cell hashing experiment (2 billion/(20,000x4)). This is slightly higher than 10X Genomics minimum recommendation of 20,000 reads per cell [(10X Genomics)](https://kb.10xgenomics.com/hc/en-us/articles/115002022743-What-is-the-recommended-sequencing-depth-for-Single-Cell-3-and-5-Gene-Expression-libraries-). We were curious if we could potentially increase the number of pools per lane even further, perhaps to say 10 pools per lane (for an average number of reads per cell of 10,000).

In this post, we seek to use published data to compare cell type identification from a very deeply sequenced pool of cells (over 75,000 average reads per cell) to various levels of more shallow sequencing. For this analysis, we assume that random downsampling of the raw reads (FASTQs) provides an adequate simulation of lower sequencing depth. Then, we compare clustering results from Seurat (Stuart & Butler 2019) to quantify how similar or distinct cell type identification is between different levels of sequencing.

A major research area of the New York Genome Center is neuropsychiatric and neurodegenerative disease, so we chose a dataset of neurons to best match the single-cell or single-neuron samples that we might encounter as part of this research in the future. The dataset we used contains over 5,000 cells from a combined cortex, hippocampus, and subventricular zone of an E18 mouse. Libraries were prepared using 10X Genomics latest "Next GEM" kit. This data is available for download in the "datasets" section of the 10X Genomics website [(5k mouse neurons, Next GEM kit)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_neuron_v3_nextgem?src=search&lss=google&cnm=sem-goog-2020-website-page-ra_g-p_brand-amr&cid=7011P000000oWx1).

# Analysis

## Downsampling the raw data

After downloading and unpacking the FASTQ tarball, the following FASTQs were located under sub-directory 5k_neuron_v3_nextgem_fastqs.

	5k_neuron_v3_nextgem_S1_L001_I1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L001_R1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L001_R2_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L002_I1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L002_R1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L002_R2_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L003_I1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L003_R1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L003_R2_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L004_I1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L004_R1_001.fastq.gz
	5k_neuron_v3_nextgem_S1_L004_R2_001.fastq.gz

To subsample these files, we used dcjones' subsample package, "a small, Unix-y program to efficiently subsample lines from a file without replacement" [(subsample package Github)](https://github.com/dcjones/subsample). 

The original data had 78,597 mean reads per cell. Since we were interested in looking at average read numbers per cell of 10,25, and 50,000, this would mean that we would want to test downsample levels of 13, 32, and 64%. We ended up testing 65% instead of 64%, but otherwise we applied these three levels of downsampling.

Count exact number of reads per lane before start actual downsampling. Would run this for each of L001,2,3, and 4.

    zcat 5k_neuron_v3_nextgem_fastqs/5k_neuron_v3_nextgem_S1_${lane}_I1_001.fastq.gz | wc -l | awk '{print $1/4}' > ${lane}_read_count.txt

Unzip all FASTQs under 5k_neuron_v3_nextgem_fastqs.

	for file in 5k_neuron_v3_nextgem_fastqs/*.gz;do gunzip $file;done

After that, here is the shell script code to run downsampling given a type of FASTQ to downsample (I1, R1, or R2), a lane number (L001,2,3, or 4), and a percent downsample (13,32,or 65).

	#Get variables from command line arguments.

	type=$1
	lane=$2
	downsample=$3

	#Get exact number of reads to select for downsample.

	full_lane_reads=`cat ${lane}_read_count.txt`
	num_reads=`echo $downsample $full_lane_reads | awk '{printf("%.0f",($1/100)*$2)}'`

	#Create output directory if it does not exist already.
	
	if [ ! -e downsample_pct_${downsample}_fastqs ]
	then
	mkdir downsample_pct_${downsample}_fastqs
	fi

	#Run subsample.
	#Use seed so that same reads are selected for I1, R1, and R2 from the same lane.
	#Seed 1392 was used for this analysis.

	subsample -k 4 -n $num_reads -s 1392 5k_neuron_v3_nextgem_fastqs/5k_neuron_v3_nextgem_S1_${lane}_${type}_001.fastq > downsample_pct_${downsample}_fastqs/downsample_pct_${downsample}_S1_${lane}_${type}_001.fastq

	#Zip up output file.
	
	gzip downsample_pct_${downsample}_fastqs/downsample_pct_${downsample}_S1_${lane}_${type}_001.fastq

Finally, the full data FASTQs were zipped back up as well.

	for file in 5k_neuron_v3_nextgem_fastqs/*.fastq;do gzip $file;done

## Get gene-cell count matrix for each downsample.

After this, 10X's CellRanger pipeline was used to get the gene-cell UMI count matrix for each downsample from these FASTQs.

Since we used CellRanger v3.1.0 (vs. the 3.0.2 used for the gene-cell matrix available from the 10X website), we also ran CellRanger on the full data.

We used reference version v3.0.0 (mm10) downloaded from the 10X website as well.

This resulted in the gene-cell UMI count matrix being available under the following paths:

	CellRanger_output/subsamp_13/subsamp_13/outs/filtered_feature_bc_matrix
	CellRanger_output/subsamp_32/subsamp_32/outs/filtered_feature_bc_matrix
	CellRanger_output/subsamp_65/subsamp_65/outs/filtered_feature_bc_matrix
	CellRanger_output/full_data/full_data/outs/filtered_feature_bc_matrix

Unfortunately the 10X software does not provide the option to set a seed, so your output from this step may differ slightly from ours given identical raw data input.

## Run Seurat standard analysis on each downsample.

After this, we ran fairly standard Seurat analysis on each downsample level.

For this analysis, we used Seurat v3.0.1.9015.

Would run the code below using arguments "subsamp_13", "subsamp_32", "subsamp_65", and "full_data".

	#Accept downsample level as command line argument.

	args <- commandArgs(trailingOnly=TRUE)
	downsample_level <- args[1]

	#Load libraries.

	library(Seurat)

	#Read in gene-cell UMI count matrix.

	counts <- Read10X(paste0("CellRanger_output/",downsample_level,"/",downsample_level,"/outs/filtered_feature_bc_matrix"))

	#Create Seurat object.
	#Require genes to be expressed in 10+ cells.

	seurat.obj <- CreateSeuratObject(counts,min.cells=10)

	#Get mitochondrial rates.

	seurat.obj$percent.mito <- PercentageFeatureSet(seurat.obj,pattern="^mt-")

	#Remove cells with >30% mitochondrial rate.

	seurat.obj <- subset(seurat.obj,percent.mito < 30)

	#Remove cells in the bottom 5% by nGene.
	#This results in a cutoff of ~500 using the full dataset.

	nGene_cutoff <- quantile(seurat.obj$nFeature_RNA,probs=.05)

	seurat.obj <- subset(seurat.obj,nFeature_RNA > nGene_cutoff)

	#Removal of low-quality cells complete.
	#Ready to run standard pipeline up to PCA analysis.
	#Again we set seed as appropriate for reproducibility.

	seurat_seed = 1392

	seurat.obj <- NormalizeData(seurat.obj)
	seurat.obj <- FindVariableFeatures(seurat.obj,selection.method = "mean.var.plot")
	seurat.obj <- ScaleData(seurat.obj,features = VariableFeatures(seurat.obj))
	seurat.obj <- RunPCA(seurat.obj,features = VariableFeatures(seurat.obj),seed.use=seurat_seed)

	#Normally here, would explore the variance explained per PC to choose which PCs to run.
	#However this manual review is not very conducive to automated analysis.
	#So here, we run with the same number of PCs (first 15) for all downsample levels.
	#This should be more than sufficient given the relatively low number of cells here.

	PCs_to_use = 1:15

	seurat.obj <- FindNeighbors(seurat.obj,reduction = "pca",dims = PCs_to_use)
	seurat.obj <- FindClusters(seurat.obj,resolution = 0.6,random.seed=seurat_seed) 
	seurat.obj <- RunUMAP(seurat.obj,reduction = "pca",dims = PCs_to_use,seed.use=seurat_seed)

	#Save Seurat object.

	save(seurat.obj,file=paste0("seurat_object_",downsample_level,".Rdata"))

	#Output clustering results and QC metrics (nUMI and nGene and mito rates) per cell as well.

	write.csv(seurat.obj@meta.data,file=paste0("seurat_metadata_",downsample_level,".csv"),row.names=TRUE,quote=TRUE)

We were also curious how much clustering results would vary running Seurat multiple times on the same dataset. Therefore, we also ran the full data using the same pipeline above, but with different seeds (1393-1397).

These results were saved in files seurat_metadata_full_data_repeat_iteration_${1-5}.csv.

Now, let's begin comparing clustering results by using the adjusted rand index score to summarize how well the clusters correspond using a single number.

For the adjusted rand index score, we use a command from the fossil package.

Let's get the scores for all pairwise comparisons of 13%, 32%, 65%, and 100% downsample levels.

	library(fossil)

	downsample_levels <- c(13,32,65,100)

	adj_rand_index_scores <- matrix(NA,nrow=4,ncol=4,
		dimnames=list(as.character(downsample_levels),as.character(downsample_levels)))

	for(i in 1:4){adj_rand_index_scores[i,i] <- 1}

	for(larger_sample in c(32,65,100))
	{
		larger_sample_string <- ifelse(larger_sample == 100,
			"full_data",paste0("subsamp_",larger_sample))
		for(smaller_sample in downsample_levels[downsample_levels < larger_sample])
		{
				downsample <- read.csv(paste0("seurat_metadata_subsamp_",smaller_sample,".csv"),
					header=TRUE,row.names=1,stringsAsFactors=FALSE)
				full <- read.csv(paste0("seurat_metadata_",larger_sample_string,".csv"),
					header=TRUE,row.names=1,stringsAsFactors=FALSE)
				matched_cells <- intersect(rownames(downsample),rownames(full))
				downsample <- downsample[matched_cells,"seurat_clusters"]
				full <- full[matched_cells,"seurat_clusters"]
				adj_rand_index_scores[as.character(smaller_sample),as.character(larger_sample)] <- adj.rand.index(downsample,full)
		}
	}

Visualize using pheatmap.

	library(pheatmap)

	rownames(adj_rand_index_scores) <- paste0(downsample_levels,"%")
	colnames(adj_rand_index_scores) <- rownames(adj_rand_index_scores)

	pheatmap(adj_rand_index_scores,
		cluster_rows=FALSE,cluster_cols=FALSE,
		display_numbers=TRUE,
		main="Adj rand index scores between downsample levels")

We find that the clusters from the 13% downsample appear very different from not only the full data, but also from the 32% and 65% downsamples according to this metric.

We also find that the 32% and 65% downsamples appear to be relatively comparable in their clustering results according to this metric.

How do these adjusted rand index scores compare to different runs of Seurat on the exact same data?

	library(fossil)

	iteration0 <- read.csv("seurat_metadata_full_data.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
	
	for(i in 1:5)
	{
		other_iteration <- read.csv(paste0("seurat_metadata_full_data_repeat_iteration_",i,".csv"),header=TRUE,row.names=1,stringsAsFactors=FALSE)
		iteration0 <- iteration0[rownames(other_iteration),] #Cells will all match, just reorder as needed.
		print(round(adj.rand.index(iteration0$seurat_clusters,other_iteration$seurat_clusters),digits=2))
	}

Different runs of Seurat run on the exact same data, compared to an initial run, produce adjusted rand index scores from 0.91-0.96.

As we might expect, this is substantially higher than the 0.86-0.88 we see when comparing between the 32, 65%, and 100% downsamples.

However, it is notable that even different runs of the exact same data may be nowhere near identical.

Now, let's compare clustering results based on percentage of cells from each cluster assigned to new clusters for the full data vs. the 13, 32, and 65% downsamples, as well as the fourth iteration of running Seurat with a different seed between runs.



# References

Stoeckius, M., Zheng, S., Houck-Loomis, B. et al. Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome Biol 19, 224 (2018). https://doi.org/10.1186/s13059-018-1603-1

Stuart, T., Butler, A., Hoffman P., et al. Comprehensive Integration of Single-cell Data. Cell 177, 7 (2019). https://doi.org/10.1016/j.cell.2019.05.031
