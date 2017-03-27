# regulators_integrating_RNAseq_ChIPseq
This approach assumes that lineage regulators tend to show consistent changes in both RNA-seq and ChIP-seq data, e.g., increase in expression, increase in H3K4me3 and/or decreases in H3K27me3, while structural/non-regulatory genes tend to show very significant changes only in RNA-seq.  
There are three lineages HA (i.e., C-EC), HP (i.e., H-EC) and CVP (i.e., CPC).  
First, the method calculates an integrated score for each pairwise lineage comparison (e.g., HP vs. CVP). The integrated score is defined as: S<sub>RNAseq</sub> + S<sub>H3K4me3</sub> - S<sub>H3K27me3</sub>. S itself is defined as sign(log2 fold change) x(-log10(FDR)) in each data source. This approach therefore ranks genes that show increase in expression, increase in active H3K4me3 mark and decrease in H3K27me3 mark, i.e., consistent changes across 3 data sources, highest. 
