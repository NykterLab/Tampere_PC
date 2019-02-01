library(pheatmap)
library(DESeq2)

setwd("/bmt-data/genomics/projects/atac_workdir/")
load("R_data/RNA_differentiallyExpressed.RData")
load("R_data/RNAdds.RData")


m1 <- assay(rna_dds, "norm.counts")[rownames(de_PCvBPH),]
m2 <- assay(rna_dds, "norm.counts")[rownames(de_CRPCvPC),]

pheatmap(m1, show_rownames = FALSE,
	annotation_col = data.frame(row.names = colnames(m1),
		type = colData(rna_dds)[,"Type"],
		protocol = colData(rna_dds)[,"RNA.Isolation.method"]),
	annotation_row = data.frame(row.names = rownames(m1),
		lfc = de_PCvBPH$log2FoldChange,
		p.val = de_PCvBPH$padj,
		absdiff = de_PCvBPH$diff),
	main = sprintf("PC vs BPH - n = %s", nrow(de_PCvBPH)),
	scale = "row",
	filename = "/bmt-data/genomics/projects/atac_workdir/pictures/DEG/PCvBPH_row.pdf",
	height = 10,
	width = 10)

pheatmap(m2, show_rownames = FALSE,
	annotation_col = data.frame(row.names = colnames(m2),
		type = colData(rna_dds)[,"Type"],
		protocol = colData(rna_dds)[,"RNA.Isolation.method"]),
	annotation_row = data.frame(row.names = rownames(m2),
		lfc = de_CRPCvPC$log2FoldChange,
		p.val = de_CRPCvPC$padj,
		absdiff = de_CRPCvPC$diff),
	main = sprintf("CRPC vs PC - n = %s", nrow(de_CRPCvPC)),
	scale = "row",
	filename = "/bmt-data/genomics/projects/atac_workdir/pictures/DEG/CRPCvPC_row.pdf",
	height = 10,
	width = 10)


d1 <- cbind(log2FoldChange=de_PCvBPH$log2FoldChange, absolute_median_difference=de_PCvBPH$diff, p.value=de_PCvBPH$padj, m1)
d2 <- cbind(log2FoldChange=de_CRPCvPC$log2FoldChange, absolute_median_difference=de_CRPCvPC$diff, p.value=de_CRPCvPC$padj, m2)

write.csv(d1, file = "/bmt-data/genomics/projects/atac_workdir/tables/DEG/de_PCvBPH.csv")
write.csv(d2, file = "/bmt-data/genomics/projects/atac_workdir/tables/DEG/de_CRPCvPC.csv")