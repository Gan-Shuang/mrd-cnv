library(optparse)
library(karyoploteR)
library(biomaRt)
library(regioneR)
library(zoo)
##parse_setting
option_list <- list(make_option("--genes", type = "character", help = "Mark annotated genes."),
                    make_option("--input",type="character",help = "Path to CNV file Required."),
                    make_option("--ploidy",type="numeric",help = "Sample Ploidy."),
                    make_option("--outdir",type="character",help = "Output Dir."))
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
anno_genes=opt$genes
cnv_file=opt$input
ploidy=opt$ploidy
outputdir=opt$outdir
##gene position
setwd(outputdir)
#anno_genelist<-read.table(anno_genes,header = F,sep = "\t")
#colnames(anno_genelist)=c("gene")
#gene.symbols <- anno_genelist$gene
#ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#load('/mnt/13d1/ganshuang/containers/ensembl.RData')
#genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
#seqlevelsStyle(genes) <- "UCSC"
##CNV_load
data <- read.table(cnv_file,header = T,sep = "\t")
colnames(data)=c("chr","pos","end","copy_number","event","value","status","corrected_copy_number","corrected_call","corrected_copy_number")
data <- as.data.frame(data)
data$value=(data$value+log2(ploidy/2))/4+0.5
##plot_chr
for(i in c(1:22)){
  plot_data=data[which(data$chr==i),]
  png(paste0("chr",i,".png"), width = 8000, height =5000,res=800)
  kp<-plotKaryotype(genome="hg19",chromosomes=paste0("chr",i))
  kpAxis(kp,ymin=-2,ymax=2,numticks = 5)
  kpAddLabels(kp, labels="log2 Copy Number",srt=90,pos =1,data.panel = 1,label.margin=0.05)
  kpAddBaseNumbers(kp)
  kpPoints(kp, chr = paste0("chr",i), x=plot_data$pos, y=plot_data$value)
  kpAbline(kp, h=0.5, col="gray", r0=0, r1=1)
  kpAddCytobandLabels(kp, force.all=TRUE, srt=90, cex=0.6)
  #kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol,r0=0,r1=1.35)
  dev.off()
  print(paste0("plot chr",i))
}


