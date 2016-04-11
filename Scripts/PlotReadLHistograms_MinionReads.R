# Plot individual and overlayed histograms for two nanopore long read fasta files (1 = Pass, 2 = Fail). 

# Dependencies
library(ggplot2)
library(cowplot)

# Read Length Distributions
# Script to generate read lengths in bash = awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' corrected_pass_S17.fasta > Length_combined_corrected.txt

# Files via command line
args <- commandArgs(trailingOnly = TRUE)

# File 1 - behind/top of bar charts
h_data1<-read.csv(args[1])
#h_data1<-read.csv("/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperAnalysis/CompleteFinal/Reads/minION/Read_length_fail_str14.txt")
colnames(h_data1)<-"Read_Length"

# File 2
h_data2<-read.csv(args[2])
#h_data2<-read.csv("/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperAnalysis/CompleteFinal/Reads/minION/Read_length_pass.txt")
colnames(h_data2)<-"Read_Length"

# Output Directory
output=args[3]
#output="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperAnalysis/CompleteFinal/Figures"

# Combined
combined=rbind(h_data1, h_data2)

# Summary
sprintf("Summary File 1")
summary(h_data1$Read_Length)
sprintf("Summary File 2")
summary(h_data2$Read_Length)
sprintf("Combined File 2")
summary(combined$Read_Length)

# Number of reads
no_reads1=length(h_data1$Read_Length)
sprintf("Number of Reads File 1 = %s",no_reads1)
no_reads2=length(h_data2$Read_Length)
sprintf("Number of Reads File 2 = %s",no_reads2)
co_reads=length(combined$Read_Length)
sprintf("Number of Reads Combined = %s",co_reads)

# ggplot Histogram of File 1
g_hist1<-ggplot(h_data1, aes(Read_Length))+geom_histogram(binwidth=500, fill="skyblue", colour="black")+
   ylab("Count")+
  ggtitle("Read Length Distribution of 2D Pass Nanopore Reads")+
  theme(text = element_text(size=20, face="bold"), plot.title = element_text(size=30, face="bold"),
      axis.text.x = element_text(size=16),
       axis.text.y = element_text(size=16))
  
# Save as tiff
tiff(filename = sprintf("%s/Pass_ReadL_Hist.tiff", output), width = 1024, height = 960, units = "px", pointsize = 12)
g_hist1
dev.off()


# ggplot Histogram of File 2
g_hist2<-ggplot(h_data2, aes(Read_Length))+geom_histogram(binwidth=500, fill="firebrick3", colour="black")+
  ylab("Count")+
  ggtitle("Read Length Distribution of 2D Fail Nanopore Reads")+
  theme(text = element_text(size=20, face="bold"), plot.title = element_text(size=30, face="bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
g_hist2

# Save as figure.
tiff(filename = sprintf("%s/Fail_ReadL_Hist.tiff", output), width = 1024, height = 960, units = "px", pointsize = 12)
g_hist2
dev.off()


# Combined figure
c_data<-data.frame(Read_Length=rbind(h_data1, h_data2), FileID=c(rep("Pass",length(h_data1$Read_Length)), rep("Fail",length(h_data2$Read_Length))))
c_hist<-ggplot(c_data, aes(Read_Length, fill = factor(FileID, levels = c("Pass","Fail")))) + geom_histogram(binwidth=500, color="black") +
  ylab("Count")+
  ggtitle("Read Length Distribution of 2D Pass and Fail Nanopore Reads")+
  theme(text = element_text(size=20, face="bold"), plot.title = element_text(size=30, face="bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16), 
        legend.position=c(0.8,0.8))+
        #legend.title=element_text("Read Group"))+
  scale_fill_manual(name="Read_Group", values = c("skyblue", "firebrick3"))
c_hist


# Save as figure.
tiff(filename = sprintf("%s/ReadL_Hist.combined.tiff", output), width = 1024, height = 960, units = "px", pointsize = 12)
c_hist
dev.off()

