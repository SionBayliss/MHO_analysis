# Plot output of BLASR search (-m 4)

# Dependencies
library(ggplot2)
library(gridExtra)
library(cowplot)
library(plyr)

## Input data
args <- commandArgs(trailingOnly=TRUE)

# Output paths
figure_out=sprintf("%s.figure.pdf",args[5])
table_out=sprintf("%s.table.pdf",args[5])

# Group Names 
G1=args[6]
G2=args[7]

# BLASR results
blasr.Pass<-read.delim( args[1] , sep=" ", header=TRUE)
blasr.Fail<-read.delim( args[2] , sep=" ", header=TRUE)

# Original read lengths
lengths.Pass<-read.delim( args[3] , sep="\t", header=FALSE) 
lengths.Fail<-read.delim( args[4] , sep="\t", header=FALSE) 

## Generate some basic statistics

# Table summary for output
summary<-data.frame(Pass=rep(0,6), Fail=rep(0,6))
row.names(summary)<-c("# Reads","# BLASR Hits (% # Reads)", "Mean Hit Length (%)", "Mean Percentage Match (%)", "# Hits < 75% Length (%)",  "# Hits >= 75% Length (%)"  )

# No unique reads in original fasta file
no_reads.Pass=length(lengths.Pass$V2)
no_reads.Fail=length(lengths.Fail$V2)
summary$Pass[1]<-no_reads.Pass
summary$Fail[1]<-no_reads.Fail

# BLASR hits against reference
hits.Pass=length(blasr.Pass$qName)
hits.Fail=length(blasr.Fail$qName)
per.unique.Pass=(hits.Pass/no_reads.Pass)*100
per.unique.Fail=(hits.Fail/no_reads.Fail)*100
summary$Pass[2]<-sprintf("%s (%.2f%%)",hits.Pass,per.unique.Pass)
summary$Fail[2]<-sprintf("%s (%.2f%%)",hits.Fail,per.unique.Fail)

# Find length of alignment as a percentage of original read length.
blasr.Pass$perl<-((blasr.Pass$qEnd-blasr.Pass$qStart)/blasr.Pass$qLength)*100
blasr.Fail$perl<-((blasr.Fail$qEnd-blasr.Fail$qStart)/blasr.Fail$qLength)*100

# Average hit length as percentage 
mean_hit.Pass<-mean(blasr.Pass$perl)
mean_hit.Fail<-mean(blasr.Fail$perl)
summary$Pass[3]<-sprintf("%.2f",mean_hit.Pass)
summary$Fail[3]<-sprintf("%.2f",mean_hit.Fail)

# Average hit percantage 
mean_per.Pass<-mean(blasr.Pass$percentSimilarity)
mean_per.Fail<-mean(blasr.Fail$percentSimilarity)
summary$Pass[4]<-sprintf("%.2f",mean_per.Pass)
summary$Fail[4]<-sprintf("%.2f",mean_per.Fail)

# Number of reads with hit of 25-50% length original read
no_hit.lt75.Pass<-sum( (blasr.Pass$perl<75) )
no_hit.lt75.Fail<-sum( (blasr.Fail$perl<75) )
summary$Pass[5]<-sprintf("%s (%0.2f%%)",no_hit.lt75.Pass,(no_hit.lt75.Pass/hits.Pass)*100)
summary$Fail[5]<-sprintf("%s (%0.2f%%)",no_hit.lt75.Fail,(no_hit.lt75.Fail/hits.Fail)*100)

# Number of reads with hit of 25-50% length original read
no_hit.gt75.Pass<-sum( (blasr.Pass$perl>=75) )
no_hit.gt75.Fail<-sum( (blasr.Fail$perl>=75) )
summary$Pass[6]<-sprintf("%s (%0.2f%%)",no_hit.gt75.Pass,(no_hit.gt75.Pass/hits.Pass)*100)
summary$Fail[6]<-sprintf("%s (%0.2f%%)",no_hit.gt75.Fail,(no_hit.gt75.Fail/hits.Fail)*100)

# Rename table columns for plotting
summary<-rename(summary, c("Pass"=G1, "Fail"=G2))

# Grob for plot
table<-tableGrob(summary)


# Plots 

# Set common theme 
theme_set(theme_gray())
default_theme<-theme(
      text = element_text(size=16, face="bold"), 
      plot.title = element_text(size=18, face="bold"),
      axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
      legend.text=element_text(size=12),legend.title=element_text(size=14), 
      legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border=element_rect(colour = "black",fill=NA), panel.background = element_blank(),
      legend.position=c(0.13,0.93) )
colours_2=c("#d95f02", "#1b9e77")

# Old 4 plot 
# default_theme<-theme(
#   text = element_text(size=20, face="bold"), 
#   plot.title = element_text(size=22, face="bold"),
#   axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
#   legend.text=element_text(size=12),legend.title=element_text(size=16), 
#   legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.border=element_rect(colour = "black",fill=NA), panel.background = element_blank(),
#   legend.position=c(0.15,0.90) )


# Prepare blasr results
df<-rbind(blasr.Pass,blasr.Fail)
df$status<-c(rep(G1, hits.Pass),rep(G2, hits.Fail))

# Prepare length data
lengths.Pass$status<-G1
lengths.Fail$status<-G2
df.lengths<-rbind(lengths.Pass, lengths.Fail) 

# Density plot of Pass and Fail read length distribution
l_hist<-ggplot(df.lengths, aes(V2, color=factor(status, levels = c(G1,G2)), fill=factor(status, levels = c(G1,G2)) )) + 
  geom_histogram(breaks=seq(0,max(df.lengths$V2),500), position="identity", alpha=0.1)+
  ylab("# Reads")+ xlab("Read Length (bp)")+
  default_theme+
  scale_fill_manual(name="Read Group", values = colours_2 )+
  scale_color_manual(name="Read Group", values = colours_2 )

# Box and whisker plot of percentage match to reference
box<-ggplot(df, aes(status, percentSimilarity, fill=factor(status, levels = c(G1,G2))) )+geom_boxplot(alpha=0.1)+
  default_theme+theme(legend.position='none')+
  ylab("Percentage Similarity (%)") + xlab("Read Group")+
  ylim(50, 100)

# BLASR hit length distribution
h_hist<-ggplot(df, aes(perl, color=factor(status, levels = c(G1,G2)), fill=factor(status, levels = c(G1,G2)) )) + 
  geom_histogram(breaks=seq(0,100,1), position="identity", alpha=0.1)+
  ylab("# Reads")+ xlab("Match Length (% of Read Length)")+
  default_theme+
  scale_fill_manual(name="Read Group", values = colours_2 )+
  scale_color_manual(name="Read Group", values = colours_2 )

# Plot tiled output
#p<-plot_grid(l_hist, box, h_hist, table, labels = "AUTO", label_size = 24, rel_heights = 0.98)
#save_plot(figure_out, p, base_height=14, base_width=14) ## Old 4 plot
p<-plot_grid(l_hist, box, h_hist, ncol=3, labels = "AUTO", label_size = 24, rel_heights = 0.98)
save_plot(figure_out, p, base_height=7, base_width=18)
save_plot(table_out, table, base_height=3, base_width=5)
