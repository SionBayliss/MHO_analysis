
##################
## Dependencies ##
##################
library(genoPlotR)
library(ggplot2)
library(gridExtra)

################
## Parameters ##
################

# Arguements
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=17){ die("Not enough input arguements") }

# ### TEMP
# args=NULL
# args[1]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Alignments/ChromosomeVsUSA300.backbone"
# args[2]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Annotation/Chromosome/Chromosome.gbk"
# args[3]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Annotation/USA300_FPR3757/USA300_FPR3757.gbk"
# args[4]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Scripts/Figure_Annotation.tab"
# args[5]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Mapping/ChromosomeVsIllumina.coverage"
# args[6]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Mapping/ChromosomeVsNanopore.coverage"
# args[7]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Alignments/PlasmidAVsSAP046A.backbone"
# args[8]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Annotation/PlasmidA/PlasmidA.gbk"
# args[9]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Annotation/SAP046A/SAP046A.gbk"
# args[10]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Mapping/PlasmidAVsIllumina.coverage"
# args[11]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Mapping/PlasmidAVsNanopore.coverage"
# args[12]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Alignments/PlasmidBVsSAP046B.backbone"
# args[13]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Annotation/PlasmidB/PlasmidB.gbk"
# args[14]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Annotation/SAP046B/SAP046B.gbk"
# args[15]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Mapping/PlasmidBVsIllumina.coverage"
# args[16]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Mapping/PlasmidBVsNanopore.coverage"
# args[17]="/mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/PaperDrafts/GigascienceSubmission/Github/MHO_analysis/Figures/Figure2_Corrected.pdf"
#  
 

# ggplot theme details for all plots
g_theme<-theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               plot.title = element_text(lineheight=.8, face="bold"),
               axis.title.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.text.x=element_blank(),
               legend.position = "none",
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.text.y=element_text(size = 10),
               axis.title.y=element_blank(),
               axis.text.x=element_text(size = 10),
               plot.margin = unit(c(0,0,0,0), "npc"))

# Rotation for genoplotR annotation
rot<-"0"

# nanopore colouf fill
nan_col<-"#56B4E9"

# illumina colour fill
il_col<-"#D55E00"

################  
## Chromosome ##
################

## Backbone plot ##

# Data locations
chr_bbone<-args[1]
chr_gb1<-args[2]
chr_gb2<-args[3]
annotation_data<-read.table(args[4], sep="\t")

# Read data
bbone <- read_mauve_backbone(chr_bbone)
GB1 <- try(read_dna_seg_from_file(chr_gb1,tagsToParse=c("CDS","rRNA")))
GB2 <- try(read_dna_seg_from_file(chr_gb2,tagsToParse=c("CDS","rRNA")))

# Make 16s and 23s rRNA operon red. 
GB1$col[grep("S ribosomal RNA",GB1$product)] <- "red"
GB2$col[grep("S ribosomal RNA",GB2$product)] <- "red"

# Prepare Annotation
annot_1 <- annotation(rot=rep(rot,length(annotation_data$V5[annotation_data$V5=="USA300"])),x1=annotation_data[annotation_data$V5=="USA300",1], annotation_data[annotation_data$V5=="USA300",3], text=as.character(annotation_data[annotation_data$V5=="USA300",4]))
annot_2 <- annotation(rot=rep(rot,length(annotation_data$V5[annotation_data$V5=="MHO_001"])),x1=annotation_data[annotation_data$V5=="MHO_001",1], annotation_data[annotation_data$V5=="MHO_001",3], text=as.character(annotation_data[annotation_data$V5=="MHO_001",4]))

## Coverage Plots ##

# Illumina
chr_cov_data<-args[5]
chr_cov <- read.table(chr_cov_data)
#l <- as.integer(unlist(strsplit(as.character(chr_cov[1,1]), "_"))[4])
l <-length(chr_cov[,1]) # 
bin_width <-seq(0,l,l/1000)
m_cov_chr<-tapply(chr_cov[,3], cut(chr_cov[,2], bin_width), mean )

# Plot
df.chr<-data.frame(Coverage=m_cov_chr, Position=bin_width[1:1000])
p.chr<-ggplot(df.chr,aes(Position,Coverage))+ geom_area(stat = "identity", fill=il_col)+
  coord_cartesian(xlim = c(0, 2970000),ylim=c(0,200))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(labels=function(label) sprintf('%5.0f', label))+ # Pads y axis so axis margins are equal
  g_theme
#p.chr

# Nanopore
chr_nan_data<-args[6]
chr_cov_nan <- read.table(chr_nan_data)
#l_n <- as.integer(unlist(strsplit(as.character(chr_cov_nan[1,1]), "_"))[4])
l_n <- length(chr_cov_nan$V1) #
bin_width_n <-seq(0,l_n,l_n/1000)
m_cov_chr_nan<-tapply(chr_cov_nan[,3], cut(chr_cov_nan[,2], bin_width_n), mean )

# Plot
df.chr.nan<-data.frame(Coverage=m_cov_chr_nan, Position=bin_width_n[1:1000])
p.chr.nan<-ggplot(df.chr.nan,aes(Position,Coverage))+geom_area(stat = "identity", fill=nan_col)+
  coord_cartesian(xlim = c(0, 2970000),ylim=c(0,20))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(labels=function(label) sprintf('%5.0f', label))+ # Pads y axis so axis margins are equal
  g_theme
#p.chr.nan

##############
## PlasmidA ##
##############

## Backbone plot ##

# Data locations
p1_bbone<-args[7]
p1_gb1<-args[8]
p1_gb2<-args[9]

# Read data
bbone_p1 <- read_mauve_backbone(p1_bbone)
GB1_p1 <- try(read_dna_seg_from_file(p1_gb1,tagsToParse=c("CDS","rRNA")))
GB2_p1 <- try(read_dna_seg_from_file(p1_gb2,tagsToParse=c("CDS","rRNA")))

# Prepare Annotation
annot_1_p1 <- annotation(rot=rep(rot,length(annotation_data$V5[annotation_data$V5=="SAP046A"])),x1=annotation_data[annotation_data$V5=="SAP046A",1], annotation_data[annotation_data$V5=="SAP046A",3], text=as.character(annotation_data[annotation_data$V5=="SAP046A",4]))
annot_2_p1 <- annotation(rot=rep(rot,length(annotation_data$V5[annotation_data$V5=="PlasmidA"])),x1=annotation_data[annotation_data$V5=="PlasmidA",1], annotation_data[annotation_data$V5=="PlasmidA",3], text=as.character(annotation_data[annotation_data$V5=="PlasmidA",4]))

## Coverage Plots ##

# Illumina
p1_cov_data<-args[10]
p1_cov <- read.table(p1_cov_data)
#l_p1 <- as.integer(unlist(strsplit(as.character(p1_cov[1,1]), "_"))[4])
l_p1 <- length(p1_cov$V1)
bin_width_p1 <-seq(0,l_p1,l_p1/1000)
m_cov_p1<-tapply(p1_cov[,3], cut(p1_cov[,2], bin_width_p1), mean )
m_cov_p1[is.na(m_cov_p1)]<-0

# Plot
df.p1<-data.frame(Coverage=m_cov_p1, Position=bin_width_p1[1:1000])
p.p1<-ggplot(df.p1,aes(Position,Coverage))+ geom_area(stat = "identity", fill=il_col)+
  coord_cartesian(xlim = c(0, l_p1), ylim=c(0,350))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))+ # Pads y axis so axis margins are equal
  g_theme
#p.p1

# Nanopore
p1_nan_data<-args[11]
p1_cov_nan <- read.table(p1_nan_data)
#l_n_p1 <- as.integer(unlist(strsplit(as.character(p1_cov_nan[1,1]), "_"))[4])
l_n_p1 <- length(p1_cov_nan$V1)
bin_width_n_p1 <-seq(0,l_n_p1,l_n_p1/1000)
m_cov_p1_nan<-tapply(p1_cov_nan[,3], cut(p1_cov_nan[,2], bin_width_n_p1), mean)
m_cov_p1_nan[is.na(m_cov_p1_nan)]<-0

# Plot
df.p1.nan<-data.frame(Coverage=m_cov_p1_nan, Position=bin_width_n_p1[1:1000])
p.p1.nan<-ggplot(df.p1.nan,aes(Position,Coverage))+geom_area(stat = "identity", fill=nan_col)+
  coord_cartesian(xlim = c(0, l_p1),ylim=c(0,20))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))+ # Pads y axis so axis margins are equal
  g_theme
#p.p1.nan


##############
## PlasmidB ##
##############

## Backbone plot ##

# Data locations
p2_bbone<-args[12]
p2_gb1<-args[13]
p2_gb2<-args[14]

# Read data
bbone_p2 <- read_mauve_backbone(p2_bbone)
GB1_p2 <- try(read_dna_seg_from_file(p2_gb1,tagsToParse=c("CDS","rRNA")))
GB2_p2 <- try(read_dna_seg_from_file(p2_gb2,tagsToParse=c("CDS","rRNA")))

# Prepare Annotation
annot_1_p2 <- annotation(rot=rep(rot,length(annotation_data$V5[annotation_data$V5=="SAP046B"])),x1=annotation_data[annotation_data$V5=="SAP046B",1], annotation_data[annotation_data$V5=="SAP046B",3], text=as.character(annotation_data[annotation_data$V5=="SAP046B",4]))
annot_2_p2 <- annotation(rot=rep(rot,length(annotation_data$V5[annotation_data$V5=="PlasmidB"])),x1=annotation_data[annotation_data$V5=="PlasmidB",1], annotation_data[annotation_data$V5=="PlasmidB",3], text=as.character(annotation_data[annotation_data$V5=="PlasmidB",4]))

## Coverage Plots ##

# Illumina
p2_cov_data<-args[15]
p2_cov <- read.table(p2_cov_data)
#l_p2 <- as.integer(unlist(strsplit(as.character(p2_cov[1,1]), "_"))[4])
l_p2 <- length(p2_cov$V1)
bin_width_p2 <-seq(0,l_p2,l_p2/1000)
m_cov_p2<-tapply(p2_cov[,3], cut(p2_cov[,2], bin_width_p2), mean )
m_cov_p2[is.na(m_cov_p2)]<-0

# Plot
df.p2<-data.frame(Coverage=m_cov_p2, Position=bin_width_p2[1:1000])
p.p2<-ggplot(df.p2,aes(Position,Coverage))+ geom_area(stat = "identity", fill=il_col)+
  coord_cartesian(xlim = c(0, l_p2),ylim=c(0,8000))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))+ # Pads y axis so axis margins are equal
  g_theme
#p.p2

# Nanopore
p2_nan_data<-args[16]
p2_cov_nan <- read.table(p2_nan_data)
#l_n_p2 <- as.integer(unlist(strsplit(as.character(p2_cov_nan[1,1]), "_"))[4])
l_n_p2 <- length(p2_cov_nan$V1)
bin_width_n_p2 <-seq(0,l_n_p2,l_n_p2/1000)
m_cov_p2_nan<-tapply(p2_cov_nan[,3], cut(p2_cov_nan[,2], bin_width_n_p2), mean)
m_cov_p2_nan[is.na(m_cov_p2_nan)]<-0

# Plot
df.p2.nan<-data.frame(Coverage=m_cov_p2_nan, Position=bin_width_n_p2[1:1000])
p.p2.nan<-ggplot(df.p2.nan,aes(Position,Coverage))+geom_area(stat = "identity", fill=nan_col)+
  coord_cartesian(xlim = c(0, l_p2), ylim=c(0,20))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))+ # Pads y axis so axis margins are equal
  g_theme
#p.p2.nan


##############
## Plot all ##
##############

# Using viewports

#tiff(filename = args[17], width = 2500, height = 1000, units="px", res=125, compression="lzw")
#tiff(filename = args[17], width = 170, height = 68, units="mm", res=3175, compression="lzw")
pdf(file=args[17], width=19, height=6.8)

# Root vp #
root.vp <-viewport(width=1,height=1)
pushViewport(root.vp)
upViewport()

# Axis labels 

lab1.vp<-viewport(x=0.015, y=0.5, width=0.02, height=0.25, just=c("left","bottom"))
pushViewport(lab1.vp)
grid.text("Illumina",  x=.5, y=.5, gp=gpar(fontsize=12), rot=90)
upViewport()

lab2.vp<-viewport(x=0.015, y=0.75, width=0.02, height=0.25, just=c("left","bottom"))
pushViewport(lab2.vp)
grid.text("MinION",  x=.5, y=.5, gp=gpar(fontsize=12), rot=90)
upViewport()

cov.vp<-viewport(x=0.005, y=0.62, width=0.01, height=0.25, just=c("left","bottom"))
pushViewport(cov.vp)
grid.text("Read Coverage", x=.5, y=.5, gp=gpar(fontsize=14), rot=90)
upViewport()

cov.vp<-viewport(x=0.01, y=0.15, width=0.01, height=0.25, just=c("left","bottom"))
pushViewport(cov.vp)
grid.text("MAUVE Alignment", x=.5, y=.5, gp=gpar(fontsize=14), rot=90)
upViewport()



# Chromosome #

top_chr.vp<-viewport(x=0.032, y=0.74, width=0.712, height=0.24, just=c("left","bottom"))  
print(p.chr.nan,vp=top_chr.vp)

middle_chr.vp<-viewport(x=0.03, y=0.5, width=0.713, height=0.24, just=c("left","bottom"))
print(p.chr,vp=middle_chr.vp)

t1_chr.vp<-viewport(x=0.02, y=0.486, width=0.717, height=0.01, just=c("left","bottom"))
pushViewport(t1_chr.vp)
grid.text("MHO 001",  x=.5, y=.5, gp=gpar(fontsize=12))
upViewport()

bottom_chr.vp<-viewport(x=0.028, y=0.03, width=0.715, height=0.47, just=c("left","bottom"))
pushViewport(bottom_chr.vp)
#plot_gene_map(annotation_cex=0.4,dna_segs=bbone$dna_segs, comparisons=bbone$comparisons,annotations=list(annot_2, annot_1), plot_new=FALSE)
plot_gene_map(dna_segs=list(GB1, GB2), comparisons=bbone$comparisons,plot_new=FALSE,n_scale_ticks=8, annotation_height=1.5,
              annotation_cex=0.8,gene_type="side_blocks",dna_seg_scale=TRUE, scale=FALSE,annotations=list(annot_2, annot_1),scale_cex=0.6)
upViewport()

t2_chr.vp<-viewport(x=0.02, y=0.0075, width=0.727, height=0.03, just=c("left","bottom"))
pushViewport(t2_chr.vp)
grid.text("USA300 FPR3757", x=.5, y=.5, gp=gpar(fontsize=12))
upViewport()

# Plasmid A #

top_p1.vp<-viewport(x=0.7375, y=0.74, width=0.135, height=0.23, just=c("left","bottom"))
print(p.p1.nan,vp=top_p1.vp)

middle_p1.vp<-viewport(x=0.735, y=0.5, width=0.136, height=0.23, just=c("left","bottom"))
print(p.p1,vp=middle_p1.vp)

t1_p1.vp<-viewport(x=0.745, y=0.486, width=0.14, height=0.01, just=c("left","bottom"))
pushViewport(t1_p1.vp)
grid.text("Plasmid A",  x=.5, y=.5, gp=gpar(fontsize=12))
upViewport()

bottom_p1.vp<-viewport(x=0.745, y=0.03, width=0.136, height=0.47, just=c("left","bottom"))
pushViewport(bottom_p1.vp)
#plot_gene_map(dna_segs=bbone_p1$dna_segs, comparisons=bbone_p1$comparisons,plot_new=FALSE, # No gene annotation 
#              annotation_cex=0.4,gene_type="side_blocks",dna_seg_scale=TRUE, scale=FALSE,)
plot_gene_map(dna_segs=list(GB1_p1, GB2_p1), comparisons=bbone_p1$comparisons,plot_new=FALSE, annotations=list(annot_2_p1, annot_1_p1),
              annotation_cex=0.8,gene_type="side_blocks",dna_seg_scale=TRUE, scale=FALSE,scale_cex=0.6,annotation_height=1.5)
upViewport()

t2_p1.vp<-viewport(x=0.745, y=0.0075, width=0.14, height=0.03, just=c("left","bottom"))
pushViewport(t2_p1.vp)
grid.text("SAP046A", x=.5, y=.5, gp=gpar(fontsize=12))
upViewport()


# Plasmid B #

top_p2.vp<-viewport(x=0.882, y=0.74, width=0.107, height=0.23, just=c("left","bottom"))
print(p.p2.nan,vp=top_p2.vp)

middle_p2.vp<-viewport(x=0.8775, y=0.5, width=0.104, height=0.23, just=c("left","bottom"))
print(p.p2,vp=middle_p2.vp)

t1_p2.vp<-viewport(x=0.89, y=0.486, width=0.1, height=0.01, just=c("left","bottom"))
pushViewport(t1_p2.vp)
grid.text("Plasmid B",  x=.5, y=.5, gp=gpar(fontsize=12))
upViewport()

bottom_p2.vp<-viewport(x=0.89, y=0.03, width=0.1025, height=0.47, just=c("left","bottom"))
pushViewport(bottom_p2.vp)
#plot_gene_map(dna_segs=bbone_p2$dna_segs, comparisons=bbone_p2$comparisons,plot_new=FALSE, # No gene annotation 
#              annotation_cex=0.4,gene_type="side_blocks",dna_seg_scale=TRUE, scale=FALSE,)
plot_gene_map(dna_segs=list(GB1_p2, GB2_p2), comparisons=bbone_p2$comparisons,plot_new=FALSE, annotation_cex=0.8, annotation_height=1.5,
              gene_type="side_blocks",scale=FALSE,dna_seg_scale=TRUE,scale_cex=0.6,annotations=list(annot_2_p2, annot_1_p2))
upViewport()

t2_p2.vp<-viewport(x=0.89, y=0.0075, width=0.1, height=0.03, just=c("left","bottom"))
pushViewport(t2_p2.vp)
grid.text("SAP046B", x=.5, y=.5, gp=gpar(fontsize=12))
upViewport()

# Panel Labels
cov.vp<-viewport(x=0.02, y=0.842, width=0.01, height=0.25, just=c("left","bottom"))
pushViewport(cov.vp)
grid.text("A", x=.5, y=.5, gp=gpar(fontsize=18))
upViewport()

cov.vp<-viewport(x=0.7275, y=0.842, width=0.01, height=0.25, just=c("left","bottom"))
pushViewport(cov.vp)
grid.text("B", x=.5, y=.5, gp=gpar(fontsize=18))
upViewport()

cov.vp<-viewport(x=0.8725, y=0.842, width=0.01, height=0.25, just=c("left","bottom"))
pushViewport(cov.vp)
grid.text("C", x=.5, y=.5, gp=gpar(fontsize=18))
upViewport()


dev.off()