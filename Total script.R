###########################################
#Converting sff files to fastq files#
#required packages
BiocManager::install("R453Plus1Toolbox")
library(R453Plus1Toolbox)
setwd("~/Research/Renamed") #Change to wherever the renamed files are. 

#read the sff files into an object of SFFContainer
a <- readSFF("J1_R1.sff")
#to ensure that it reads as an SFFContainer object, you'll have to do each file individually. reading as a strong c(files) will only change it to a list, which will then not be able te be read by sff2fastq
##To change the sff to a fastq file
sff2fastq(a, "~/Research/fastq") #file name, outdirectory
#the outdirectory should already be made as a folder where all of these files go. if you want to rename the files (which I did not) You'll have to put the renamed file name after the outdir
###the readSFF and sff2fastq will have to be done individually for each file
#######just make sure they're all going to the same outdir.

#####################################
#The following tutorial is modified from:
#https://benjjneb.github.io/dada2/tutorial.html
#with edits by Carly D. Kenkel & Alizah Ali & Nicola Kriefall & Sarah Davies

# source("https://bioconductor.org/biocLite.R")
# biocLite("dada2")
# biocLite('vegan')
#####################################
library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

#Set path to unzipped, renamed fastq files
path <- "~/Dropbox/BU/BU_Teaching/EcolEvolGenomics/Metabarcoding/Community_Data/Comm_Data" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path)
#Let's make sure that all of our files are there
fns

################################
##### Trimming/Filtering #######
################################

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures reads are in same order
#fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files- these are old 454 data but most data are paired end

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; OTHERWISE MODIFY
sample.names <- sapply(strsplit(fastqs, ".fastq"), `[`, 1) #the last number will select the field for renaming
sample.names
# Specify the full path to the fnFs
fnFs <- file.path(path, fastqs)
fnFs

#########Visualize Raw data

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs[c(1,2,3,4,5,6,7,8,9)])
plotQualityProfile(fnFs[c(10,11,12,13,14,15,16,17,18)])
plotQualityProfile(fnFs[c(19,20,21,22,23,24,25,26,27)])
plotQualityProfile(fnFs[c(28,29,30,31,32,33,34,35)])
plotQualityProfile(fnFs[c(36,37,38,39,40,41)])

#Recommend trimming where quality profile crashes - in this case, forward reads mostly fine up to 300
#For common ITS amplicon strategies with paired end reads, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus. That is OK, just leave out truncLen. Make sure you removed the forward and reverse primers from both the forward and reverse reads though! 

#Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

# Filter
out <- filterAndTrim(fnFs, filtFs, truncLen= 300, #end of single end reads = approx. 300 bp
                     maxN=0, #DADA does not allow Ns
                     maxEE=1, #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=20, #N nucleotides to remove from the start of each read: ITS2 primer = F 20bp
                     rm.phix=TRUE, #remove reads matching phiX genome
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)
tail(out)

#A word on Expected Errors vs a blanket quality threshold
#Take a simple example: a read of length two with quality scores Q3 and Q40, corresponding to error probabilities P=0.5 and P=0.0001. The base with Q3 is much more likely to have an error than the base with Q40 (0.5/0.0001 = 5,000 times more likely), so we can ignore the Q40 base to a good approximation. Consider a large sample of reads with (Q3, Q40), then approximately half of them will have an error (because of the P=0.5 from the Q2 base). We express this by saying that the expected number of errors in a read with quality scores (Q3, Q40) is 0.5.
#As this example shows, low Q scores (high error probabilities) dominate expected errors, but this information is lost by averaging if low Qs appear in a read with mostly high Q scores. This explains why expected errors is a much better indicator of read accuracy than average Q.

################################
##### Learn Error Rates #######
################################
#DADA2 learns its error model from the data itself by alternating estimation of the error rates and the composition of the sample until they converge on a jointly consistent solution (this is similar to the E-M algorithm)
#As in many optimization problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
#Maximum cycles was set to 30, but Convergence was found after 4 rounds

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 

################################
##### Dereplicate reads #######
################################
#Dereplication combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

################################
##### Infer Sequence Variants #######
################################

#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S
setDadaOpt(BAND_SIZE=32)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaFs[[25]]

#construct sequence table
seqtab <- makeSequenceTable(dadaFs)
head(seqtab)

################################
##### Remove chimeras #######
################################
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# Identified 1 bimeras out of 117 input sequences.

sum(seqtab.nochim)/sum(seqtab)
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)
#For our sample, this ratio was 0.9998201, there was only 1 bimera

write.csv(seqtab,file="Alizah_seqtab.csv")
write.csv(seqtab.nochim,file="Alizah_nochim.csv")
################################
##### Track Read Stats #######
################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_AllData_final.csv",row.names=TRUE,quote=FALSE)


################################
##### Assign Taxonomy #######
################################

#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
#DADA2 provides a native implementation of the RDP's naive Bayesian classifier. The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.
#Here, I have supplied a modified version of the GeoSymbio ITS2 database listing more taxonomic info as phyloseq requires (Franklin et al. 2012)
#For example: GeoSymbio data (taken from "all clades" at https://sites.google.com/site/geosymbio/downloads):
#>A1.1
#modified version for phyloseq looks like this instead:
#>Symbiodinium; Clade A; A1.1

taxa <- assignTaxonomy(seqtab.nochim, "GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta", minBoot=5,multithread=TRUE,tryRC=TRUE,outputBootstraps=FALSE)
#minboot should be higher
#Obtain a csv file for the taxonomy so that it's easier to map the sequences for the heatmap.
write.csv(taxa, file="taxa.csv",row.name=TRUE,quote=FALSE)
unname(head(taxa, 30))
unname(taxa)

#Now, save outputs so can come back to the analysis stage at a later point if desired
saveRDS(seqtab.nochim, file="final_seqtab_nochim.rds")
saveRDS(taxa, file="final_taxa_blastCorrected.rds")

#If you need to read in previously saved datafiles
seqtab.nochim <- readRDS("final_seqtab_nochim.rds")
taxa <- readRDS("final_taxa_blastCorrected.rds")
head(taxa)

################################
##### handoff 2 phyloseq #######
################################

#import dataframe holding sample information
#have your samples in the same order as the seqtab file in the rows, variables as columns
samdf<-read.csv("variabletable.csv")
head(samdf)
head(seqtab.nochim)
head(taxa)
rownames(samdf) <- samdf$sample

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

#replace sequences with shorter names (correspondence table output below)
ids<-taxa_names(ps)
ids <- paste0("sq",seq(1, length(colnames(seqtab.nochim))))
colnames(seqtab.nochim) <- ids

#Bar-plots
top90 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:90]
ps.top90 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top90 <- prune_taxa(top90, ps.top90)

plot_bar(ps.top90, x="Sample", fill="Class") 

#visusalize via counts rather than abundances:
plot_bar(ps, x = "sample", fill= "Class") #+ facet_wrap("tank")
#
#Obtain a csv file for the phyloseq data. Will give you the abundances for each sample and class. Useful for constructing the heatmap. Also, enables you to use ggplot, and construct more aesthetically pleasing plot.
psz <- psmelt(ps.top90)
write.csv(psz, file="Phyloseqoutputfinal.csv")
p <- ggplot(psz, aes(x=Sample, y=Abundance, fill=Class))
p + geom_bar(stat="identity", colour="black")

################################
##### output 'OTU' table #######
################################

#seqtab.nochim is the 'OTU' table...but is a little unwieldy
#For Symbiodinium, sequence classification is not so great...
#want fasta file of 'OTUs' and table designated by 'OTU'

#First, output fasta file for 'OTUs'
path <-'~/Research/final.fasta'
uniquesToFasta(seqtab.nochim, path)

#then, rename output table and write it out
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
colnames(seqtab.nochim)<-ids

write.csv(seqtab.nochim,file="final_OutputDADA_AllOTUs_FocusYesOnly.csv",quote=F)

str(seqtab.nochim)

#subset data
# focus = subset_samples(ps, Focus= "Yes")
# seqtab<-otu_table(focus)
# ids <- paste0("sq", seq(1, length(colnames(seqtab))))
# colnames(seqtab)<-ids
# head(seqtab)
# write.csv(seqtab,file="final_OutputDADA_AllOTUs_FocusYesOnly.csv",quote=F)

#############################################
##Create the Diversity Plots#################
##############################################

library("Rmisc")

#Replace "measures" for whatever you want to use, for our purposes we used Shannon
df=data.frame(estimate_richness(ps, split=TRUE, measures = "Shannon"), sample_data(ps)$tank)
df
#Remove the host sample, as it is a single individual and does not give a valuable diversity measure values.
df_nohost=subset(df, sample_data(ps)$tank!= "host") 
df_nohost$tank=df_nohost$sample_data.ps..tank

diver=summarySE(df_nohost,measurevar="Shannon",groupvars=c("tank"), na.rm=T)

#We add the +0.001 as it allows the anova to work on it - doesn't work without for some reason. Also, considering our values +0.001 will make no substantial difference in the output.

lm2=aov((Shannon+0.001)~tank, data=df_nohost)

summary(lm2)

TukeyHSD(lm2)

ggplot(df_nohost, aes(x=tank, y=Shannon)) + 
  geom_violin(trim=FALSE,scale='width', aes(x=tank,y=Shannon, color=tank, fill=tank)) + 
  geom_boxplot(width=0.15) + 
  scale_fill_manual(values=c("#009ACD" , "blue2","#CD3700", "goldenrod1"))+ 
  scale_colour_manual(values=c("#191970", "blue4", "brown4", "darkgoldenrod4")) +
  geom_jitter(mapping = aes(x=tank, y=Shannon, colour=tank), df_nohost,width=0)+
  theme_bw()+
  xlab("Treatment") +
  ylab("Shannon Diversity")

#########################################################
#########################################################
###########Creating the Heatmap#########################
##########################################################
#########################################################

library(vegan)
library(MCMC.OTU)
setwd("/Users/daviessw/Desktop")
dat <- read.csv(file="rawcountswithseq.csv", sep=",", header=TRUE, row.names=1)
summary(dat)
str(dat)
colnames(dat)
dat$type
dat$B1=dat$sq1...B1+dat$sq11...B1+dat$sq12...B1+dat$sq13...B1+dat$sq16...B1+dat$sq17...B1+dat$sq18...B1+dat$sq20.B1+dat$sq23.B1+dat$sq25.B1+dat$sq29.B1+dat$sq37.B1+dat$sq44.B1
dat$B2=dat$sq2...B2+dat$sq5...B2+dat$sq8...B2+dat$sq10...B2+dat$sq14...B2+dat$sq21.B2+dat$sq28.B2+dat$sq30.B2+dat$sq33.B2+dat$sq35.B2+dat$sq43.B2
dat$B3=dat$sq3...B3+dat$sq4...B3+dat$sq6...B3+dat$sq7...B3+dat$sq9...B3+dat$sq19.B3+dat$sq22.B3+dat$sq31.B3+dat$sq32.B3+dat$sq36.B3
dat$A4a=dat$sq24.A4+dat$sq26.A4a
dat$A2=dat$sq34.A2    
dat$A3=dat$sq15...A3+dat$sq27.A3+dat$sq38.A3+dat$sq39.A3 
dat$A4.3=dat$sq40.A4.3+dat$sq41.A4.3+dat$sq42.A4.3 

colnames(dat)

head(dat)
goods2=dat[,48:54]
head(goods2)
# creating a log-transfromed normalized dataset for PCA:
nl=startedLog(data=goods2,count.columns=1:length(names(goods2)), logstart=1)

trans=as.data.frame(t(nl))
head(trans)
names(trans)=c( "FAV1", "FAV2", "FAV3", "FAV4", "FAV5", "FAV6", "DIP1", "DIP2",  "DIP3", "DIP4", "DIP5", "DIP6", "hostsed1", "hostsed2", "hostsed3", "hostsed4", "hostsed5",  "hostsed6", "hostsed7", "hostsed8", "hostsed9", "hostsed10", "hostsed11", "hostsed12", "hostsed13", "hostsed14",  "hostsed15", "hostsed16", "hostsed17", "hostsed18", "hostsed19", "hostsed20",  "hostsed21", "hostsed22", "host1", "sed1", "sed2", "sed3", "sed4", "sed5", "sed6")

library(pheatmap)
pheatmap(as.matrix(trans),cluster_cols=F,scale="row", show_rownames = F)

# Make annotation table for pheatmap
ann = data.frame(cond = c("FAV", "FAV", "FAV", "FAV", "FAV", "FAV", "DIP", "DIP",  "DIP", "DIP", "DIP", "DIP", "hostsed", "hostsed", "hostsed", "hostsed", "hostsed",  "hostsed", "hostsed", "hostsed", "hostsed", "hostsed", "hostsed", "hostsed", "hostsed", "hostsed",  "hostsed", "hostsed", "hostsed", "hostsed", "hostsed", "hostsed",  "hostsed", "hostsed", "host", "sed", "sed", "sed", "sed", "sed", "sed"))
rownames(ann) <- names(explc)

# Set colors
Var1        <- c("deepskyblue3",  "blue", "firebrick", "olivedrab3", "goldenrod2")
names(Var1) <- c("FAV", "DIP", "hostsed", "host", "sed")
anno_colors <- list(cond = Var1)

pheatmap(as.matrix(trans),annotation_col=ann,annotation_colors=anno_colors,cex=1.2,border_color=NA,clustering_distance_cols="correlation", show_rownames=T)

means=apply(trans,1,mean) # means of rows
explc=trans-means # subtracting them
head(explc)

pheatmap(trans,annotation_col=ann,annotation_colors=anno_colors,cex=1.2,border_color="grey",show_rownames=T, cluster_cols=F, color = colorRampPalette(c("white", "coral2"))(25))

head(nl)
head(dat)
check=subset(dat[,1:3])
head(check)
new=cbind(check, nl)
head(new)

new2=subset(new[c(13:34,36:41),])
plot(B1~tank, data=new2)
lm1=lm(B1~tank, data=new2)
summary(lm1)

#######################################################
##########################################################
#################Creating the Phylogenic Tree############
#########################################################
########################################################
#making phylogenetic tree using ggtree package
#April 26th, 2018
#Nicola Kriefall

getwd()
setwd("~/Desktop")

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")

#read in tree data in Newick format 
#this is the Newick output from Phylogeny.fr if you use the one-click tool from that site
                                
treeb <- read.tree("treenewick_cladeB.txt")
treea <- read.tree("~/Google Drive/Uptake Project/Tree/treenewick_cladeA.txt")

#color assignments:

# B19: "#4F5090" 
# B2: "#80C8F8" 
# B1: "#307068" 
# B10: "#98B8E8" 
# B3: "#31396F"

# A4a: #F057A9
# A4.3: #F8D0D8
# A4: #B83857
# A3: #870F10
# A2: #E06888
# A1.1: F9A8C0

#create the tree with custom colors & custom scale & custom labels
quartz()
ggtree(treeb)+geom_treescale(family="Times")+geom_tiplab(hjust=-0.05,color=c("black","#98B8E8","#307068","black","black","black","#4F5090","black","#80C8F8","black","#4F5090"),size=10,family="Times")+ggplot2::xlim(0, 0.04)+
  geom_text2(aes(subset = !isTip, label=label,hjust=1,vjust=1.3),color="red")+
  theme(text=element_text(family="Times"))

quartz()
ggtree(treea)+geom_treescale(family="Times")+geom_tiplab(hjust=-0.05,color=c("black","#F057A9","black","#F8D0D8","#B83857","black","black","#870F10","black","#F9A8C0","black","black"),size=10,family="Times")+ggplot2::xlim(0, 0.15)+
  geom_text2(aes(subset = !isTip, label=label,hjust=1,vjust=1.3),color="red")+
  theme(text=element_text(family="Times"))

########################################################
########################################################
################Uptake Plot#############################
#######################################################
######################################################

library(survival)
install.packages("survminer",dependencies=T)
install.packages("ggpubr")
library(survminer)
install.packages("cmprsk")
library(cmprsk)
library(ggplot2)
setwd("~/Dropbox/Documents/zooxUptake_2018/zooxUptake_response_jan2019/")
data <- as.data.frame(read.csv('uptakedata.csv'))
data
data$treat=factor(data$treat,levels=c("control","host","sed","hostsed"))

# creating individual factors
data$host=as.factor(as.numeric(data$treat %in% c("host","hostsed")))
data$sed=as.factor(as.numeric(data$treat %in% c("sed","hostsed")))

# adding survival
data$s <- Surv(time = data$day, event = data$uptake)
data$s1 <- Surv(time = data$day, event = rep(1,nrow(data))) # this one assumes that everyone took up zoox by day 35. 


##COX PROPORTIONAL HAZARD MODEL FOR SURVIVAL##
model0 <- coxph(s~sed*host, data) # problem because zero uptake in control
summary(model0)

model <- coxph(s1~sed*host, data)  # model assuming everyone took up zoox eventually. Day 35 or 1000, does not matter - somehow (I checked)
summary(model)
# sediment alone marginally increases uptake (2.5 -fold, p = 0.0544)
# host alone does not increase uptake significantly (1.1 fold, p = 0.82)
# interaction term is significant (p = 0.0218): 5.7 fold increase in uptake if sediment and host are combined


##significant!
##The hazard values gives you the "hazard" or the risk of that particular event happening. therefore a significant value would mean that there is a significant risk that the recruit would take up symbiodinium
##it's similar to odds ratios
##Also, it's better to compare the data either by the least amount of uptake/death or the most because of the way the model is. to ensure this, label the one you want to compare it by with "a" in the beginning of the name
####because R goes alphabetically. For our model, we had to do the most uptake comparision because of the zero uptake in the control which gave inconclusive results.


##CREATE OBJECTS FOR USE IN survfit FUNCTION##
tod <- data$day
symb <-data$uptake
treat <-data$treat

#####chartreuse3
##CREATE DATA FRAME CONTAINING ONLY VARIABLES OF INTEREST##
df <-data.frame(treat,tod,symb)
df


##SURVIVAL FITTED CURVE##
sfit <-survfit(Surv(tod, symb)~treat, data = df)


##PLOT IT##
##The survplot is very useful for mortality
ggsurvplot(sfit, position = "dodge") 

#TO MAKE THE CI PLOT
#KM or CI plots are useful to see the uptake or increase in abundance alongside time.
fit.CIF <- cuminc(ftime=tod,fstatus=symb,group=treat,cencode=0)
fit.CIF
#Gives you the values at that time, the anova value for the model (still significant!) and the variances at each time point

plot(fit.CIF,fun=function(x) 1-x, ylab="Proportion of Uptake",xlab="Days",lwd=3,lty=1, ljoin = 1, col = c("navyblue", "chartreuse3", "brown", "tan3"), font = 2, lend = 2)
