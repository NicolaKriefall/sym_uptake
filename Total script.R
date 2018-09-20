###########################################
#Converting sff files to fastq files#
#required packages
install.packages("R453Plus1Toolbox")
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

# installing packages: execute just once when first using a script:

source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
#OR
#
install.packages("~/Downloads/dada2-1.4", #to install from source, just indicate pkg download location
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))

#####################################

library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

#Set path to unzipped, renamed fastq files
path <- "~/Research/RenamedFastq" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path)
fns

################################
##### Trimming/Filtering #######
################################

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files



# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; OTHERWISE MODIFY
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #the last number will select the field for renaming
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)


#########Visualize Raw data

#First, lets look at quality profile of R1 reads
 
plotQualityProfile(fnFs[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)])
plotQualityProfile(fnFs[c(21,22,23,24,25,26,27,28,29)])


#Recommend trimming where quality profile crashes - in this case, forward reads mostly fine up to 300


#For common ITS amplicon strategies, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus. That is OK, just leave out truncLen. Make sure you removed the forward and reverse primers from both the forward and reverse reads though! 


# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter
out <- filterAndTrim(fnFs, filtFs, truncLen= 300, #leaves ~30bp overlap
              maxN=0, #DADA does not allow Ns
              maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
              truncQ=2, 
              trimLeft=c(20,21), #N nucleotides to remove from the start of each read: ITS2 primers = F 20bp; R 21bp 
              rm.phix=TRUE, #remove reads matching phiX genome
              matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
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



#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) #some issues with C2G and G2C variants being underestimated, but not terrible

################################
##### Dereplicate reads #######
################################
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
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

dadaFs <- dada(derepFs, err=errF, multithread=TRUE,HOMOPOLYMER_GAP_PENALTY=-1)

#now, look at teh dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]


################################
##### Merge paired reads #######
################################

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[2]])

summary((mergers[[2]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

######################################
##### Construct sequence table #######
######################################
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 294-304 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Do merged sequences all fall in the expected range for amplicons? ITS2 Pochon ~340bp-41bp primers; accept 294-304
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(200,371)] #again, being fairly conservative wrt length

table(nchar(getSequences(seqtab2)))
dim(seqtab2)

#doesn't change

################################
##### Remove chimeras #######
################################
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)
#For our sample, this ratio was 1, therefore there were no chimeras

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
#Here, I have supplied a modified version of the GeoSymbio ITS2 database (Franklin et al. 2012)
#
taxa <- assignTaxonomy(seqtab.nochim, "~/Research/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta", minBoot=5,multithread=TRUE,tryRC=TRUE,outputBootstraps=FALSE)
#Obtain a csv file for the taxonomy so that it's easier to map the sequences for the heatmap.
write.csv(taxa, file="Y.csv",row.name=TRUE,quote=FALSE)
unname(head(taxa, 30))
unname(taxa)

#Lowered bootstrap threshold from 50 to 5. Was not returning hits for many sequences. But reducing to 5 improved sequence return and identities largely match separate blastn search against the same database

#Now, save outputs so can come back to the analysis stage at a later point if desired
saveRDS(seqtab.nochim, file="final_seqtab_nochim.rds")
saveRDS(taxa, file="final_taxa_blastCorrected.rds")

#If you need to read in previously saved datafiles
seqtab.nochim <- readRDS("~/Research/final_seqtab_nochim.rds")
taxa <- readRDS("~/Research/final_taxa_blastCorrected.rds")




################################
##### handoff 2 phyloseq #######
################################

#import dataframe holding sample information
samdf<-read.csv("variabletablex.csv")
head(samdf)
rownames(samdf) <- samdf$Sam

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
ps.top90 <- prune_taxa(top305, ps.top305)

plot_bar(ps.top90, x="Sample", fill="Class") 

#visusalize via counts rather than abundances:
plot_bar(ps, x = "Sam", fill= "Class") #+ facet_wrap("tank")
#
#Obtain a csv file for the phyloseq data. Will give you the abundances for each sample and class. Useful for constructing the heatmap. Also, enables you to use ggplot, and construct more aesthetically pleasing plot.
psz <- psmelt(ps.top305)
write.csv(psz, file="Phyloseqoutputfinal.csv")
p <- ggplot(psz, aes(x=Sample, y=Abundance, fill=Class)) + colorblind_pal()
p + geom_bar(stat="identity", colour="black", aes(fill(values=c()))) #+ facet_wrap("tank")

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

write.csv(seqtab.nochim,file="final_OutputDADA_AllOTUs.csv",quote=F)

str(seqtab.nochim)
0.4951895
0.455
#subset data
focus = subset_samples(ps, Focus= "Yes")
seqtab<-otu_table(focus)
ids <- paste0("sq", seq(1, length(colnames(seqtab))))
colnames(seqtab)<-ids
head(seqtab)
write.csv(seqtab,file="final_OutputDADA_AllOTUs_FocusYesOnly.csv",quote=F)

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
getwd()
setwd("~/Desktop")
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")

library("ape")
biocLite("Biostrings")
library("ggplot2")
library("ggtree")

#read in tree data in Newick format 
#this is the output from Phylogeny.fr if you use the one-click tool from that site, or it is #the same format from Geneious. I used Geneious in order to do a custom alignment (GTR+I). The #default Muscle alignment used by Phylogeny.fr was not the most effective after using
#Jmodeltest to check the delta AIC values for the different alignment models.  

tree <- read.tree("finaltree2.txt")

#create the tree with custom colors & custom scale & custom labels
ggtree(tree,ladderize=TRUE,size=0.8)+geom_treescale(x=0.9,y=0)+geom_tiplab(hjust=-0.05,color=c("black","#98B8E8","#307068","black","#80C8F8","black","#31396F","#4F5090","black","#870F10","#F9A8C0","#F057A9","#B83857","#F8D0D8","black","#E06888"),size=5,family="Times")+ggplot2::xlim(0, 1.4)

########################################################
########################################################
################Uptake Plot#############################
#######################################################
######################################################

library(survival)
library(survminer)
library(cmprsk)
library(ggplot2)
setwd("~/Research")
data <- as.data.frame(read.csv('data.csv'))
data

##CREATE A SURVIVAL OBJECT, CHECK THAT IT WORKED##
s <- Surv(time = data$day, event = data$uptake)
class(s)


##COX PROPORTIONAL HAZARD MODEL FOR SURVIVAL##
model <- coxph(s~treat, data)
summary(model)
anova(model)
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



