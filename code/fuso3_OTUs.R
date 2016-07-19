#code to make fuso OTU level plots by treatment
#7 9 16

library(RColorBrewer)

meta_file <- read.table(file='fuso3_metadata.txt', sep='\t', header = TRUE, fill = TRUE, row.names = 1)
shared_file <- read.table(file='real_fuso3.an.0.03.subsample.shared', header = T, sep= '\t', row.names=2)
tax_file <- read.table(file='real_fuso3.an.cons.taxonomy', header = T, sep ='\t', row.names=1)
meta_76 <- read.table(file='fuso3_metadata_minus4.txt', sep='\t', header=T, row.names=1, fill =T)

####Use Nick's awesome code to get the OTU abundance table sorted out right###

#Create df with relative abundances
rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))

#Create vector of OTUs with median abundances >1%
med_rel_abund_cage <- aggregate(rel_abund, by=list(meta_76$cage),median)
cage_IDs <- med_rel_abund_cage[,"Group.1"]
med_rel_abund_cage <- med_rel_abund_cage[,!colnames(med_rel_abund_cage) %in% 'Group.1']
OTUs_1 <- apply(med_rel_abund_cage, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]

###some df code may go here, skipping for now to see if I really need it...

#get taxonomy - df with columns for each level - 
# level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom)
## - convert taxonomy text list to df, remove percentages
## - subset df by desired tax level, then replace any unclassified with next level up
get_tax <- function(tax_level=1){
  if (tax_level %in% c(1:5)){
    taxonomy <- tax_file[OTU_list,]
    taxonomy <-  data.frame(do.call('rbind', strsplit(as.character(taxonomy$Taxonomy),';',fixed=TRUE)))
    taxonomy <- data.frame(sapply(taxonomy,gsub,pattern="\\(\\d*\\)",replacement=""))
    level <- 7-tax_level
    tax_out <- as.character(taxonomy[,level])
    for (i in level:2){
      next_level <- i-1
      tax_out[tax_out=='unclassified'] <- 
        as.character(taxonomy[tax_out=='unclassified',next_level])
    }
    return(data.frame(tax= tax_out, row.names=OTU_list))
  } else {print(
    'Error: Level must be 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum)')
  }
}

taxonomy_genus <- get_tax(1)
taxonomy_family <- get_tax(2)
taxonomy_phylum <- get_tax(5)

taxonomy_family <- na.omit(taxonomy_family)

sum_OTU_by_tax_level <- function(TAX_DF,OTU_DF){
  tax_levels <- as.character(unique(TAX_DF$tax))
  OUTPUT_DF <- data.frame(rep(0,length(rownames(OTU_DF))), row.names=rownames(OTU_DF))
  for (i in 1:length(tax_levels)){
    OTU_by_level <- rownames(TAX_DF)[TAX_DF$tax %in% tax_levels[i]]
    if (length(OTU_by_level)>1){
      level_column <- apply(OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]],1,sum)
    } else {
      level_column <- OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]]
    }     
    OUTPUT_DF[,i] <- level_column
    colnames(OUTPUT_DF)[i] <- tax_levels[i]
  }
  return(OUTPUT_DF)
}


#df of OTUs w abundances >1% for day 3 and day 9 
rel_abund_d3 <- rel_abund[meta_76$day == 3, OTUs_1]
rel_abund_d9 <- rel_abund[meta_76$day == 9, OTUs_1]

#get family breakdown by sample 
family_d3 <- sum_OTU_by_tax_level(taxonomy_family, rel_abund_d3)
family_d9 <- sum_OTU_by_tax_level(taxonomy_family, rel_abund_d9)


#makes plot of all cages for day 3 timepoint 
n=27
barplot(t(family_d3), col=rainbow(n), ylim=c(0,100), cex.lab=0.9, cex.axis=0.7, cex.names=0.7)

#barplot for all day 9 cages 
barplot(t(family_d9), col=rainbow(n), ylim=c(0,100), cex.lab=0.9, cex.axis=0.7, cex.names=0.7)

#Now let's combine the OTUs by cage/treatment and report median of each for day 3
#add a column of cages to the df

#color brewer for colors
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

cages <- meta_76$cage[meta_76$day==3]
cages <- na.omit(cages)
rel_abund_d3[47] <- cages
colnames(rel_abund_d3)[47] <- "cage"

rel_abund_d3_med <- aggregate(rel_abund_d3[,1:46], list(rel_abund_d3$cage), median)

d3_fam_med <- sum_OTU_by_tax_level(taxonomy_family, rel_abund_d3_med)
d3_tmt <- c("DSS", "CONV", "Vanc high", "Vanc low")
d3_fam_med[27] <- d3_tmt
#this code plots the d3 communities by treatment 
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(d3_fam_med[,1:26]), ylab='Relative Abundance', main="Microbial communities in mice on day 3", 
        col = getPalette(n), axis.lty=1, ylim=c(0,80), cex.axis=0.8, xlab="treatment group", names.arg=d3_fam_med$V27)
fam_labels3 <- as.character(unique(taxonomy_family$tax))
legend('topright', fill = getPalette(n), inset=c(-0.2,0), fam_labels3[1:26], cex=0.5)

#let's do this also for day 9 

cages9 <- meta_76$cage[meta_76$day==9]
cages9 <- na.omit(cages9)
rel_abund_d9[47] <- cages9
colnames(rel_abund_d9)[47] <- "cage"

rel_abund_d9_med <- aggregate(rel_abund_d9[,1:46], list(rel_abund_d9$cage), median)

d9_fam_med <- sum_OTU_by_tax_level(taxonomy_family, rel_abund_d9_med)
d9_tmt <- c("DSS", "Uninfected", "CONV", "Vanc high", "Vanc low")
d9_fam_med[27] <- d9_tmt

#this code plots the d10 communities by treatment 
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(d9_fam_med[,1:26]), ylab='% Relative Abundance', main="Microbial communities in mice on day 10", 
        col = getPalette(n), axis.lty=1, ylim=c(0,80), cex.axis=0.8, xlab="treatment group", names.arg=d9_fam_med$V27)
fam_labels3 <- as.character(unique(taxonomy_family$tax))
legend('topright', fill = getPalette(n), inset=c(-0.25,0), fam_labels3[1:26], cex=0.5)





