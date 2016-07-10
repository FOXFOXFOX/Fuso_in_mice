# 7 9 16 fuso experiment 3 analysis

metadata <- read.delim(file='fuso3_metadata.csv', sep=',', header = T)
metadata <- metadata[,1:7]
shared <- read.table(file='real_fuso3.an.0.03.subsample.shared', header = T)
fuso_OTUs <- merge(metadata, shared, by.x="sample", by.y="Group")
taxonomy_file <- read.table(file='real_fuso3.an.cons.taxonomy', header = T, row.names=1)

#take taxonomy column full of taxonomy names and assign it to taxonomy
taxonomy <- taxonomy_file$Taxonomy

#assign OTU row names to 'names'
names(taxonomy) <- rownames(taxonomy_file)

#use gsub as regex command to trim taxonomy OTUs to highest level of classification possible
#and remove punctuation

taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub(";unclassified", "", taxonomy)
taxonomy <- gsub("/.*", "", taxonomy)
taxonomy <- gsub(";$", "", taxonomy)
taxonomy <- gsub(".*;", "", taxonomy)

#change OTU names to taxonomy names. push this to file now?
otu <- gsub("Otu0*", "OTU~", names(taxonomy))
names(otu) <- names(taxonomy)

#this takes out label and numOTUs column from shared file 
shared <- shared[,!(colnames(shared) %in% c("label", "numOtus"))]

#fuso OTU numbers. need to figure out which is the right one
#153, 258, 259, 751, 929, 930, 1066
#lets assume what we want is OTU153. should use bin.seqs later to confirm

fuso153 <- shared[,(colnames(shared) %in% c("Group", "Otu00153"))]
fuso_153abund <- merge(fuso153, metadata, by.x="Group", by.y="sample")

#to get % relative abundance, divide Otu00110 column by 30.33 (because subsampled to 3033 reads per sample)
fuso_153abund[,9] <- (fuso_153abund[,2]/30.33)
names(fuso_153abund)[9] <- "Otu153_relAbund" 

#now we plot

stripchart(fuso_153abund$Otu153_relAbund~fuso_153abund$infected, vertical = 'TRUE', method = 'jitter', main = 'fusobacterium OTU153 rel abund', xlab = "infected?", ylab="relabund", ylim=c(0,1))

stripchart(fuso_153abund$Otu153_relAbund~fuso_153abund$location, vertical = 'TRUE', method = 'jitter', main = 'F. nucleatum relative abundance by sample site', xlab = "sample site", ylab="% Relative abundance", ylim=c(0,1))

stripchart(fuso_153abund$Otu153_relAbund~fuso_153abund$treatment, vertical = 'TRUE', method = 'jitter', main = 'F. nucleatum relative abundance in each treatment group', xlab = "treatment group", ylab="% Relative abundance", ylim=c(0,1))

stripchart(fuso_153abund$Otu153_relAbund~fuso_153abund$day, vertical = 'TRUE', method = 'jitter', main = 'fusobacterium OTU153 rel abund', xlab = "day", ylab="relabund", ylim=c(0,1))

stripchart(fuso_153abund$Otu153_relAbund~fuso_153abund$cage, vertical = 'TRUE', method = 'jitter', main = 'fusobacterium OTU153 rel abund', xlab = "day", ylab="relabund", ylim=c(0,1))





