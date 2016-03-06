#code for analyzing fuso experiment 1 and 2 data

fuso_nmds <- read.table("fuso_mouse.final.tx.unique_list.thetayc.1.lt.ave.nmds.axes", header = T)

#grep out to get axes for each organ

cecum <- fuso_nmds[grep('cecum', fuso_nmds$group), c(2,3)]
cecal_contents <- fuso_nmds[grep('cecon', fuso_nmds$group), c(2,3)]
colon <- fuso_nmds[grep('colon', fuso_nmds$group), c(2,3)]
colon_contents <- fuso_nmds[grep('colcon', fuso_nmds$group), c(2,3)]
stool <- fuso_nmds[grep('poop', fuso_nmds$group), c(2,3)]
spleen <- fuso_nmds[grep('spleen', fuso_nmds$group), c(2,3)]

#take out data from Aranoff's samples. grep looks for any number, dash, 8, anthing else. invert
# changes to just give false results 
fuso_nmds<-fuso_nmds[grep("[0-9]+-8.+", fuso_nmds$group, invert=T),]

#build the NMDS graph

plot(fuso_nmds$axis1, fuso_nmds$axis2)

points(cecum, pch=16, col="blue")
points(cecal_contents, pch=16, col="red")
points(colon, pch=16, col="dark green")
points(colon_contents, pch=16, col="purple")
points(stool, pch=16, col="brown")
points(spleen, pch=16, col="pink")

legend <- c("cecum", "cecal_contents", "colon", "colon_contents", "stool", "spleen")
legend(x = "topright", legend, col = c("blue", "red", "dark green", "purple", "brown", "pink"), pch=16)

##maybe just do the below in python 
#start making metadata file and sorting things
#to get just stuff from experiment 2 

#this only gets samples ending in 2. need to get 96 samples somehow too 
expt2 <-fuso_nmds[grep("[a-z].2", fuso_nmds$group), c(1,2,3)]

#get 96 samples, then add to expt 2 list 
expt96 <- fuso_nmds[grep("[a-z]96", fuso_nmds$group), c(1,2,3)]
expt2data <- rbind(expt2, expt96)

#then do the same, do NOT on grep to get expt 1
expt1 <- fuso_nmds[grep("[a-z].2", fuso_nmds$group, invert = T), c(1,2,3)]

#the above ended up being kind of pointless but it was a good exercise

#now that I have a metadata file i can load it in with other stuff and make graphs
metadata <- read.delim("metadata.tsv", sep = "\t", header = T)

#merge metadata with axes file
full_data <- merge(fuso_nmds, metadata, by.x = 'group', by.y='sampleID')

subsampled <- read.table("fuso_mouse.final.an.unique_list.0.03.subsample.shared", header = T)

#merge metadata and subsampled file
fusoOTUs <- merge(full_data, subsampled)

#build taxonomy file
taxonomy_file <- read.table(file="fuso_mouse.final.an.unique_list.0.03.cons.taxonomy", header =T, row.names = 1)

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
subsampled <- subsampled[,!(colnames(subsampled) %in% c("label", "numOtus"))]

#don't need axis stuff for this section of analysis, just merge subsampled and metadata 
fusoOTUs <- merge(metadata, subsampled)

#now want to add all of the fuso columns together and just plot abundance. 
#so from taxonomy file we know the fuso species are Otu00110, Otu00243, Otu00313
# Otu00961, Otu00994, Otu01880

#pull out just fuso columns and sample names first. now i have a table of sample and fuso abundance
only_fuso <- subsampled[, (colnames(subsampled) %in% c("Group", "Otu00110", "Otu00243", "Otu00313", "Otu00961", "Otu00994", "Otu01880"))]

#can add columns together now. but maybe i need to make this rel abundance first? 

#added columns together to make a 7th sum column
only_fuso[,7] <- only_fuso[,2] + only_fuso[,3] + only_fuso[,4] + only_fuso[,5] + only_fuso[,6]

#rename sum column 
names(only_fuso)[7] <- "fuso_total"

#trim table to make easier for merging
fuso_avg <- only_fuso[, colnames(only_fuso) %in% c("Group", "fuso_total")]

fuso_abund <- merge(fuso_avg, metadata, by.x='Group', by.y='sampleID')

#load the metadata in separately by experiment 
expt1metadata <- read.delim("expt1_metadata.tsv", sep = "\t", header = T)
expt2metadata <- read.delim("expt2_metadata.tsv", sep = "\t", header = T)
names(expt2metadata)<- c("sampleID", "mouse", "cage", "day", "location", "experiment", "infected")

fuso_abund1 <- merge(fuso_avg, expt1metadata, by.x="Group", by.y="sampleID")
fuso_abund2 <- merge(fuso_avg, expt2metadata, by.x="Group", by.y="sampleID")

#subset for just infected for exp 1
fuso_1infected <- subset(fuso_abund1, infected == 'TRUE')
fuso_2infected <- subset(fuso_abund2, infected == 'TRUE')

#plot fuso abundance by location (only for infected)
plot(fuso_1infected$location, fuso_1infected$fuso_total, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium abundance by sample location, expt 1")
plot(fuso_2infected$location, fuso_2infected$fuso_total, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium abundance by sample location, expt 2")

#to write plot 1 to a PDF
pdf(file='~/Documents/Schloss_Lab/Data/Fuso_mouse/expt1abundancelocation.pdf', width=14, height=10)
plot(fuso_1infected$location, fuso_1infected$fuso_total, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium abundance by sample location, expt 1")
dev.off()

#to write plot 2 to a PDF
pdf(file='~/Documents/Schloss_Lab/Data/Fuso_mouse/expt2abundancelocation.pdf', width=14, height=10)
plot(fuso_2infected$location, fuso_2infected$fuso_total, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium abundance by sample location, expt 2")
dev.off()

#plot abundance for uninfected first, probs a good idea
fuso_1uninfected <- subset(fuso_abund1, infected == 'FALSE')
plot(fuso_1uninfected$location, fuso_1uninfected$fuso_total, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium abundance by sample location, uninfected, expt 1")

#only pull out OTU00110 as this appears to be the innoculum OTU
fuso_110 <- only_fuso[, colnames(only_fuso) %in% c("Group", "Otu00110")]
fuso_110_1 <- merge(fuso_110, metadata, by.x='Group', by.y='sampleID')

#expt 1 stuff
fuso_110_abund <- merge(fuso_110, expt1metadata, by.x="Group", by.y="sampleID")
fuso_100_infected <- subset(fuso_110_abund, infected == 'TRUE')

#expt 2 stuff 
fuso_110_abund2 <- merge(fuso_110, expt2metadata, by.x = "Group", by.y = "sampleID")
fuso_110_2infected <- subset(fuso_110_abund, infected == 'TRUE')

#this plots only the amount of OTU110 in the infected samples
plot(fuso_100_infected$location, fuso_100_infected$fuso_total, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium 110 abundance by sample location, infected mice, expt 1")

#uninfected
fuso_100_uninfected <- subset(fuso_110_abund, infected == 'FALSE')
plot(fuso_100_uninfected$location, fuso_100_uninfected$Otu00110, type = "p", xlab = "sample location", ylab = "n fusobacterium reads", main = "fusobacterium 110 abundance by sample location, uninfected expt 1")

#a strip chart makes it easier to tell whats going on 
#expt 1, looks like some neg controls were infected
stripchart(fuso_110_abund$Otu00110~fuso_110_abund$infected, vertical = 'TRUE', method = 'jitter', main = "n samples with fusobacterium OTU00110 present, expt 1", xlab = "infected?", ylab = "n OTU00110 fuso reads")

#expt 2, less problem with neg controls infected 
stripchart(fuso_110_abund2$Otu00110~fuso_110_abund2$infected, vertical = 'TRUE', method = 'jitter', main = "n samples with fusobacterium OTU00110 present, expt 2", xlab = "infected?", ylab = "n OTU00110 fuso reads")

# make stripchart for infected by location at 48hr for both experiments

#to get % relative abundance, divide Otu00110 column by 12 (because subsampled to 1200 reads per sample)
#expt 1 table
fuso_110_abund[,10] <- (fuso_110_abund[,2]/12)
names(fuso_110_abund)[10] <- "Otu00110_relAbund"

#xpt 2
fuso_110_abund2[,10] <- (fuso_110_abund2[,2]/12)
names(fuso_110_abund2)[10] <- "Otu00110_relAbund"

#expt 1 graphs
#fuso abundance in all samples, infected vs not.
stripchart(fuso_110_abund$Otu00110_relAbund~fuso_110_abund$infected, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 1", xlab = "infected?", ylab = "n OTU00110 fuso abundance")
#by location
stripchart(fuso_110_abund$Otu00110_relAbund~fuso_110_abund$location, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 1", xlab = "location", ylab = "n OTU00110 fuso abundance")
#by mouse
stripchart(fuso_110_abund$Otu00110_relAbund~fuso_110_abund$mouse, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 1", xlab = "mouse #", ylab = "n OTU00110 fuso abundance")
#by time
stripchart(fuso_110_abund$Otu00110_relAbund~fuso_110_abund$day, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 1", xlab = "day", ylab = "n OTU00110 fuso abundance")

#expt2 graphs 
#by infected or not
stripchart(fuso_110_abund2$Otu00110_relAbund~fuso_110_abund2$infected, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 2", xlab = "infected?", ylab = "n OTU00110 fuso abundance")
#by location
stripchart(fuso_110_abund2$Otu00110_relAbund~fuso_110_abund2$location, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 2", xlab = "location", ylab = "n OTU00110 fuso abundance")
#by mouse
stripchart(fuso_110_abund2$Otu00110_relAbund~fuso_110_abund2$mouse, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 2", xlab = "mouse", ylab = "n OTU00110 fuso abundance")
#by time
stripchart(fuso_110_abund2$Otu00110_relAbund~fuso_110_abund2$day, vertical = 'TRUE', method = 'jitter', main = "fusobacterium OTU00110 abundance, expt 2", xlab = "day", ylab = "n OTU00110 fuso abundance")

#to add colors feed each plot a list
#example for storage
stripchart(fuso_110_abund2$Otu00110_relAbund~fuso_110_abund2$day, vertical = 'TRUE', method = 'jitter', col = c("red", "orange", "yellow", "green", "blue", "purple"), main = "n samples with fusobacterium OTU00110 present, expt 1", xlab = "day", ylab = "n OTU00110 fuso abundance")


#To do:
#plot abundance by time point for infected expt 1 
#export as graphs, update Rmd file w graphs and thoughts and plan for Monday 
