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

#this takes out label and numOTUs column from shared file 
subsampled <- subsampled[,!(colnames(subsampled) %in% c("label", "numOtus"))]

#don't need axis stuff for this section of analysis, just merge subsampled and metadata 
fusoOTUs <- merge(metadata, subsampled)

