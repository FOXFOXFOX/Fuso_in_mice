# Fuso_in_mice
This repo contains data, metadata, analysis code and figures for the first two experiments in which I infected mice with Fusobacterium nucleatum 

####most analytic code is in fuso nmds.R file####

Summary update, 3/4/16:
I have analyzed my mouse fuso infection 16S results (more plots and code to come over the weekend) and have findings to discuss. The two experiments differed in the ability of fuso to colonize and the levels at which it colonized. Perhaps the most interesting and consistent finding is that fuso is more abundant in cecal and colon tissue than stool, by at least 1 order of magnitude. I need to finish the analysis over the weekend so there's not much to look at just yet in terms of plots but hang tight. 

### **Project next steps:**
Even though I love all of the information we can get from it, going forward it would ideal to have a faster/cheaper way than 16S to quantify colonization. now that we have "controls" we could design an assay
BUT:
- qPCR won't be able to get good specificity w primers, even in "negative" controls there are other species. this is consisted with Eric Marten's findings
	+ try taqman probes, we might already have primers, need to check w Nick 
	+ could consider barcoding, then primers are more specific. Fuso is genetically tractable. There is a transposon system but that introduces camR which other species are certainly resistant to and have same cat gene.
- similar issues with culturing
- can still try the broth subculturing experiments in antibiotic, compare output side by side, fuso WT and fuso stool
- GFP labeled fuso or antibiotic tagged with a transposon? This would be cool
	+ anaerobic GFP is expensive (evoglow). 
	+ talked with Eric and even when they do GFP fusions with constitutively active promoter, cant see GFP in gut or in stool or in plated stool. only in vitro.
	+ quantify in plate reader based on fluorescence, would be cool
	+ gets around background issue 
- FISH maybe best, quantitative FISH, sectioning of tissues. But this is no longer "quick and dirty"
- same with monoclonal antibodies
	
### **Several options going forward (while also troubleshooting other assays)**
- repeat this/these experiments as done before, 16S to see if this was just some infection error
- redo 8week gavage model, shorter even, try to prove that the other model itself even works
- look into cloning, GFP experiments
- try other antibiotics for mice to be on
- germ free mouse experiments? 
- other analytics I can/should get out of this data? 

##Analysis and experiment notes as of 3/3/16: 

OTU00110 appears to be the innoculum, BLASTs to an animalis strain, Emma's lab papers
show that the EAVG002 are closest to animalis. 

The V4 region of the EAVG-002 innoculum strain is:
ACGTATGTCACGAGCGTTATCCGGATTTATTGGGCGTAAAGCGCGTCTAGGTGGT
TATGTAAGTCTGATGTGAAAATGCAGGGCTCAACTCTGTATTGCGTTGGAAACTG
TATAACTAGAGTACTGGAGAGGTAAGCGGAACTACAAGTGTAGAGGTGAAATTCG
TAGATATTTGTAGGAATGCCGATGGGGAAGCCAGCTTACTGGACAGATACTGACG
CTAAAGCGCGAAAGCGTGGGTAGCAAACAGG 

The full 16S sequence is NCBI 4_8 or 21_1A from the EAV paper

The other OTUs are different subspecies of F. nucleatum and so we can treat those as 
background and exclude them from analysis. 
But wait- are these found in experiment 2 as well?
and why is there so much human Fuso in the mouse background?
Fuso at t = 0???

 
Analysis plans before meeting with Pat:
- make a stripchart of OTU110 rel abundance at 48 hr per location 
colon, uninfected/infected side by side for each experiment 
  + OTU110/12 to get % abundance 



