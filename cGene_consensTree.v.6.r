# This script will be used to create a series of tree files from alignments of individual genes. Then we will generate consensus trees.  Further we will
# Create trees of the concatenations of 21 genes, and also review their consensus.  In this way we can review some of Rokas' 2003 work and comment on its
# accuracy.  NOTE: that consense will be run manually outside of this pipeline.

# After version 1's results from 28-05-2015 it was found that the trees had terribly low support values.  Thus we will modify this process
# by looking only at phylogenetically informative non-synonymous sites for our tree creation.  We will no compare the consensus tree of each gene
# and a concatenation of the total.

# load library to read alignments
library(ape)
library(seqinr)
library(doMC)
library(foreach)

# First clean the session to make certain there are no other objects or items loaded
system("rm .RData")
# This also removes any objects loaded into R
rm(list=ls())

###################################################################################################
################# BEGIN OF PREDEFINITIONS #########################################################
###################################################################################################

######## Folder and working definitions ##########
# This is the basal location which we can shift depending on the server I'm using for work.
rootDir <- "/Users/Jdench/Desktop/"
baseDir <- paste(rootDir,"PARuns/Sequences/PA_wholeGenomes/",sep="")
# This is the directory into which we want to save information and operate with files
workDir <- paste(baseDir,"newTrials/",sep="")
# This is where we will find the sequences with which to work
seqDir <- paste(workDir,"rawAln/",sep="")
# This is the final directory into which final output would be placed
outputDir <- workDir
# This is the concatenated alignment files, first raw then cleaned
cleanAln_file <- paste(workDir,"wholeAln.fas",sep="")
# This is the name of the output folders based on if it is single gene, concatenated genes, goTerm genes, or synonymous SNPs
sGene_out <- paste(workDir,"perGene/",sep="")
sGene_file <- "allSng_trees.tree"
cGene_out <- paste(workDir,"cGene/",sep="")
cGene_file <- "allConcat_trees.tree"
cogGene_out <- paste(workDir,"cogGene/",sep="")
cogGene_file <- "cogTree.tree"
snpGene_out <- paste(workDir,"snpGene/",sep="")
snpGene_file <- "allSNP_trees.tree"
# This is the filename of the perl script for concatenating alignments, NOTE: it must be within cleanDir to work
cat2File <- "catfasta2phyml.pl"

########### Program file locations #################
# This is the location of fastTree and it's call name
fasttreeLoc <- paste(rootDir,"Jonathan_Dench/3rd_Party_Software/FastTree2.1.7/fasttree",sep="")
# This is the location of the program consense from PHYLIP and a control file we associate with it
consenseLoc <- paste(rootDir,"Jonathan_Dench/3rd_Party_Software/PHYLIP/phylip-3.695/exe/consense.app",sep="")
consenseCtl <- " "

raxmlLoc <- paste(rootDir,"/Jonathan_Dench/3rd_Party_Software/RAxML/raxmlHPC-PTHREADS-AVX",sep="")

######### Constants and Defaults ############
# This registers our cores for the foreach paralellisation
registerDoMC(cores=4)
# This is the number of genes to have in each concatenation, NOTE: we want at least 100 trees generated!
nGenes <- 40
# This defines whether or nor we are using a core genome alignment genes previously identified
coreAln <- FALSE

# This is the location and file name of the PA14 database table
load(paste(rootDir,"/PARuns/Sequences/PA14DBase/pa14_genicDB.RData",sep=""))

# This is the location and filename of a table of COG's from PA01 which we will use to cross reference the PA14 information
cogTable <- read.csv(file="/Users/Jdench/Desktop/PARuns/Sequences/PA01DBase/PA01_COG.csv",header=TRUE)
# Then we adjust the types for the columns to make handling easier, the calls which include substring due to the shape of the table being known
# to have a leading and trailing blank space character
for (x in c("org","pro")){ cogTable[,x] <- sapply(as.character(cogTable[,x]),USE.NAMES=FALSE,function(y){ as.numeric(substr(y,2,nchar(y)-1)) }) }
for (x in c("COG","cat")){ cogTable[,x] <- sapply(as.character(cogTable[,x]),USE.NAMES=FALSE,function(y){ substr(y,2,nchar(y)-1) }) }
cogTable$annotation <- as.character(cogTable$annotation)
# This is a list of the COG types that are of interest for us to pursue, these relate to information/conserved types of COG classes
cogList <- c("A","B","J","K","L")

###################################################################################################
################# END OF PREDEFINITIONS ###########################################################
###################################################################################################


##############################################################################################################
######################################### BEGIN FUNCTIONS ####################################################
##############################################################################################################

	
##############################################################################################################
########################################### END FUNCTIONS ####################################################
##############################################################################################################



###################################################################################################
########################### BEGIN OF BODY #########################################################
###################################################################################################


##### This portion is for the generaiton of the single and concatenated gene trees #######
# This will identify all of the cleaned seqeuence alignment files we've stored in our seqDir, this requires knolwedge of the naming format used.
# Here we use the conditional of coreAlign to see if we look at all the sequence files or only those defined by the coreGenes object
if (coreAln){
	# This is a vector of the files in seqDir folder which are part of the core genome as identified by Core_Genome_Assess.v.1.r
	load(paste(seqDir,"coreGenes.RData",sep=""))
	setFiles <- cleanFiles <- coreGenes
} else {
	# We identify all the alignment files that exist, modify the grep pattern to suit your naming style
	setFiles <<- cleanFiles <<- list.files(seqDir,pattern="_seqAln.fas")
}

# step 1 we create all the tree files for the single genes
sngGene_trees <- foreach(thisFile = cleanFiles, .combine="rbind") %dopar% {
	# Now we simply call fasttree on each of these files, we call it for a nucleotide alignment using gtr model
	system(paste(fasttreeLoc, " -nt -gtr -gamma < ", seqDir,thisFile," > ", sGene_out,substr(thisFile,1,nchar(thisFile)-4),".tree",sep=""))
	return(data.frame("File"=thisFile,"Complete"=TRUE))
}

# Now we need to create a single tree file with all the trees held inside so we can pas this to PHYLIP consense
treeFiles <- grep(".tree",list.files(sGene_out),value=TRUE)
# Now we will write a growing file with all of the trees identified by tree files
cat(readLines(paste(sGene_out,treeFiles[1],sep="")),file=paste(sGene_out,sGene_file,sep=""),append=FALSE,sep="\n")
for (x in 2:length(treeFiles)) { cat(readLines(paste(sGene_out,treeFiles[x],sep="")),file=paste(sGene_out,sGene_file,sep=""),append=TRUE,sep="\n") }

# Step 2 we will create catches of 21 files from the entire set, concatenate them and then pass this on to be run in fasttree

# This basically scrambles the order of the files so that when we draw sets they are random
setFiles <- setFiles[sample(1:length(setFiles),length(setFiles),replace=FALSE)]
geneBatches <- lapply(1:floor(length(setFiles)/nGenes),function(x){ 
	# This will return series' of 21 genes excepting on the last instance in which all remaining genes are returned
	return(setFiles[(1+(nGenes*(x-1))):if (x != floor(length(setFiles)/nGenes)){ (nGenes*x) } else { length(setFiles) } ]) 
	})
		
# Sanity check on geneBatches, where abc should then be equal to length(geneBatches) since this is the last entry
#abc <- sapply(geneBatches,function(x){ return(length(x)) })
#which(abc != nGenes)

cGene_trees <- foreach(thisBatch = 1:length(geneBatches), .combine="rbind") %dopar% {
	# Now we will load and concatenate all the genes listed in the index
	tmpConcat <- read.dna(paste(seqDir,geneBatches[[thisBatch]][1],sep=""),format="fasta")
	for (thisGene in 2:length(geneBatches[[thisBatch]])){ tmpConcat <- cbind(tmpConcat,read.dna(paste(seqDir,geneBatches[[thisBatch]][thisGene],sep=""),format="fasta")) }
	# Now we write out a temporary concatenation file to pass to fasttree with no gaps and each sequence on one line
	tmpOut_name <- paste(cGene_out,thisBatch,"_tmp.fas",sep="")
	write.dna(tmpConcat,tmpOut_name,format="fasta",nbcol=-1,colsep="")
	system(paste(fasttreeLoc, " -nt -gtr -gamma < ", tmpOut_name," > ",cGene_out,"cGenes_Batch_",thisBatch,".tree",sep=""))
	# Remove the tmp file
	unlink(tmpOut_name)
	# Return the data of this job's work
	return(data.frame("Batch"=thisBatch,"Complete"=TRUE,"Files"=list(geneBatches[[thisBatch]])))
}

# Now we need to create a single tree file with all the trees held inside so we can pas this to PHYLIP consense
treeFiles <- grep(".tree",list.files(cGene_out),value=TRUE)	
# Now we will write a growing file with all of the trees identified by tree files
cat(readLines(paste(cGene_out,treeFiles[1],sep="")),file=paste(cGene_out,cGene_file,sep=""),append=FALSE,sep="\n")
for (x in 2:length(treeFiles)) { cat(readLines(paste(cGene_out,treeFiles[x],sep="")),file=paste(cGene_out,cGene_file,sep=""),append=TRUE,sep="\n") }

# This will save the sngGene_trees and cGene_trees objects for later review
save(sngGene_trees,cGene_trees,file=paste(workDir,"geneTrees_report.RData",sep=""))
	

########### This portion is for creation of gene trees based off of information level COG term genes ##########
# Using our COG types and categories we will try and find all those Locus_Tag's related to COG's of interest

# If this has been done then we will have saved an RData object of our concatenated sequence, so we can more easily just load this 
# rather than bothering to rebuild it
if (file.exists(paste(cogGene_out,"cogGene_seqAln.fas",sep=""))){
	load(paste(cogGene_out,cogGene_file,sep=""))
} else {
	# We find all unique COG terms in the PA14 database
	pa14COG <- unique(as.character(pa14_genicDB$COG))
	# This will look at all the pa14COG's and return the Locus_Tag associated if that COG is of the group type in cogList
	useLocus_Tag <- unlist(sapply(pa14COG,USE.NAMES=FALSE,function(x){ 
		# we find x in our cogTable, and grab it's category (cat), NOTE this does not handle multiple returns, though that should not be possible
		tmpCat <- cogTable$cat[which(cogTable$COG == x)]
		# Now we evaluate if something was found
		if (length(tmpCat) > 0){
			# we break this into a character string and look for intersects with cogList, if one exists we return the pa14_genicDB$Locus_Tag
			if (length(intersect(strsplit(tmpCat,"")[[1]],cogList)) > 0) { 
				return( as.character(pa14_genicDB$Locus_Tag[which(as.character(pa14_genicDB$COG) == x)]) ) 
			} else { 
				return( NULL ) }
		} else {
			return( NULL )
		}
	}) )
	# We check if there were any null return values and remove them
	#if (length(which(useLocus_Tag == "RemoveMe")) > 0) { useLocus_Tag <- useLocus_Tag[-which(useLocus_Tag == "RemoveMe")] }
	
	#### This portion need not be run if the files from cleanSet have had genes removed which are not sufficiently core
	#### Otherwise this will review all files from useLocus_Tag in seqDir, find those with at least 70% of sequences having 70% of sequence length
	useLocus_Tag <- foreach(thisFile = useLocus_Tag, .combine="c") %dopar% {
		# We load each file, and review the alignment to see if most sequences have most of the sequence
		tmpAln <- read.dna(file=paste(seqDir,thisFile,"_seqAln.fas",sep=""),format="fasta",as.character=TRUE)
		# We will make a table of the characters in each sequence and review if "a","c","g","t" characters are 70% or more
		tmpReview <- sapply(1:nrow(tmpAln),function(x){ 
			tmpTable <- table(tolower(tmpAln[x,]))
			# Now we ask if the gap or N character elements are more than 70% of total
			return( (sum(tmpTable[unlist(sapply(c("a","c","g","t"),function(y){ which(dimnames(tmpTable)[[1]] == y) }))])/sum(tmpTable) >= 0.7) )
		})
		# Now if mean(tmpReview) < 0.7 we know that most sequences are just gappy space
		if (mean(tmpReview) < 0.7){
			return( NULL )
		} else {
			return( thisFile )
		}
	}
	########################
	########################
	
	# Using the clean Files previously identified we find those sequence aligments which exist in our cleaned set
	useFiles <- sapply(useLocus_Tag,USE.NAMES=FALSE,function(x){ if(length(which(grepl(x,cleanFiles))) > 0){ return(grep(x,cleanFiles,value=TRUE)) } else { return(NULL) } })
	# Now we re-order the set just for practical reasons of consistency in our concatenation
	useFiles <- useFiles[order(useFiles)]
	
	# Now we perform the concatenation of these genes
	for (thisFile in useFiles){
		system(paste("cp ", seqDir,thisFile," ",cogGene_out,sep="")) 
	}
	# Now we execute the perl script to concatenate the files
	system(paste("cd ",cogGene_out,"; ./",cat2File," -f *.fas > ",cogGene_out,"cogGene_seqAln.fas",sep=""))
	# Now we save some drive space and remove the moved gene files
	unlink(paste(cogGene_out,useFiles,sep=""))
}

# Now we create the tree file for this alignment if not already done
if (!file.exists(paste(cogGene_out,cogGene_file,sep=""))){
	# Call to fasttree
	system(paste(fasttreeLoc, " -nt -gtr -gamma < ", cogGene_out,"cogGene_seqAln.fas > ",cogGene_out,cogGene_file,sep=""))
}
# Now we make a call to RAxML to root our unrooted trees of interest
system(paste("cd ",outputDir," ; ",raxmlLoc," -m GTRCAT -f I -n testRooting -t ",cogGene_out,cogGene_file,sep=""))

# Now the tree information will be written within the outfile and must be manually extracted as I have not written anything to parse this or perform rooting.



# This was a quick and dirt script I created in order to extract some specific information about the genes used to create my cogTree
# This requires the useFiles object created by "cGene_consensTree.v.4.r" and requires that the genetic database "pa14_genicDB" be loaded.

# This creates some objects
cogTree_filesUsed <- useFiles
cogTree_genesUsed <- substr(cogTree_filesUsed,1,gregexpr("_", cogTree_filesUsed)[[1]][2]-1)
# This starts up the data-frame.
cogTree_refFrame <- data.frame("Locus_Tag"=cogTree_genesUsed,"File"=cogTree_filesUsed,"LocusName"=0,"File_ntLength"=0,stringsAsFactors=FALSE)
cogTree_refFrame$LocusName <- vapply(cogTree_refFrame$Locus_Tag, USE.NAMES = FALSE, FUN.VALUE=vector(mode="character",length=1),
								function(x){ as.character(pa14_genicDB$geneID[which(pa14_genicDB$Locus_Tag == x)]) })

# Now we simply look at the number of nucleotides in each of the files reported by this ref frame and extract that information
for (i in 1:nrow(cogTree_refFrame)){
	tmpFile <- read.dna(paste(cogGene_out,cogTree_refFrame$File[i],sep=""),format="fasta",as.character=TRUE)
	cogTree_refFrame$File_ntLength[i] <- ncol(tmpFile)
}

# As a sanity check I made certain that the size of the alignment created was the same as the sum of File_ntLength, I got a result of TRUE from the script below
### tmpFile <- read.dna("~/Desktop/PARuns/Sequences/PA_wholeGenomes/consensTree/cogGene/cogGene_seqAln.fas",format="fasta",as.character=TRUE)
### sum(cogTree_refFrame$File_ntLength) == ncol(tmpFile)
write.csv(cogTree_refFrame,file=paste(cogGene_out,"cogTree_refFrame.csv",sep=""),row.names=FALSE)