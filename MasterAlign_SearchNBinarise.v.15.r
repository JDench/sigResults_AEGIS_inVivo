# This is a new file for the working of my PA14 analysis, this is to find the gene within the whole PA14 database and then create 
# an alignment with a target genome to scan thus finding the true start and stop positions within our alignment 

########## DEPENDENCIES ##########
# This is an updated version, and is dependent on cGene.aln.Clean.Index.v.1.r and cGene_aln.Blast.v.12.r

# First clean the session to make certain there are no other objects or items loaded
# This removes any hidden previous session info
system("rm .RData")
# This also removes any objects
rm(list=ls())
# This is to accomodate when we are writting files with particularly small branch lengths
# We are removing this for the moment as it is more practical that we write out small numbers in Sci notation
#options("scipen"=100, "digits"=5)

# Load our libraries
library(ape)
library(seqinr)
library(doMC)
library(foreach)

###################################################################################################
################# BEGIN OF PREDEFINITIONS #########################################################
###################################################################################################

###### System Folder settings ##########
# This is the basic directory within which we will find all other information
baseDir <- "/Users/Jdench/Desktop/"
# This is the main directory under which multiple sub directories may be found
mainDir <- paste(baseDir,"PARuns/Sequences/",sep="")
# Then we define the database directory where we find stored resources
dbDir <- paste(mainDir,"PA14DBase/",sep="")
# Then we define our working directory where changes will occur and be written
workDir <- paste(mainDir,"PA_FixedSite/newTrials/polyAln_0.8_retry/",sep="")
# This is the folder in which our blast program's bin folder is found
blastDir <- paste(baseDir,"Blast/ncbi-blast-2.2.30+/",sep="")
# This is the location and name of our blast databases against which we will identify gene locations
blastDB <- paste(blastDir,"db/cGene_Alignment/wholeAln_Kos",sep="")
# This secondary database is used to asses the internal nucleotide positions of sites based on the PA14 reference 
blastDB_intPos <- paste(dbDir,"blastDB_pa14/pa14",sep="")
# These are filesnames use for tmp files introduced in the polyFrame construction stage for int_ntPos identification
fastaOut <- "tmpCutout.fas"
blastOut <- "int_ntPos_blast.table"
# These are blast call setups I use to standardise the process and make changes easier.
outFmt <- ' -outfmt "6 qseqid sseqid length qstart qend sstart send evalue pident" '

######## File Information ##########
# This is the name of the data file which is a table of the gene, accessions and other information about the genomes we will work with.
pa14dbData <- "pa14_genicDB.RData"
# This identifies the tree file standard name if we are to perform searches on the tree
treeFile <- "wholeAln.tree"
# This is the location of all data tables that contain information for sites of interest we only need to extract the 
# Gene names when there is a SNP of interest, these should be RData objects which contain a foldMat object 
dataFiles <- c("PARuns/AnitaSites/Anita_foldMat.RData","PARuns/AWongSites/Wong_foldMat.RData")
# These are the names of alignment files and processed information used downstream.  Changing these should allow all downstream parts
# to handle modified information formats.
alnName	<- "polyAln"


######## Constants ##########     ###### MANY IMPORTANT CHANGES TO BE DONE IN THIS SECTION REVIEW CAREFULLY ######
# This setups how many cores we are calling for parallellisation
registerDoMC(cores=4)
# We initialise the geneType boolean
geneType <- 1
# This is a percent to identify the minimum number of differences that we will require to consider a site pattern
percDiffs <- 0.05
# This is the majority consensus threshold we are setting for our assignment of character states, typically 50%
consensVal <- 0.8
# These are the sequences names in our alignment that we define as outgroup or ancestral

# Here we are using information from a previous run, relating to phenotype information, all that must remain is the outGroup_names object
if (file.exists(paste(workDir,"susceptStrains.RData",sep=""))){
	print("Outgroups loaded from file susceptStrains.RData")
	load(paste(workDir,"susceptStrains.RData",sep=""))
	outGroup_names <- susceptStrains
	rm(susceptStrains)
} else {
	outGroup_names <- NULL
	print(paste("Outgroups are defined from the .r file values entered:", outGroup_names,sep=""))
}
	
# This is a size threshold for the number of site patterns before we fork
sizeThresh <- 20000		

# These are sequences that we know need to be removed from alignments before processing, they may be enterred here
# but could be updated when loading an index.RData file, depends on the method.
removeSeq <- NULL
# This is to define is we are manually entering the genes of interest, if NULL then we will load from dataFiles.
manGenes <- c("dnaA","dnaN","morA","gyrA","gyrB","nfxB","parC","parE","lpd3","ribD","rpoB","serC")

###################################################################################################
################# END OF PREDEFINITIONS ###########################################################
###################################################################################################


###################################################################################################
###################### BEGIN OF FUNCTIONS #########################################################
###################################################################################################

# This is to take long vectors and compress them into a string with a separator
mkString <- function(aVec,aSep){
	float <- NULL
	for (a in 1:length(aVec)){
		float <- paste(float,aVec[a],sep=aSep)
	}
	# If the separation character is not blank we trim the first character 
	if (aSep != ""){ float <- substr(float,2,nchar(float)) }
	return(float)
}

# This takes in a vector of indexes and returns the number and size of all continuous index series'
# This returns a matrix as the result where each row gives the index start,stop of blocks and length
# It takes as a second argument the step size for sequence
seqGroup <- function(aVec,stepSize){
	# initialise some objects
	tLength <- 1
	tRow <- 1
	tStart <- 1
	returnMat <- matrix(NA,nrow=length(aVec),ncol=3)
	colnames(returnMat) <- c("Start","Stop","Length")
	
	for (x in 2:length(aVec)){
		# We will search at every position now to see if the index is sequential (by =1)
		if (aVec[x] == (aVec[x-1] + stepSize)){
			# Then this is sequential
			tLength <- tLength + 1
		} else {
			# Otherwise we have found a discontinuity we will store the start, stop and length
			returnMat[tRow,] <- c(aVec[tStart],aVec[x-1],tLength)
			tLength <- 1
			tStart <- x
			tRow <- tRow + 1
		}
	}
	# We need to add the last sequential sequence found
	returnMat[tRow,] <- c(aVec[tStart],aVec[x-1],tLength)
	
	# We remove any of the remaining NA rows and return
	if (length(attr(na.omit(returnMat),which="na.action")) > 0){ returnMat <- returnMat[-attr(na.omit(returnMat),which="na.action"),] }
	return(returnMat)
}

# This will identify if all the objects in a vector fed in as v, are identical
all_identical <- function(v) all(sapply( as.list(v[-1]),FUN=function(z) {identical(z, v[1])}))

# This will identify if there at least a defined number of differences in the data series
# By differences we mean values that are not the dominant symbol in the vector
minDiffs <- function(inVec, nDiff){
	# We remove any instances other than 0,1 BECAUSE in this script we only call this for our traitMat
	# meaning that we are passing items that have been binarised where a value 2 is for gaps of unclear nucleotides
	inVec <- inVec[ which(inVec != 2) ] 
	# This takes in a vector and identifies if there at least n Differences
	tmpTable <- table(inVec)
	# If difference between the sum of characters and the amount of dominant character can be greater or equal to nDiff
	return( (sum(tmpTable) - max(tmpTable)) >= nDiff )  
}

# This is a function we send to foreach that will allow us to return lists and combine them
comb <- function(...) {
	# We create a list of lists out of the input given to us, this assumes the arguments input can be forced to a list.
	return( lapply(list(...),function(x){ as.list(x) }) )
}

###################################################################################################
###################### END OF FUNCTIONS ###########################################################
###################################################################################################


###################################################################################################
################# BEGIN OF OBJECT LOADING #########################################################
###################################################################################################
# We will load our alignment cleaning index, it holds the following key objects:
# removeSeq - this tells us which sequences can be removed before actualy data processing.
# wholeName - this is the filename which is used for our whole alignment file, suffixes not included
if (file.exists(paste(workDir,alnName,"_index.RData",sep=""))){
	load(paste(workDir,alnName,"_index.RData",sep=""))
} else {
	print(" No _index file loaded, may not be in error, but warning has been posted ")
}
# This loads a previously saved database from an R object (pa14_genicDB)
load(paste(dbDir,pa14dbData,sep=""))

# Then we check if we have tidied the alignment of removeBranhces type sequences 
if (file.exists(paste(workDir, alnName,"_reduced.RData",sep=""))){
	# We load and report
	print("loading the reduced alignment from previous saved data")
	load(paste(workDir, alnName,"_reduced.RData",sep=""))
} else {
	# In the event the alignment has never been saved as an RData object we read in the .fas and then write out an RData object for later use
	if (!file.exists(paste(workDir, alnName,".RData",sep=""))){
		# This file path is hard set, it sohuld be updated if being rerun on other data.
		assign(alnName,read.dna(paste(workDir,alnName,".fas",sep=""),as.character=TRUE,format="fasta"),pos=".GlobalEnv")
		save(list=(alnName), file=paste(workDir, alnName,".RData",sep=""))
	} else {
		# Otherwise we can simply load our file
		load(paste(workDir, alnName,".RData",sep=""))
	}
	# This removes the rows of sequences which pertain to our remove branches type files, IF they actual exist....
	if (length(intersect(removeSeq, rownames(eval(as.name(alnName))))) > 0){
		assign(alnName,eval(as.name(alnName))[- unlist(sapply(removeSeq,USE.NAMES=FALSE,function(x){ which(grepl(x,rownames(eval(as.name(alnName))))==TRUE) })) ,],pos=".GlobalEnv")
	}
	save(list=(alnName),file=paste(workDir, alnName,"_reduced.RData",sep=""))
}

# We evaluate if we are using dataFiles to load our geneVec, or if we are manually entering this information
if (length(manGenes) > 0){
	geneVec <- manGenes[order(manGenes)]
} else {
	# We load the foldMat information from our data sources and extract all unique gene names
	geneVec <- NULL
	for (x in dataFiles){
		# We load the RData which should include a foldMat object
		load(paste(baseDir,x,sep=""))
		# We scan the fold mat object 
		for (i in 1:nrow(foldMat)){
			# .. and extract all unique gene names found
			if (length(intersect(foldMat$Gene[i],geneVec)) == 0){
				geneVec <- c(geneVec, foldMat$Gene[i])
			}
		}
	}
}


# We load the tree file
tTree <- read.tree(paste(workDir,treeFile,sep=""))

###################################################################################################
################# END OF OBJECT LOADING ###########################################################
###################################################################################################

## This is a matrix which is large and of uncertain column uniqueness but it should be high
#playLarge <- matrix(rep(c(sample(c(0,1),500000,replace=TRUE)),15),nrow=100)
#playMega <- matrix(rep(c(sample(c(0,1),86900000,replace=TRUE)),15),nrow=237)

## This is an entirely unqiue large matrix
#playUnique <- matrix(c(sample(c(0,1),8500000,replace=TRUE)),nrow=100)

# This is a medium sized matrix to allow us to practice and improve the time requirements of the site patterns identification step
# there are 16 columns repeated 5000 times, col 1 is monomorphic, col {2-4},{5-7} identical, 8-16 unique types
#playMedium <- matrix(rep(c(rep(1,100),rep(c(1,0),150),rep(c(0,1),150),rep(c(1,0,0,1,1,0,1),100),rep(1,98),rep(0,100),1,1),5000),nrow=100)

# This is a very simple data set to inspect that the site pattern function is working
#playSmall <- matrix(c(rep(1,15),rep(0,15),rep(c(1,0,1),10),rep(c(0,1,0),10),rep(c(1,1,0,1),15),rep(1,14),0,rep(0,14),1,rep(1,15)),nrow=15)

# This is a little toy test matrix of letters for practice with the binarise function
# col1 - mono "a", col2 - equal("c","g"), col3 - dom "t", pos7-12 "minor", col4 - 1:6 gaps, 7-12 "t", col5 - no dominatn character, col 6 - "dom "g", 10-12 - "c"
# Repeat this series 10000 times
#playMat <- matrix(rep(c(rep("a",12),rep("c",6),rep("g",6),rep("t",6),rep(c("a","c"),3),rep("-",6),rep("t",6),rep("g",3),rep("c",3),rep("a",3),rep("t",3),rep("g",9),rep("c",3)),10000),nrow=12)

# This is a list of lists within which are multiple vectors
#abc <- list("abc"=list("a"=rep(1,10),"b"=rep(3,9),"c"=rep(7,50)),"def"=list("d"=rep(2,15),"e"=rep(7,11),"f"=rep(11,17)))

###################################################################################################
################# Begin Of Tree Handling ##########################################################
###################################################################################################

if (length(intersect(removeSeq, tTree$tip.label)) > 0){   
	# remove any multifurcations, tidy up very small branch lengths, remove node labels (beacuse Bayes traits can't handle them)
	tTree <- root(tTree,removeSeq)
	# We use the vector of indexes rather than the names themselves becuase it has been found that the naming 
	# shceme may cause matches to not be found where they exist when using strings
	tTree <- drop.tip(tTree, unlist(sapply(removeSeq,USE.NAMES=FALSE,function(x){ which(grepl(x,tTree$tip.label)==TRUE) })) )
}
tTree <- multi2di(tTree)
tTree$node.label <- NULL
tTree$edge.length[tTree$edge.length < 0.0001] <- 0.0001  

# Now we write out the tree file in basic and nexus format for downstream use:
write.tree(tTree, file=paste(workDir,sub(".tree","_reduced.tree",treeFile),sep=""))
write.nexus(tTree, file=paste(workDir,sub(".tree","_reduced.nexus",treeFile),sep=""))

# We re-load the tTree file to ensure that the corrOrder is properly synced with the tTree object written out
# and later fed to BayesTraits
tTree <- read.nexus(paste(workDir,sub(".tree","_reduced.nexus",treeFile),sep=""))

###################################################################################################
################# End Of Tree Handling ############################################################
###################################################################################################



###################################################################################################
################# BEGIN OF TRAIT MAT CREATION #####################################################
###################################################################################################

# This verifies if we have yet created a traitMat RData file, if not we should since this will be utilised by
# downstream files to avoid them having to load and perform this action themselves
if (!file.exists(paste(workDir,alnName,"_traitMat_5perc_",sub(".tree","Tree",treeFile),".RData",sep=""))){
	
	# This is a verification such that if some of our outGroup are within the removeSeq we will create a temp alignment
	# For this section while we binarise, this means that at the end we must then remove those rows from our binary alignment
	if (length(intersect(removeSeq, outGroup_names)) > 0){
		# We create a global variable called tmpRemove and assign the removeSeq which are in outGroup_names
		tmpRemove <- removeSeq[- sapply(intersect(removeSeq, outGroup_names),function(x){ which(removeSeq == x)}) ]
		# We create a global variable called tmpRemove and assign the removeSeq which are in outGroup_names
		assign("recallRemove",intersect(removeSeq, outGroup_names),  pos=".GlobalEnv")
		# We load the unreduced wholeAln and use the names from removeSeq which are not part of outGroup_names
		load(paste(workDir,alnName,".RData",sep=""))
		# We remove the tmpRemove rows of this alignment but also make certain tmpRemove is part of the alignment 
		# (due to posible pre-processing it not necessarly be)
		assign(alnName,eval(as.name(alnName))[setdiff(rownames(eval(as.name(alnName))),tmpRemove),], pos=".GlobalEnv")
	}
	
	# This will take the outGroup_names and then find their positions in the alignment for later use
	outGroup_rows <- unlist(sapply(outGroup_names,USE.NAMES=FALSE,function(x){ which(grepl(x,rownames(eval(as.name(alnName)))) == TRUE) }))
	
	# We will binarise our alignment object previously loaded then save it out
	# This is to track how long this process is taking just for interest
	startTime <- proc.time()
	# This assumes the alignment of sequences is going to be given as a matrix with nrow = number of sequences and ncol is length
	nSeq <- nrow(eval(as.name(alnName)))
	seqLen <- ncol(eval(as.name(alnName)))
	
	# Here we allow for a short changing of this step as it can be time intensive for long sequences
	if (file.exists(paste(workDir,alnName,"_traitMat.RData",sep=""))){
		print("Loading traitMat from saved data")
		load(paste(workDir,alnName,"_traitMat.RData",sep=""))
	} else {
		traitChunk <- list()
		# It has been found in practice that at 5.5 million parts foreach is not able to reassemble the pieces so we break this into smaller parts
		for (x in 1:100){
			# We establish which parts of the job that this section will handle
			thisChunk <- (1+ ceiling(seqLen/100)*(x-1)):(min(ceiling(seqLen/100)*x,seqLen))
			
			# For each column (site in the alignment), we will identify the majority consensus characters and assign them as out group nucleotides
			traitChunk[[x]] <- foreach(j = thisChunk, .combine="cbind") %dopar% {
				# This is a tracking mechanism which can be read in the out file
				if (j %% 100000 == 0){ print(paste("We have gotten to column number: ",j," of ",seqLen,sep="")) }
				# We count the numbers of characters present in this column of the alignment and set it as a proportion
				ntCounts <- table(eval(as.name(alnName))[outGroup_rows,j])
				ntCounts <- ntCounts/length(outGroup_rows)
				# Now we define the outgroup characters based on those which are present 
				outNT <- attr(ntCounts,"names")[which(ntCounts >= consensVal)]
				# Now we handle special cases where gap characters or degenerat nt code may exist
				if (length(intersect(outNT,"-")) > 0 || length(intersect(outNT,"N")) > 0) {
					# In this case if a gap is ancestral we assume it means we know nothing so all is ancestral, we incldue degenerate code letters
					outNT <- c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V")
				}
				# We are going to ignore cases of the degenerate two letters codes such as R, Y, S, M, etc...
				# the reason being if they are likely present at very low number, and are representative of sequencing error more often than not
				return( unlist(sapply(eval(as.name(alnName))[,j],USE.NAMES=FALSE,function(x){
					if (length(intersect(toupper(x),"-")) > 0 || length(intersect(toupper(x),"N")) > 0) { return(2) }
					# If the nt in inputAln[i,j] does match at least one of our outgroup nucleotides it is outgroup = 0
					else if (length(intersect(toupper(x),toupper(outNT))) > 0) { return(0) }
					# If the nt in inputAln[i,j] does not match our outgroup nucleotides this is derived = 1
					else if (length(intersect(toupper(x),toupper(outNT))) == 0) { return(1) }
					# In all other cases we treat it as uniformative and assign it the value of 2
					else { return(2) }
				}) ))	
			}
		}
		# Now we reassemble our traitChunk parts
		traitMat <- traitChunk[[1]]
		for (x in 2:100){ traitMat <- cbind(traitMat,traitChunk[[x]]) }
		
		# Now that we have built the traitMat having kept the outGroup sequences that intersect with those we are removing
		# we will simply re-assign the traitMat to be only the rows that are not elements of removeSeq
		if (length(intersect(removeSeq, outGroup_names)) > 0){
			# We use the removeSeq object and keep all those rows that exsit but are not part of the set from removeSeq
			recallKeep_rows <- unlist(sapply(setdiff(rownames(eval(as.name(alnName))),removeSeq),function(x){ which(rownames(eval(as.name(alnName))) == x) }))
			traitMat <- traitMat[recallKeep_rows,]
			# And we need to update our alignment object as being the reduced equivalent which means keeping only elements we recall to keep
			assign(alnName,eval(as.name(alnName))[recallKeep_rows,],pos=".GlobalEnv")
			# We recreate the outGroup_rows based on which outGroup_names still exist
			outGroup_rows <- unlist(sapply(setdiff(outGroup_names,removeSeq),USE.NAMES=FALSE,function(x){ which(grepl(x,rownames(eval(as.name(alnName)))) == TRUE) }))
			# We also fix the nSeq and seqLen objects
			nSeq <- nrow(eval(as.name(alnName)))
			seqLen <- ncol(eval(as.name(alnName)))
		}
		# And save this object for reloading 
		save(traitMat,outGroup_rows,list=(alnName),file=paste(workDir,alnName,"_traitMat.RData",sep=""))
	}
	# This is to record the time taken on this process, a matter of interest
	#print("foreach binariseDNA")
	print(proc.time() - startTime)
	# sanity checks to be printed into the out file 
	#print(str(eval(as.name(alnName))))
	#print(str(eval(as.name(alnName))))
	# We strip the colnames of traitMat as this is not required and may just cause confusion
	colnames(traitMat) <- NULL
    
    
	# This identifies the number of differences we must observe to be interested in a site pattern, to a minimum of 1
	numDiff <- max(1,round((nSeq * percDiffs),0))
	# If we have a previous allSites object we can simply load that
	if (file.exists(paste(workDir,"allSites_",alnName,"_5perc.RData",sep=""))){
		# We load those objects instead of remaking them
		print("Loading allSites_base from previous save")
		load(paste(workDir,"allSites_",alnName,"_5perc.RData",sep=""))
	} else {
		# This vector will store the positions of all site patterns found on first sweep
		allSites <- rep(NA,seqLen)
	    # These will track the total number of patterns and novel ones respectively
	    nbPatterns <- 1
	    # This is to track how long this process is taking just for interest
		startTime <- proc.time()
	    # We make a first pass and identify all those sites which have polymorphic patterns
	    for (j in 1:seqLen){
	    	# Reporting mechanism
	    	if (j %% 100000 == 0) { print(paste("1st Pass: We have gotten to column number ",j,sep="")) }
	    	# Verify that this site fits within our definition of polymorphic based on min number of differences
	    	if(minDiffs(traitMat[,j],numDiff)){ 
	    		# We store this location as being a site pattern, then increase our counter of number of patterns
	    		allSites[nbPatterns] <- j
	    		nbPatterns <- nbPatterns + 1
	    	}
	    }
	    # We adjust the number of patterns to remove that initialising "blank" required in our loop
	    nbPatterns <- nbPatterns - 1
	    # We identify the number of monomorphic sites as the difference between the sequence length and the number of patterns
	    numMono <- seqLen - nbPatterns
	    # We trim the allSites vector by removing NA's
	    allSites <- allSites[!is.na(allSites)]
	    print(paste("First pass took this long and we found ",length(allSites)," sites of interest.",sep=""))
	    print(proc.time() - startTime)
	    save(allSites, numMono, nbPatterns, file=paste(workDir,"allSites_",alnName,"_5perc.RData",sep=""))
	}
    
    # This is to track how long this process is taking just for interest
	startTime <- proc.time()
	
	# Since the identification or repeated site patterns is one of the longer processes herein we will break this into chunks should there be a critical number
	# of total site patterns
	if (length(allSites) >= 2*sizeThresh){
		# If we have a midpoint allLists object we may as well load it
		if (file.exists(paste(workDir,"allLists_",alnName,".mid.RData",sep=""))){
			print("Loading from previous allLists_",alnName,".mid.RData")
			load(paste(workDir,"allLists_",alnName,".mid.RData",sep=""))
		} else {
			# We need to create this initial allLists item
			
			# we identify chucnks that will be handled
			thisChunk <- list()
			for (x in (1: ceiling(length(allSites)/sizeThresh))){
				thisChunk[[x]] <- allSites[(1+((x-1)* sizeThresh)):min(length(allSites), sizeThresh*x)]
			}
	
			allLists <- foreach( workingChunk = 1:length(thisChunk), .combine="comb", .multicombine=TRUE) %dopar% {
				# We identify those sites we will look at so we can collapse our search as we progress
				tmpSearch_sites <- thisChunk[[workingChunk]]
				
				# This will be the column ID's of our sites of interest, since we feed this a cleaned alignment this will mean the ID's are not true genomic ones
				# This necessitates that we still blast our found sites of interest to recover where in the genome they truly lie.
				tmpSites_list <- list()
				# This matrix will be used to store those novel site patterns we found
				tmpMat <- matrix(nrow = nSeq, ncol = length(tmpSearch_sites))
				# This will keep trakc of which site pattern is being evaluated for matches
				thisSite <- 1
				# This is a safety check in the event that there were at least 2 site patterns that fit within out definition of polymorphic
				if (length(tmpSearch_sites) > 1){
					# Now we make a second pass and find all similar sites patterns, group and store them
					while (thisSite <= length(tmpSearch_sites) ){
						# This is a reporting mechanism
						if (thisSite %% 1000 == 0 || thisSite == 1) { 
							print(paste("We are now evaluating thisSite number ", tmpSearch_sites[thisSite],sep="")) 
							print(paste("There are this many sites that remain to be matched: ",length(tmpSearch_sites)-length(tmpSites_list),sep=""))
							print("So far this second pass has been running for this long")
							print(proc.time() - startTime)	
						}
						
						# We now transfer over the next stored site pattern to our pattern matrix
						tmpMat[,thisSite] <- traitMat[, tmpSearch_sites[thisSite]]
						# Then we evaluate to find all similar site patterns
						
						# If this is not our last site pattern in the allSites vector then we look for redundant patterns
						if (thisSite != length(tmpSearch_sites)) {
					    	tmpSites_list[[thisSite]] <- tmpSearch_sites[thisSite]
					    	for (j in tmpSearch_sites[(thisSite +1): length(tmpSearch_sites)]){
					    		if(sum(abs(tmpMat[,thisSite] - traitMat[,j]),na.rm=TRUE) == 0){
					    			# Now we want to remove the allSites indexes that are redundant positions, so not the first thisSite value
					    			tmpSearch_sites[which(tmpSearch_sites == j)] <- NA
					    			# We return the j index value to be added to the siteID list index of thisSite
					    			tmpSites_list[[thisSite]] <- c(tmpSites_list[[thisSite]],j)
					    		}
					    	}
						} else {
							# If we have gotten to the last entry and it has not been matched elsewhere then this is a unique site we add
							tmpSites_list[[thisSite]] <- tmpSearch_sites[thisSite]
							tmpSearch_sites[thisSite] <- NA
						}
						
						# Now we can remove all the NA values we might have introduced into allSites during the siteID build
						tmpSearch_sites <- tmpSearch_sites[!is.na(tmpSearch_sites)]
						# And we increase the value of thisSite before considering our next while loop
						thisSite <- thisSite + 1
					}
				} else if (length(tmpSearch_sites) == 1) {
					# This means that there was only 1 site with a pattern found so we store it but report this instance
					tmpMat[,thisSite] <- traitMat[, tmpSearch_sites[thisSite]]
					tmpSites_list[[thisSite]] <- tmpSearch_sites[thisSite]
					# We report this instance for posterity.
					print("Please note there is only 1 site pattern found and no replicates")
				} else {
					# This means there are zero or somehow fewer than zero site Patterns, thus an error
					print("There were no site patterns identified, review how site patterns are defined or alignment for more")
				} 
				# This returns the list of site patterns we have identified in thisChunk of work
				return(tmpSites_list)
			# This closes out our forearch loops, 
			}
			# Save a midpoint item in case
			save(allLists, allSites, file=paste(workDir,"allLists_",alnName,".mid.RData",sep=""))
		}
		
		# Now we will reduce the allLists object into a single list
		print("Now performing the combination of lists")

		# We will not sequentially reduce the allLists object by comparing pairs of lists and keeping only unqiue site patterns found among each.
		# Once all the lists are merged into 1 then the first object within the list allLists will no longer be a list, it should be numeric or character
		while (length(allLists) > 1 && class(allLists[[1]]) == "list"){
			print(paste("There are ",length(allLists)," lists that remain to be merged",sep=""))
			print(proc.time() - startTime)
			# This will keep updating our allLists item as the unique elements between two sequential lists
			allLists <- foreach( thisList = seq(1, length(allLists),by=2), .combine="comb", .multicombine=TRUE) %dopar% {
				# This will be the object we are building for return
				returnList <- list()
				
				# If thisList, thus the leading list, is the last list in allLists, it won't have a pair so return it entirely
				if (thisList == length(allLists)) {
					return( allLists[[thisList]] )
				} else {
					# Then we have a pair to merge
					for (thisList_item in 1:length(allLists[[thisList]])){
						# Each new item of the leading list is by default a unique item
						returnList[[length(returnList)+1]] <- allLists[[thisList]][[thisList_item]]
						# We start the otherList counter as thisList + 1 to check the next other list
						otherList <- thisList + 1
						# this will store any items that need to be removed due to matches found
						removeItems <- NULL
						# Now we compare every list item to the leading list's patternCounter based item
						for (otherList_item in 1:length(allLists[[otherList]])){
							if(sum(abs(traitMat[,allLists[[thisList]][[thisList_item]][1] ] - traitMat[,allLists[[otherList]][[otherList_item]][1]]),na.rm=TRUE) == 0){
				    			# This means these two objects are the same so we add this otherList_item to thisList_item
				    			returnList[[length(returnList)]] <- c(allLists[[thisList]][[thisList_item]],allLists[[otherList]][[otherList_item]])
				    			removeItems <- c(removeItems,otherList_item)
				    		}
						}
						# Now if we foubnd matches we can delete these list items to as they do not need to be added to our returnList
						if (length(removeItems) > 0 ){ 
							for (x in removeItems){ allLists[[otherList]][[x]] <- NULL }
						}
					}
					# Now we evaluate if otherList has become entirely empty
					if (length(allLists[[otherList]]) == 0) { 
						# if so then we've found all the matches that exist and we return our returnList
						return( returnList )
					} else { 
						# we need to add all remaining otherList items to the return list then return that completed list
						for (otherList_item in 1:length(allLists[[otherList]])){
							returnList[[length(returnList)+1]] <- allLists[[otherList]][[otherList_item]]
						}
						return ( returnList )
					}
				}
			}
		}	
		
	
		# Now we create our siteID from the allLists object which should be one list, if it is a list of lists this logical will return FALSE
		if (class(allLists[[1]]) != "list") {
			siteID <- allLists
			# This matrix will be used to store those novel site patterns we found
			patternMat <- matrix(nrow = nSeq, ncol = length(siteID))
			for (j in 1: length(siteID)) { patternMat[,j] <- traitMat[,siteID[[j]][1]] }
		} else {
			print(paste("Something has gone wrong in creating the allLists  (length = ",length(allLists),") and thus siteID object",sep=""))
		}
		
	# This is the alternative if allSites is not longer than 2* sizeThreshold	
	} else {
		# we can use our updated serial sapply method of site pattern analysis 
		startTime <- proc.time()
		# This will be the column ID's of our sites of interest, since we feed this a cleaned alignment this will mean the ID's are not true genomic ones
		# This necessitates that we still blast our found sites of interest to recover where in the genome they truly lie.
		siteID <- list()
		# This matrix will be used to store those novel site patterns we found
		patternMat <- matrix(nrow = nSeq, ncol = nbPatterns)
		# This will keep trakc of which site pattern is being evaluated for matches
		thisSite <- 1
		# This is a safety check in the event that there were at least 2 site patterns that fit within out definition of polymorphic
		if (length(allSites) > 1){
			# Now we make a second pass and find all similar sites patterns, group and store them
			while (thisSite <= length(allSites) ){
				# This is a reporting mechanism
				if (thisSite %% 1000 == 0 || thisSite == 1) { 
					print(paste("We are now evaluating thisSite number ",allSites[thisSite],sep="")) 
					print(paste("There are this many sites that remain to be matched: ",length(allSites)-length(siteID),sep=""))
					print("So far this second pass has been running for this long")
					print(proc.time() - startTime)	
					save(siteID, allSites, thisSite, patternMat, file=paste(workDir,"allSites_",alnName,".mid.RData",sep=""))
				}
				
				# We now transfer over the next stored site pattern to our pattern matrix
				patternMat[,thisSite] <- traitMat[,allSites[thisSite]]
				# Then we evaluate to find all similar site patterns
				
				# If this is not our last site pattern in the allSites vector then we look for redundant patterns
				if (thisSite != length(allSites)) {
			    	siteID[[thisSite]] <- allSites[thisSite]
			    	for (j in allSites[(thisSite +1): length(allSites)]){
			    		if(sum(abs(patternMat[,thisSite] - traitMat[,j]),na.rm=TRUE) == 0){
			    			# Now we want to remove the allSites indexes that are redundant positions, so not the first thisSite value
			    			allSites[which(allSites == j)] <- NA
			    			# We return the j index value to be added to the siteID list index of thisSite
			    			siteID[[thisSite]] <- c(siteID[[thisSite]],j)
			    		}
			    	}
				} else {
					# If we have gotten to the last entry and it has not been matched elsewhere then this is a unique site we add
					siteID[[thisSite]] <- allSites[thisSite]
					allSites[thisSite] <- NA
				}
				
				# Now we can remove all the NA values we might have introduced into allSites during the siteID build
				allSites <- allSites[!is.na(allSites)]
				# And we increase the value of thisSite before considering our next while loop
				thisSite <- thisSite + 1
			}
		} else if (length(allSites) == 1) {
			# This means that there was only 1 site with a pattern found so we store it but report this instance
			patternMat[,thisSite] <- traitMat[,allSites[thisSite]]
			siteID[[thisSite]] <- allSites[thisSite]
			# We report this instance for posterity.
			print("Please note there is only 1 site pattern found and no replicates")
		} else {
			# This means there are zero or somehow fewer than zero site Patterns, thus an error
			print("There were no site patterns identified, review how site patterns are defined or alignment for more")
		}    
	# This closes our our decision to split site patterns analysis by number of sites with patterns    
	}
	
	print("Now done creating the siteID")
	print(proc.time() - startTime)
	
	# Backup in case of failure as this would be the longest portion and unfortunate to lose.
	save.image(paste(workDir,"sitePat_",alnName,".mid.RData",sep=""))
	
   	# This is to prevent a crash out in the circumstances where all sites are monomorphic
	# It will cause us to simply check the next replicate
    if(length(siteID) >= 2){
	    # Now we convert our reduced matrix into one that Bayes traits can use, so we assign the gap
	    # characters back which BayesTraits simply uses as either 0 or 1, and we make all value characters.
	    # NOTE: This will also remove all the NA columns pre-generated in patternMat's initiation.
	    builtMat<- matrix(mapply(FUN=function(x){
	    	if(x == 2){
	        	return("-")
	        }else{
	        	return(as.character(x))
	        }
	    },patternMat[,1:length(siteID)]),nrow=nSeq,ncol=length(siteID))
	    ## Sanity check for debugging
	    #builtMat
    }  else {
    	print(paste("There were only ",length(siteID)," unique site patterns!",sep=""))
    }
    # This is to record the time taken on this process, a matter of interest
    print("It took this long to perform the second pass")
	print(proc.time() - startTime)
	
	# Now we will use the tree related to this alignment and find the columnar order of sequences by comparing the wholeAln names and the tip labels
	# This correct order will inform us how to reorder a site column before feeding it to BayesTraits in downstream applications since we want the column 
	# values to match up with the tip positions.
	corrOrder <- sapply(tTree$tip.label,function(x){ 
		# Because some names may start as numerics we need to handle how they will be returned with string quotes
		if (substr(x,1,1) == "'"){
			which(rownames(eval(as.name(alnName))) == substr(x,2,nchar(x)-1))
		} else {
			which(rownames(eval(as.name(alnName))) == x) 
		} })	
	
	# To verify the correlation between our coding method and our traitMat sites patterns we create a binary vector for each
	# sequence if it is, or not, an outGroup element
	rowTypes <- rep(1,nrow(eval(as.name(alnName))))
	rowTypes[outGroup_rows] <- 0
	# We want to review the correlation of our recoding method and the binarised traitMat site patterns generated 
	# we will ignore monomorphic positions since they are NULL to be considered as there will be by default 100% correlation
	# Our sites in traitMat will relate to all those indicated by siteID
	reviewCorrelation <- foreach(thisCol = traitMat[,unlist(siteID)], .combine="c") %dopar% {
		# We will return the correlation between recoding types and the site pattern, which while itself is somewhat
		# meaningless as many sites may not be correlated to our coding method, the recoding method that maximises,
		# this value will at least have the best representation of our binary re-coding and nucleotide values
		return( mean(rowTypes == thisCol) )
	}
	
	# We will return the correlation between our outGroups and recoding, we want a strong correlation
	# between our outGroup type and the traitMat, this means that our outGroup coding method is representative	
	outCorrelation <- foreach(thisCol = traitMat[outGroup_rows,unlist(siteID)], .combine="c") %dopar% {
		return( mean(rowTypes[outGroup_rows] == thisCol) )
	}

	# Now we save the objects we've created and unlink the intermediate files if they exist
	save(traitMat,patternMat,builtMat,numMono,siteID,outGroup_rows,corrOrder, reviewCorrelation, outCorrelation ,file=paste(workDir,alnName,"_traitMat_5perc_",sub(".tree","Tree",treeFile),".RData",sep=""))
	if (file.exists(paste(workDir,"allSites_",alnName,"_5perc.RData",sep=""))) { unlink(paste(workDir,"allSites_",alnName,"_5perc.RData",sep="")) }
	if (file.exists(paste(workDir,alnName,"_traitMat.RData",sep=""))) { unlink(paste(workDir,alnName,"_traitMat.RData",sep="")) }
	if (file.exists(paste(workDir,"allLists_",alnName,".mid.RData",sep=""))) { unlink(paste(workDir,"allLists_",alnName,".mid.RData",sep="")) }
	if (file.exists(paste(workDir,"allSites_",alnName,".mid.RData",sep=""))) { unlink( paste(workDir,"allSites_",alnName,".mid.RData",sep="")) }
	if (file.exists(paste(workDir,"sitePat_",alnName,".mid.RData",sep=""))) { unlink( paste(workDir,"sitePat_",alnName,".mid.RData",sep="")) }
} else {
	# otherwise we simply load the previously made traitMat information
	load(paste(workDir,alnName,"_traitMat_5perc_",sub(".tree","Tree",treeFile),".RData",sep=""))
}
###################################################################################################
###################### END OF TRAIT MAT CREATION ##################################################
###################################################################################################


###################################################################################################
###################### BEGIN OF REPORT and POLYMORPHISM CREATE ####################################
###################################################################################################
# We will now take advantage of the index files or colnames of our alignment in order to extract the location
# of geneVec items and the polymorphic sites within our alignments.  
##### NOTE:  Previously we had a robust method that used blast, and it functions perfectly BUT since we've 
#####        Updated our upstream mehtods to be more explicit and track better we do not need to follow a blast
#####        method.
##### NOTE:  This assumes the structure of the index object is data frame and contains certain colnames such as Gene

# This is a boolean to report if we have index file information or if we must default to colnames of our alnName object
indexFile <- FALSE

# We look for an index file related to our alnName, and if it exists we load it.
if (file.exists(paste(workDir,alnName,"_index.RData",sep=""))){
	load(paste(workDir,alnName,"_index.RData",sep=""))
	indexFile <- TRUE
} else {
	print(paste(" Could not find index file to load, this may affect downstream functionality if positions are not in colnames of ",alnName,sep=""))
}

# Now we will built an object that reports on the geneVec positions in our alnName object
##### NOTE: This assumes that the positions are linearly order and non-interupted! - Expected by upstream methods.
reportFrame <- foreach(thisGene = geneVec, .combine="rbind") %dopar% {
	# Our purpose is to return a data.frame row which reports on thisGene's positions in alnName
	tmpReturn <- data.frame("Gene"=thisGene,"Alignment"=alnName,"Index_Used"=indexFile,"Start"=0,"End"=0)
	if (indexFile){
		# This means we use the _index object that should have been loaded, this will indicate the positions for thisGene
		return( data.frame("Gene"=thisGene,"Alignment"=alnName,"Index_Used"=indexFile,
				"Start"= min(which(eval(as.name(paste(alnName,"_index",sep="")))$Gene == thisGene)),
				"End"= max(which(eval(as.name(paste(alnName,"_index",sep="")))$Gene == thisGene)), stringsAsFactors=FALSE) )
	} else {
		# This means we use the colnames of the alnName object on order to report on the positions 
		return( data.frame("Gene"=thisGene,"Alignment"=alnName,"Index_Used"=indexFile,
				"Start"= min(which(grepl(thisGene,colnames(eval(as.name(alnName)))))),
				"End"= max(which(grepl(thisGene,colnames(eval(as.name(alnName)))))),stringsAsFactors=FALSE ) )
	}
}


# Since alnName was used to create traitMat we have a correlation of nucleotide and binary format
# As a sanity check we will evaluate if the rownames of both are identical (NOTE: if we use an index file then
# the rownames may be stored as a different object as as.name(alnName)_names
matchingMats <- FALSE
if (indexFile){
	# We look for an object of as.name(alnName)_names if rownames of the alnName object is null
	if (length(rownames(eval(as.name(alnName)))) == 0){
		matchingMats <- all_identical(eval(as.name(paste(alnName,"_names",sep=""))) == rownames(traitMat))
	} else {
		matchingMats <- all_identical(rownames(eval(as.name(alnName))) == rownames(traitMat))
	}
} else {
	matchingMats <- all_identical(rownames(eval(as.name(alnName))) == rownames(traitMat))
}

# So if the traitMat and alnName objects do not match row per row we have a problem and don't proceed
if (matchingMats){
	# This stores the number of differences that will be required for a site to be polymorphic more than the percDiffs
	numDiff <- max(1,round((nrow(eval(as.name(alnName))) * percDiffs),0))
	# Now we create a data frame that houses the nucleotide and binary site pattern for each polymorphic site within the span of thisGene
	polyFrame <- foreach(thisGene = geneVec, .combine="rbind") %dopar% {
		# We establish a return data frame object 
		returnFrame <- data.frame("Gene"= thisGene,"int_ntPos"=NA, "nt_Seq"=NA, "binary_Seq"=NA, stringsAsFactors=FALSE)
		# Now we cycle through each position bounded by tmpStart and tmpEnd to look for polymorphic sites
		# in the collapsed binary traitMat form, which have at least a min of percDiffs (percent differences)
		for (thisSite in colnames(eval(as.name(alnName)))[which(grepl(thisGene,colnames(eval(as.name(alnName)))))]){
			polyPoint <- which(colnames(eval(as.name(alnName))) == thisSite)
			# Since we only keep certain sites, we will not bother to record a position which is in our alignment YET is not within out siteID
			if(is.element(polyPoint,unlist(siteID))){
				tmp_refSite <- 1; found = FALSE
				while(!found && tmp_refSite <= length(siteID)){
					 if(is.element(polyPoint,siteID[[tmp_refSite]])){
					 	found = TRUE
					 } else {
					 	tmp_refSite <- tmp_refSite + 1
					 }
				}
				# This will exit our script if there is a problem and report it
				if(!found){ 
					print(paste("There was a problem finding the site for ",thisGene,thisSite,sep=" "))
					q(save="no")
				}
				# This is a sanity check
				if(!all(traitMat[which(traitMat[,polyPoint] != 2),polyPoint] == builtMat[which(builtMat[,tmp_refSite] != "-"),tmp_refSite])){
					print(paste("There was a problem associating our builtMat and traitMat for ",thisGene,thisSite,sep=" "))
					q(save="no")
				}
				# We evaluate if the traitMat has a minimum number of differences, this is identified by nSeq and the percDiffs
				# objects which the max function returns the number of rows in the alnName object that must be polymorphic
				# This is a sanity check sa it should be by defacto true
				if ( minDiffs(builtMat[,tmp_refSite], numDiff) ){
					# If this site is polymorphic to a degree greater than our percDiffs * nSeq we want to record it
					returnFrame <- rbind(returnFrame, c(thisGene, 
									substr(thisSite,max(gregexpr("_",thisSite)[[1]])+1,nchar(thisSite)),
									paste(eval(as.name(alnName))[,thisSite],collapse=""), 
									paste(builtMat[,tmp_refSite],collapse="_")) )	
				}
			}
		}
		# We return the returnFrame elements which are different from initialisin
		return( returnFrame[!is.na(returnFrame$binary_Seq),] )
	}
	
	# Now we save the information we've generated to be passed to downstream functions
	save(polyFrame,file=paste(workDir,alnName,"_polyFrame.RData",sep=""))
} else {
	print(" The traitMat and alnName objects did not have identical correspondance between rownames did not proceed ")
}

# This is another sanity check
if(!all(sapply(1:nrow(polyFrame),function(x){ 
		all(eval(as.name(alnName))[,paste(polyFrame$Gene[x],polyFrame$int_ntPos[x],sep="_")] == strsplit(polyFrame$nt_Seq[x],"")[[1]]) 
	})) ){ print(" alnName and polyFrame nt_Seq information are not matching perfectly please review!") } 

###################################################################################################
######################## END OF REPORT and POLYMORPHISM CREATE ####################################	
###################################################################################################



#################################################################################################################
########################################## R BATCH CMD END FILE ARGUMENT ########################################
#################################################################################################################

quit(save="no")