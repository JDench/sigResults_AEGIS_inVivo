# This is a script to take a postAnalysis output from our PA correlated evolution BayesTraits anaylses and get some biological support

########## DEPENDENCIES ##########
# This uses output from PostAnal_BayesTraits.v.2.r and therefore also requires it's associated dependecies.

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
library(nlme)
library(doMC)
library(foreach)
library(scales) # this is to allow the use of alpha for colouring
library(geiger) # this allows us to use the node.leaves function
library(RColorBrewer)

###################################################################################################
##################################### BEGIN OF PREDEFINITIONS #####################################
###################################################################################################
farmHPCVL <- FALSE
###### System Folder settings ##########
# This is the basic directory within which we will find all other information
baseDir <- if(farmHPCVL){ "/home/hpc3058/" } else { "/Users/Jdench/Desktop/" }
# This is the main directory under which multiple sub directories may be found
mainDir <- paste(baseDir,if(farmHPCVL){ "jonathan/PARuns/Sequences/" } else { "PARuns/Sequences/" },sep="")
# Then we define the database directory where we find stored resources
dbDir <- paste(mainDir,"PA14DBase/",sep="")
# Then we define our working directories where changes will occur and be written
workDir <- paste(mainDir,"PA_FixedSite/newTrials/polyAln_0.8_retry/",sep="")
postDir <- paste(workDir,"postAnal/",sep="")
tmpDir <- paste(postDir,"tmpInfo/",sep="")
if(!dir.exists(tmpDir)){ system(paste("mkdir ", tmpDir,sep="")) }
outDir <- paste(postDir,"pValue_Only/",sep="")
if(!dir.exists(outDir)){ system(paste("mkdir ", outDir,sep="")) }
# This is only used for the isSynonymous function, it should be where there are 
# alignment files
seqDir <- paste(mainDir,"/PA_wholeGenomes/newTrials/rawAln/",sep="")

###### BLAST FOLDERS AND SETTINGS ##########
# This is the folder in which our blast program's bin folder is found
blastDir <- paste(baseDir,"Blast/ncbi-blast-2.2.30+/",sep="")
# This is the explicit location in which we will find the blast programs we might want to run
blastLoc <- paste(blastDir,"bin/",sep="")
# This is the location and name of our blast databases against which we will identify gene locations
blastDB <- paste(blastDir,"db/cGene_Alignment/wholeAln_Kos",sep="")
# This secondary database is used to asses the internal nucleotide positions of sites based on the PA14 reference 
blastDB_intPos <- paste(dbDir,"blastDB_pa14/pa14",sep="")
# This is the location of blast databases and where we will write our temporary ones
blastDB_tmp <- paste(baseDir,"db/Temp/",sep="")
# These are filesnames use for tmp files introduced in the polyFrame construction stage for int_ntPos identification
fastaOut <- "tmpCutout.fas"
blastOut <- "int_ntPos_blast.table"
# These are blast call setups I use to standardise the process and make changes easier.
outFmt <- ' -outfmt "6 qseqid sseqid length qstart qend sstart send evalue pident" '


######## File Information ##########
# This is the name of the data file which is a table of the gene, accessions and other information about the genomes we will work with.
pa14dbData <- "pa14_genicDB.RData"
# This identifies the tree file standard name if we are to perform searches on the tree
treeFile <- "wholeAln_reduced.nexus"
# This is the q-value threshold acceptance criteria and the filename of the associated output file we'll need for analysis.
ourAccept <- 1e-04
file_postAnal <- paste("sitesInfo_",ourAccept,".RData",sep="")
# These are the names of alignment files and processed information used downstream.  Changing these should allow all downstream parts
# to handle modified information formats.
alnName	<- "polyAln"
# This is the name of a file where we will store all previously explored gene interactions as taken from out geneInteraction_detailed object
knownInteractions <- "knownInteractions"
# This is a save filename that we'll use for a large object reviewing the signficance of site pairs based on the number of expected genes in the pair
fileName_expGene_numSigs <- paste("numSigs_byExpectation_",ourAccept,".RData",sep="")
# This is where we can find a copy of the catfasta2phyml.pl script
catfasta_loc <- paste(baseDir,"Jonathan_Dench/3rd_Party_Software/catfasta2phyml.pl",sep="")
# This is where fastTree can be found
fastTree_loc <- paste(baseDir,"Jonathan_Dench/3rd_Party_Software/FastTree2.1.7/fasttree",sep="")

######## Constants ##########
# This setups how many cores we are calling for parallellisation
registerDoMC(cores=4)
# We initialise the geneType boolean
geneType <- 1
# This is a percent to identify the minimum number of differences that we will require to consider a site pattern
percDiffs <- 0.05
# This is the majority consensus threshold we are setting for our assignment of character states, typically 50%
consensVal <- 0.2
# This is the minimum gene interaction score that we require for a pair to be considered as having biological support
min_intScore <- 600
# This is the sites in our unique_sigPairs data.frame that we want to search for relation.
siteRef_colSet <- c("siteRef_1","siteRef_2")
# This is a string used to separate pairs of sites or genes which get concatenanted into a single string
pairString <- "__:__"
# This is a logical concerning how we will apply some filtering
use_parBounds <- FALSE
# This is a logical of whether or not we want to use pre-existing biological support to subset
use_supportSTRING <- FALSE
# This is a logical to know if we are only considering pairs which are nonSynonymous mutations.
use_nonSyn_only <- TRUE

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

# These are the genes which we might expect to show a signal of correlated evolution
testedGenes <- list("expect"= c("gyrA", "gyrB", "parC", "parE", "morA", "nfxB"),
					"noExpect" = c("dnaA", "dnaN", "lpd3", "ribD", "rpoB", "serC"))

# These are particular genes and codons thereof which I'd like to assess for being polymorphic
list_geneCodons <- list("gyrA"=c(54,67,83:87),
						"gyrB"=c(429:585),
						"parC"=c(82:84,87,88,91,95),
						"parE"=c(357:503),
						"morA"=c(563,975,1056,1109,1155,1162,1213))
# These are particular genes and codons thereof which I'd like to assess for loss of function mutations
list_lossGenes <- c("nfxB","orfN")

# These are some constants for my testStatistic distribution subsetting format
# This is a logical if I will be farming this run out on HPCVL  
## NOTE: if(TRUE) This wil require that I have a script on HPCVL in this postDir than can be re-written and farmed about.
##			PLEASE REFER TO THE SECTION BELOW ####

hpcvlScript <- "siteReview_newFilter.v.2_mod2.r"
num_farmParts <- 4
num_hpcvlCores <- 12
hpcvl_farmOut <- paste(postDir,"HPCVLout/",sep="")
if(!dir.exists(hpcvl_farmOut)){ system(paste("mkdir ", hpcvl_farmOut,sep="")) }

# The number of replicate sitePattern combinations to test for a chi2 distribution
min_testReps <- 500
# From our test replicates how many of those do we wish to pass to each Kolmogorov-Smirnof and Anderson-Darling test
### NOTE: This value * numReps MUST <= min_testReps
tmp_numSamples <- 100
# How many tests shall we call and then average for confidence 
tmp_numReps <- 5
# This is the range of chi2 degrees of freedome distributions for which we will search for fit to our distribution of test stats.
dfRange <- 1:4
# This is the alpha level for acceptance in stats tests
ourDist_pAccept <- 0.01

# This is where the BayesTraits program can be found
bayesLoc <- if(farmHPCVL){ paste(baseDir,"bin/BayesTraitsV2/",sep="") } else { paste(baseDir,"Jonathan_Dench/3rd_Party_Software/BayesTraitsV2-Beta-OSX/",sep="") }
# This is an adjustment used for calls to the Unix operation wc, it is system based from experience.
wcAdj <- if(farmHPCVL){ 2 } else { 1 }


###################################################################################################
################################### END OF PREDEFINITIONS #########################################
###################################################################################################


###################################################################################################
##################################### BEGIN OF FUNCTIONS ##########################################
###################################################################################################

# This is a function to return the overlap between sets, it requires a named vector, list, matrix or data frame type object.
# The inSets should be defined as a vector of the set identifiers we want compared.
# The byCol argument tells us the dimension to search for the inSets definitions if the data type is matrix or data.frame
# the in_ mainSet and otherSet variables should be strings that match the set references which you want to compare, where the main set
# will be compared against otherSets, these two variables can have overlap.
findOverlap <- function(inData, inSets = NULL, byCol = TRUE, in_mainSet = NULL, in_otherSet = NULL){
	# First we define what is the data type, The order here matters since matrix and data frame also apply to lower types querried here.
	tmpType <- if(is.matrix(inData)){ 
		"matrix" 
	} else if (is.data.frame(inData)) { 
		"data.frame" 
	} else if (is.list(inData)) { 
		"list" 
	} else if (length(names(inData)) > 0){
		"named"
	} else {
		print("Error: findOverlap - type not suitable or identifiable")
		return ( NULL )
	}
	# Now if the inSets were not define explicitly we inherit them as all the names of the data.
	if (is.null(inSets)){
		# This handles the data type and defines the sets  
		if (is.element(tmpType,c("matrix","data.frame"))){
			if (byCol){ inSets <- colnames(inData) } else { inSets <- rownames(inData) }
		} else if(is.element(tmpType,c("list","named"))){
			inSets <- names(inData)
		}
	}
	# This checks that we've identified at least two sets otherwise we return an error
	if (length(inSets) < 2){ 
		print("Error: findOverlap - less than 1 set was defined or identifiable from inData")
		return( NULL )
	}
	
	# Now we create vectors, using the inSets identifiers, of all the elements related to that type
	setLists <- lapply(inSets,function(thisSet){ 
		# This handles the data type and returns the objects that will make it's vector 
		if (is.element(tmpType,c("matrix","data.frame"))){
			if (byCol){ return(inData[,thisSet]) } else { return(inData[thisSet,]) }
		} else if(is.element(tmpType,c("list","named"))){
			return( unlist(inData[which(names(inData) == thisSet)]) )
		}
	})
	
	# Now let's review if the user has defined any parituclar mainSet to use, if not we return full pairwise indexes
	mainIndex <- if (is.null(in_mainSet)){
						1:(length(setLists) - 1)
					} else {
						# Otherwise the mainIndex will be those values of the set elements 
						unlist(sapply(in_mainSet,function(x){ which(inSets == x) }))
					}
	
	# Now we will create all setwise comparissons and return a list of the intersecting objects
	tmpReturn <- list()
	for(mainSet in mainIndex){
		# Now we review and define the otherSet based on user Input
		otherIndex <- if (is.null(in_otherSet)){
						(mainSet + 1):length(setLists)
					} else {
						# Otherwise the mainIndex will be those values of the set elements 
						setdiff(unlist(sapply(in_otherSet,function(x){ which(inSets == x) })), mainIndex)
					}
		# Now that we know the indexes of the otherSet we use the expand.grid function to create all possible combinations of these values
		indexSets <- expand.grid(lapply(otherIndex, function(x){ return(otherIndex) }))
		# However this will create a set that has redundancy so we'll tidy it to have unique instances
		indexSets <- unique(unlist(sapply(1:nrow(indexSets),function(x){ 
			# We turn this into a string after ordering this indexSet in order to allow us to easily call the unique function
			tmpString <- unique(unlist(indexSets[x,]))
			tmpString <- tmpString[order(tmpString)]
			return( paste(tmpString, collapse="_:_") )
			})))
		# Now we unstring the indexSets and return them as lists of indexes to be evaluated
		indexSets <- lapply(indexSets, function(x){ as.numeric(unlist(strsplit(x,"_:_"))) })
		# Now we start looking at the overlap between the mainSet and all otherSets in combination defined now by the indexSets
		for(otherSet in 1:length(indexSets)){
			tmpOverlap <- setLists[[mainSet]]
			# This is our recursion of intersections
			for(tmpSet in indexSets[[otherSet]]){
				tmpOverlap <- intersect(tmpOverlap, setLists[[tmpSet]])
			}
			# Now that we've gone through all the indexes of otherSet we will return this 
			tmpReturn[[paste(inSets[c(mainSet, indexSets[[otherSet]])],collapse="_vs_")]] <- tmpOverlap
		}
	}
	# Now we can return the intersects found
	return( tmpReturn )
}

# This is the translate function from seqinr, but we needed to fix part of the code
# So I've copied it here and will make a simple change
translate <- function (seq, frame = 0, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE) {
    #print("My version")
    if (any(seq %in% LETTERS)) {
        seq <- tolower(seq)
    }
    # I added a line to cahnge NA back into gap characters, this is because any NA passed to s2n causes a crashout
    if (sens == "R") {
        seq <- comp(rev(seq), ambiguous = ambiguous)
        seq[is.na(seq)] <- "-"
    } 
    seqn <- s2n(seq, levels = s2c("tcag"))
    l <- 3 * ((length(seq) - frame)%/%3)
    c1 <- seq(from = frame + 1, to = frame + l, by = 3)
    tra <- 16 * seqn[c1] + 4 * seqn[c1 + 1] + seqn[c1 + 2] + 
        1
    code <- s2c(SEQINR.UTIL$CODES.NCBI$CODES[numcode])
    result <- code[tra]
    result[is.na(result)] <- NAstring
    if (ambiguous) {
        toCheck <- which(result == NAstring)
        for (i in toCheck) {
            codon <- seq[c1[i]:(c1[i] + 2)]
            allcodons <- as.vector(outer(as.vector(outer(amb(codon[1]), 
                amb(codon[2]), paste, sep = "")), amb(codon[3]), 
                paste, sep = ""))
            allaminoacids <- sapply(allcodons, function(x) translate(s2c(x), 
                numcode = numcode, ambiguous = FALSE))
            if (all(allaminoacids == allaminoacids[1])) 
                result[i] <- allaminoacids[1]
        }
    }
    return(result)
}

# This function is used to determine if a polymorphic nucleotide position has synoynymous or non-synonymous changes.
# tmp_ntPosition is assumed to be of the format <gene>_<ntPositon>, where nt_Position is from left to right of related 
# sequences
isSynonymous <- function(tmp_ntPosition){
	# We extract the gene and nucleotide position associated to tmp_ntPosition, this nested
	# strsplit function is done so that if a site is a single pattern assocaited to multiple sites, we'll simply take the first herein.
	# This could be handled better, but we aren't overly interested in these non-unique patterned sites.
	tmp_ntPosition <- strsplit(strsplit(tmp_ntPosition,"::")[[1]],"_")[[1]]
	# Now we can use the pa14_genibDC$geneID column to quickly grab the information of Locus_Tag and Strand and Seq_nt for this gene
	refInfo <- unlist(pa14_genicDB[which(pa14_genicDB$geneID == tmp_ntPosition[1])[1],c("Locus_Tag", "Strand", "Seq_nt")])
	# From our refInfo we create a logical of whether or not the strand is reversed, simply an easier call downstream
	tmp_isReverse <- as.numeric(refInfo[2]) == -1
	# We now load the alignment file for this gene, this ought to be unique, but R will crash likely if not.
	tmpAln <- read.dna(list.files(path=seqDir,pattern= refInfo[1], full.names=TRUE),format="fasta",as.character=TRUE)
	# We now remove all sequences in this alignment that are not part of our analysis, this is determined using the tTree object 
	tmpAln <- tmpAln[which(unlist(sapply(rownames(tmpAln),function(x){ is.element(x, tTree$tip.label) }))), ]
	
	
	
	#### THIS bit of code was proposed for adding colnames to the tmpAln, but is strictly not required....
	# We check that the information in wholeAln_index object return the same amount as the tmpAln size
	#if(ncol(tmpAln) != length(which(wholeAln_index$Locus_Tag == refInfo[1]))){
	#	stop(paste("The size of the tmpAln for ", refInfo[1]," does not match the amount of information in wholeAln_index, please review",sep=""))
	#}
	# Now we use the rows of the wholeAln_inex object to extrac the siteRef_ID information
	#colnames(tmpAln) <- sapply(which(wholeAln_index$Locus_Tag == refInfo[1]),function(indexRow){ 
	#						# We check if this is a ref_ntPos or insert_ntPos and then return the proper information
	#						return( paste(wholeAln_index$Gene[indexRow], if(is.na(wholeAln_index$Ref_ntPos[indexRow])){ 
	#																			paste("insertPos",wholeAln_index$Insert_ntPos[indexRow],sep="")
	#																		} else { 
	#																			wholeAln_index$Ref_ntPos[indexRow] 
	#																		}, sep = "_") )
	#					})
	
	
	# We will identify any positions in this gene which our index defines as being insert positions, these are removed
	# since our <ntPosition> information is only relative to the reference genome and does not account for inserts.
	tmpInserts <- which(!is.na(wholeAln_index$Insert_ntPos[which(wholeAln_index$Locus_Tag == refInfo[1])]))
	# We remove any columns that relate to insert positions
	tmpAln <- tmpAln[,setdiff(1:ncol(tmpAln),tmpInserts)]

	# This is a sanity check that we include since some genes, duplicated within the genome, cannot be sourced well and the sizes vary (ex: fabG)
	if(as.numeric(tmp_ntPosition[2]) > ncol(tmpAln)){
		return( NULL )
	}
	# We now identify the codon position based on the strand and nucleotide position
	tmpNucleotide <- if(tmp_isReverse){
						nchar(refInfo[3]) + 1 - as.numeric(tmp_ntPosition[2])
					} else {
						as.numeric(tmp_ntPosition[2])
					}
	tmpCodon <- ceiling(tmpNucleotide/3)
	
	# I'd like to return what nucletoide changes happen in this column, as well as at any other column in the codon.....
	# To simplify further steps we idenitfy if this is a 1st, 2nd or 3rd position change
	tmp_intPos <- if(tmpNucleotide %% 3 == 0){ 3 } else { tmpNucleotide %% 3 }
	# Now we extract the codon cols by using this internal position to identify which case we respond to
	# We use the strand to know whether to go forward or backward (the sign of the strand being 1 or -1 ...) through the alignment
	### NOTE: we use the tmp_ntPosition information as it refers to the linear left-right position in tmpAln
	### We add the geneID information as this should allow us to call columns as we've named them above
	tmp_codonCols <- if(tmp_intPos == 1){
							seq(as.numeric(tmp_ntPosition[2]), as.numeric(tmp_ntPosition[2]) + 2*(as.numeric(refInfo[2])))
						} else if (tmp_intPos == 2){
							seq(as.numeric(tmp_ntPosition[2]) - 1*(as.numeric(refInfo[2])),as.numeric(tmp_ntPosition[2]) + 1*(as.numeric(refInfo[2])))
						} else if (tmp_intPos == 3){
							seq(as.numeric(tmp_ntPosition[2]) - 2*(as.numeric(refInfo[2])), as.numeric(tmp_ntPosition[2]))
						} else {
							stop(paste("Nucleotide position not properly identifiable in whatAA for ", paste(tmp_ntPosition,collapse="_")))
						}
	# Sanity check that we haven't found more than three columns to work with
	if(length(tmp_codonCols) != 3){ stop(paste("Did not find three positions to assign for the codons in ",paste(tmp_ntPosition,collapse="_"),sep="")) }
	# If this strand is reverse we adjust the codon positions we'll call as well as internal ntPos
	if(tmp_isReverse){ tmp_codonCols <- tmp_codonCols[3:1]; tmp_intPos <- 4 - tmp_intPos  }
	# We now define the rows at which we observe a change in the state of our mutation of interest, this requires us to use the "traitMat" and "posNames"
	# objects so that we can find the binary states at this site, we keep the rows, corrected to be ordered like our tmpAln, are ones of mutations
	tmpOrder <- correctOrder(tmpGuide = rownames(tmpAln), tmpChange = rownames(traitMat))
	# In the event that posNames has the same instance more than once (which can happen from builds with duplicated genes... something to work on for later
	# we will only keep the first instance of a geneID being found in posNames (as we've consistently done throughout
	tmp_focalRows <- which(traitMat[tmpOrder,which(posNames == paste(tmp_ntPosition,collapse="_"))[1]] == "1")
	# Now we look at the unique nucleotide combinations which arise in the focal rows, but not at this site of interest
	tmpCodon_otherNT <- unique(tmpAln[tmp_focalRows,tmp_codonCols[-tmp_intPos]])
	# We make certain that tmpCodon_otherNT is a matrix, and if not we make it one damnit!
	if(!is.matrix(tmpCodon_otherNT)){ tmpCodon_otherNT <- t(as.matrix(tmpCodon_otherNT)) }
	# Even if this returns one element, and it must return at least that or there will be an error, quite rightfully, the shape is a matrix thus
	tmpAA <- unique(as.vector(apply(tmpCodon_otherNT, MARGIN = 1, function(tmpOthers){
				# Ok, now we need to "paste" together the other nucleotides with all unique nucleotide characters at this focal site
				return( unlist(lapply(unique(tmpAln[,tmp_codonCols[tmp_intPos]]),function(tmpFocal){
								# We now "rebuild" a codon based on the focal and other codons
								tmpVec <- rep(NA,3)
								# We first assign the position relating to our focal
								tmpVec[tmp_intPos] <- tmpFocal
								# Then since we've been using the other positions as ordered from tmpAln we simply fill in the NA values as is
								tmpVec[which(is.na(tmpVec))] <- tmpOthers
								return( translate(tmpVec, sens = if(tmp_isReverse){ "R" } else { "F" }) )
						})) ) 
				})))
	tmpReturn <- length(tmpAA[which(!is.element(tmpAA,c("-","X")))]) <= 1
	names(tmpReturn) <- paste(tmp_ntPosition,collapse="_")
	return( tmpReturn )
}


# This is basically the isSynonymous function I've written except that I want to return information like the strand sense
# nucleotide position, codon, codon position
more_mutInfo <- function(tmp_ntPosition){
	# We extract the gene and nucleotide position associated to tmp_ntPosition, this nested
	# strsplit function is done so that if a site is a single pattern assocaited to multiple sites, we'll simply take the first herein.
	# This could be handled better, but we aren't overly interested in these non-unique patterned sites.
	tmp_ntPosition <- strsplit(strsplit(tmp_ntPosition,"::")[[1]],"_")[[1]]
	# Now we can use the pa14_genibDC$geneID column to quickly grab the information of Locus_Tag and Strand and Seq_nt for this gene
	refInfo <- unlist(pa14_genicDB[which(pa14_genicDB$geneID == tmp_ntPosition[1])[1],c("Locus_Tag", "Strand", "Seq_nt")])
	# From our refInfo we create a logical of whether or not the strand is reversed, simply an easier call downstream
	tmp_isReverse <- as.numeric(refInfo[2]) == -1
	# We now load the alignment file for this gene, this ought to be unique, but R will crash likely if not.
	tmpAln <- read.dna(list.files(path=seqDir,pattern= refInfo[1], full.names=TRUE),format="fasta",as.character=TRUE)
	# We now remove all sequences in this alignment that are not part of our analysis, this is determined using the tTree object 
	tmpAln <- tmpAln[which(unlist(sapply(rownames(tmpAln),function(x){ is.element(x, tTree$tip.label) }))), ]
	
	
	
	#### THIS bit of code was proposed for adding colnames to the tmpAln, but is strictly not required....
	# We check that the information in wholeAln_index object return the same amount as the tmpAln size
	#if(ncol(tmpAln) != length(which(wholeAln_index$Locus_Tag == refInfo[1]))){
	#	stop(paste("The size of the tmpAln for ", refInfo[1]," does not match the amount of information in wholeAln_index, please review",sep=""))
	#}
	# Now we use the rows of the wholeAln_inex object to extrac the siteRef_ID information
	#colnames(tmpAln) <- sapply(which(wholeAln_index$Locus_Tag == refInfo[1]),function(indexRow){ 
	#						# We check if this is a ref_ntPos or insert_ntPos and then return the proper information
	#						return( paste(wholeAln_index$Gene[indexRow], if(is.na(wholeAln_index$Ref_ntPos[indexRow])){ 
	#																			paste("insertPos",wholeAln_index$Insert_ntPos[indexRow],sep="")
	#																		} else { 
	#																			wholeAln_index$Ref_ntPos[indexRow] 
	#																		}, sep = "_") )
	#					})
	
	
	# We will identify any positions in this gene which our index defines as being insert positions, these are removed
	# since our <ntPosition> information is only relative to the reference genome and does not account for inserts.
	tmpInserts <- which(!is.na(wholeAln_index$Insert_ntPos[which(wholeAln_index$Locus_Tag == refInfo[1])]))
	# We remove any columns that relate to insert positions
	tmpAln <- tmpAln[,setdiff(1:ncol(tmpAln),tmpInserts)]

	# This is a sanity check that we include since some genes, duplicated within the genome, cannot be sourced well and the sizes vary (ex: fabG)
	if(as.numeric(tmp_ntPosition[2]) > ncol(tmpAln)){
		return( NULL )
	# This case is when our tmpAln is no longer representing the full gene, in which case we can't handle the position and other info properly
	} else if(ncol(tmpAln) %% 3 != 0){
		return( "ERROR" )	
	}
	# We now identify the codon position based on the strand and nucleotide position
	tmpNucleotide <- if(tmp_isReverse){
						nchar(refInfo[3]) + 1 - as.numeric(tmp_ntPosition[2])
					} else {
						as.numeric(tmp_ntPosition[2])
					}
	tmpCodon <- ceiling(tmpNucleotide/3)
	
	# I'd like to return what nucletoide changes happen in this column, as well as at any other column in the codon.....
	# To simplify further steps we idenitfy if this is a 1st, 2nd or 3rd position change
	tmp_intPos <- if(tmpNucleotide %% 3 == 0){ 3 } else { tmpNucleotide %% 3 }
	# Now we extract the codon cols by using this internal position to identify which case we respond to
	# We use the strand to know whether to go forward or backward (the sign of the strand being 1 or -1 ...) through the alignment
	### NOTE: we use the tmp_ntPosition information as it refers to the linear left-right position in tmpAln
	### We add the geneID information as this should allow us to call columns as we've named them above
	tmp_codonCols <- if(tmp_intPos == 1){
							seq(as.numeric(tmp_ntPosition[2]), as.numeric(tmp_ntPosition[2]) + 2*(as.numeric(refInfo[2])))
						} else if (tmp_intPos == 2){
							seq(as.numeric(tmp_ntPosition[2]) - 1*(as.numeric(refInfo[2])),as.numeric(tmp_ntPosition[2]) + 1*(as.numeric(refInfo[2])))
						} else if (tmp_intPos == 3){
							seq(as.numeric(tmp_ntPosition[2]) - 2*(as.numeric(refInfo[2])), as.numeric(tmp_ntPosition[2]))
						} else {
							stop(paste("Nucleotide position not properly identifiable in whatAA for ", paste(tmp_ntPosition,collapse="_")))
						}
	# Sanity check that we haven't found more than three columns to work with
	if(length(tmp_codonCols) != 3){ stop(paste("Did not find three positions to assign for the codons in ",paste(tmp_ntPosition,collapse="_"),sep="")) }
	# If this strand is reverse we adjust the codon positions we'll call as well as internal ntPos
	if(tmp_isReverse){ tmp_codonCols <- tmp_codonCols[3:1]; tmp_intPos <- 4 - tmp_intPos  }
	# We now define the rows at which we observe a change in the state of our mutation of interest, this requires us to use the "traitMat" and "posNames"
	# objects so that we can find the binary states at this site, we keep the rows, corrected to be ordered like our tmpAln, are ones of mutations
	tmpOrder <- correctOrder(tmpGuide = rownames(tmpAln), tmpChange = rownames(traitMat))
	# In the event that posNames has the same instance more than once (which can happen from builds with duplicated genes... something to work on for later
	# we will only keep the first instance of a geneID being found in posNames (as we've consistently done throughout
	tmp_focalRows <- which(traitMat[tmpOrder,which(posNames == paste(tmp_ntPosition,collapse="_"))[1]] == "1")
	# Now we look at the unique nucleotide combinations which arise in the focal rows, but not at this site of interest
	tmpCodon_otherNT <- unique(tmpAln[tmp_focalRows,tmp_codonCols[-tmp_intPos]])
	# We make certain that tmpCodon_otherNT is a matrix, and if not we make it one damnit!
	if(!is.matrix(tmpCodon_otherNT)){ tmpCodon_otherNT <- t(as.matrix(tmpCodon_otherNT)) }
	# Even if this returns one element, and it must return at least that or there will be an error, quite rightfully, the shape is a matrix thus
	tmpAA <- unique(as.vector(apply(tmpCodon_otherNT, MARGIN = 1, function(tmpOthers){
				# Ok, now we need to "paste" together the other nucleotides with all unique nucleotide characters at this focal site
				return( unlist(lapply(unique(tmpAln[,tmp_codonCols[tmp_intPos]]),function(tmpFocal){
								# We now "rebuild" a codon based on the focal and other codons
								tmpVec <- rep(NA,3)
								# We first assign the position relating to our focal
								tmpVec[tmp_intPos] <- tmpFocal
								# Then since we've been using the other positions as ordered from tmpAln we simply fill in the NA values as is
								tmpVec[which(is.na(tmpVec))] <- tmpOthers
								return( translate(tmpVec, sens = if(tmp_isReverse){ "R" } else { "F" }) )
						})) ) 
				})))
	
	# We return the information now, and the changes for nt and aa are such that we look at the nt or aa characters for
	# the rows of the "1" and "0" state
	tmp_WT <- list("nt"=unique(tolower(tmpAln[-tmp_focalRows,as.numeric(tmp_ntPosition[2])])),
					"aa"=unique(toupper(apply(tmpAln[-tmp_focalRows, tmp_codonCols],MARGIN=1,function(x){
									return( translate(x, sens = if(tmp_isReverse){ "R" } else { "F" }) )
								}))) )
	tmp_Changed <- list("nt"=tolower(unique(tmpAln[tmp_focalRows,tmp_codonCols[tmp_intPos]])),
						"aa"=unique(apply(tmpAln[tmp_focalRows,tmp_codonCols],MARGIN=1,function(x){ 
															return( translate(x, sens = if(tmp_isReverse){ "R" } else { "F" }) ) 
														} )))
	# If this is the reverse strand then we change the nt information to be the compliment, this way we're consistently reporting on the coding strand
	if(tmp_isReverse){ 
		tmp_WT$nt <- comp(tmp_WT$nt)
		tmp_Changed$nt <- comp(tmp_Changed$nt)
	}
	
	tmpReturn <- data.frame("ntPos"= as.numeric(as.vector(tmpNucleotide)),
							"codon"= as.numeric(as.vector(tmpCodon)),
							"codonPos"= as.numeric(if(tmp_isReverse){ 4- tmp_intPos }else{ tmp_intPos }),
							"ntChange"= paste(paste(tmp_WT[["nt"]][which(!is.element(tmp_WT[["nt"]],c("x","-")))],collapse=","),
												as.vector(tmpNucleotide),
												paste(tmp_Changed[["nt"]][which(!is.element(tmp_Changed[["nt"]],c("x","-")))],collapse=","),sep="_"),
							"aaChange"= paste(paste(tmp_WT[["aa"]][which(!is.element(tmp_WT[["aa"]],c("X","-")))],collapse=","),
												as.vector(tmpCodon),
												paste(tmp_Changed[["aa"]][which(!is.element(tmp_Changed[["aa"]],c("X","-")))],collapse=","),sep="_"), 
							"isSyn" = length(tmpAA[which(!is.element(tmpAA,c("X","-")))]) <= 1, row.names=NULL, stringsAsFactors=FALSE)
	return( tmpReturn )
}


# This function will return the amino acid's of the relevant sequence alignment nucleotide position
# It also returns this codon and internal nucleotide position information
whatAA <- function(tmp_ntPosition){
	# We extract the gene and nucleotide position associated to tmp_ntPosition, this nested
	# strsplit function is done so that if a site is a single pattern assocaited to multiple sites, we'll simply take the first herein.
	# This could be handled better, but we aren't overly interested in these non-unique patterned sites.
	tmp_ntPosition <- strsplit(strsplit(tmp_ntPosition,"::")[[1]],"_")[[1]]
	# Now we can use the pa14_genibDC$geneID column to quickly grab the information of Locus_Tag and Strand and Seq_nt for this gene
	refInfo <- unlist(pa14_genicDB[which(pa14_genicDB$geneID == tmp_ntPosition[1])[1],c("Locus_Tag", "Strand", "Seq_nt")])
	# We now load the alignment file for this gene, this ought to be unique, but R will crash likely if not.
	tmpAln <- read.dna(list.files(path=seqDir,pattern= refInfo[1], full.names=TRUE),format="fasta",as.character=TRUE)
	# We now remove all sequences in this alignment that are not part of our analysis, this is determined using the tTree object 
	tmpAln <- tmpAln[which(unlist(sapply(rownames(tmpAln),function(x){ is.element(x, tTree$tip.label) }))), ]
	# We will identify any positions in this gene which our index defines as being insert positions, these are removed
	# since our <ntPosition> information is only relative to the reference genome and does not account for inserts.
	tmpInserts <- which(!is.na(wholeAln_index$Insert_ntPos[which(wholeAln_index$Locus_Tag == refInfo[1])]))
	# We remove any columns that relate to insert positions
	tmpAln <- tmpAln[,setdiff(1:ncol(tmpAln),tmpInserts)]

	# This is a sanity check that we include since some genes, duplicated within the genome, cannot be sourced well and the sizes vary (ex: fabG)
	if(as.numeric(tmp_ntPosition[2]) > ncol(tmpAln)){
		return( NULL )
	}

	# We now setup an aaAlignment to store the translated information
	tmp_aaAln <- t(sapply(1:nrow(tmpAln), function(thisRow){ translate(tmpAln[thisRow,], sens = if(refInfo[2] == -1){ "R" } else { "F" }) }, simplify="matrix"))  
	tmp_aaAln[is.na(tmp_aaAln)] <- "-"
	# We now identify the codon position based on the strand and nucleotide position
	tmpNucleotide <- if(as.numeric(refInfo[2]) == -1){
						nchar(refInfo[3]) + 1 - as.numeric(tmp_ntPosition[2])
					} else {
						as.numeric(tmp_ntPosition[2])
					}
	tmpCodon <- ceiling(tmpNucleotide/3)
	
	# I'd like to return what nucletoide changes happen in this column, as well as at any other column in the codon.....
	# To simplify further steps we idenitfy if this is a 1st, 2nd or 3rd position change
	tmp_intPos <- if(tmpNucleotide %% 3 == 0){ 3 } else { tmpNucleotide %% 3 }
	# Now we extract the codon cols by using this internal position to identify which case we respond to
	# We use the strand to know whether to go forward or backward (the sign of the strand being 1 or -1 ...) through the alignment
	### NOTE: we use the tmp_ntPosition information as it refers to the linear left-right position in tmpAln
	tmp_codonCols <- if(tmp_intPos == 1){
						seq(as.numeric(tmp_ntPosition[2]), as.numeric(tmp_ntPosition[2]) + 2*(as.numeric(refInfo[2])))
					} else if (tmp_intPos == 2){
						seq(as.numeric(tmp_ntPosition[2]) - 1*(as.numeric(refInfo[2])),as.numeric(tmp_ntPosition[2]) + 1*(as.numeric(refInfo[2])))
					} else if (tmp_intPos == 3){
						seq(as.numeric(tmp_ntPosition[2]) - 2*(as.numeric(refInfo[2])), as.numeric(tmp_ntPosition[2]))
					} else {
						stop(paste("Nucleotide position not properly identifiable in whatAA for ", paste(tmp_ntPosition,collapse="_")))
					}				
	# We now ask if there is polymorphism other than ambiguous ("-" and "X") characters
	tmpReturn <- list("aa"=table(tmp_aaAln[, tmpCodon]), "codonPos"=tmpCodon, "nt"=table(tmpAln[,as.numeric(tmp_ntPosition[2])]),"ntPos"= tmp_intPos, 
						"codonTable" = lapply(tmp_codonCols,function(x){table(tmpAln[,x])}),"codon" = tmpAln[, tmp_codonCols])
	return( tmpReturn )
}


# This is a function to return a vector of paste(collapsed) elements within columns of a matrix or data.frame()
# inData is the matrix/data.frame, inCols are the columns which you want concatenated, inSep is how they should be collapsed
extractPair <- function(inData, inCols, inSep = pairString){
	unique(vapply(1:nrow(inData), FUN.VALUE = vector(mode="character",length=1), FUN = function(tmpRow){
											paste(inData[tmpRow, inCols],collapse=inSep) }))
}

# This is a simple little function to help me define the correct order for elements
# This will return the ordering of tmpChange's indexes to be like tmpGuide
correctOrder <- function(tmpGuide,tmpChange){
				return( unlist(lapply(tmpGuide, function(x){ which(tmpChange == x) })) )
}

# This is a function to perform a recursive intersection call on a list of information
recursiveIntersect <- function(inList){
	if(length(inList) > 1){
		retIntersect <- inList[[1]]
		for(x in 2:length(inList)){
			retIntersect <- intersect(retIntersect,inList[[x]])
		}
		return( retIntersect )
	} else {
		# If there is only one elements than we return that list element as it intersects itself by default.
		return( inList )
	}
}

# This will identify if all the objects in a vector fed in as v, are identical
all_identical <- function(v) all(sapply( as.list(v[-1]),FUN=function(z) {identical(z, v[1])}))

# This is a function we send to foreach that will allow us to return lists and combine them
comb <- function(...) {
	# We create a list of lists out of the input given to us, this assumes the arguments input can be forced to a list.
	return( lapply(list(...),function(x){ as.list(x) }) )
}

# This is a function to find rows where the contents of two columns appear together in alternate rows
find_duplicateRows <- function(inFrame, col1, col2){
	returnList <- list()
	# We create a list of all the unions of col1 and col2 within inFrame
	allUnions <- lapply(1:nrow(inFrame), function(i){ return( union(inFrame[i,col1],inFrame[i,col2])[order(union(inFrame[i,col1],inFrame[i,col2]))] ) })
	for(thisRow in 	1:(nrow(inFrame)-1)){
		# we create an ordered union of thisRow's col1 and col2 elements
		tmpUnion <- union(inFrame[thisRow,col1],inFrame[thisRow,col2])[order(union(inFrame[thisRow,col1],inFrame[thisRow,col2]))]
		# Now we look for any of the ordered allUnions are identical to the ordered tmpUnion, this would mean it is an equivRow
		#tmpReturn <- paste(unique(unlist(lapply((thisRow +1):length(allUnions),function(thisUnion){ 
		#									if(all(allUnions[[thisUnion]] == tmpUnion)){ return( thisUnion ) } }))), collapse= combineStrings["equalRows"])
		tmpReturn <- unique(unlist(lapply((thisRow +1):length(allUnions),function(thisUnion){ 
											if(all(allUnions[[thisUnion]] == tmpUnion)){ return( thisUnion ) } })))
		if(!is.null(tmpReturn)){ returnList[[as.character(thisRow)]] <- tmpReturn }
	}
	return( returnList )
}

# This is a tiny little function which identifies a codon and returns the full name, three letter, or 
# single letter (default) amino acid which is represented.  The RNA or DNA versions
# can be passed, but if DNA is passed we're assuming the coding strand and not the template is given.
defineAA <- function(func_thisCodon, func_returnType = c("letter","three","name")){
	# We now review the func_returnType object setting it to the first value
	func_returnType <- func_returnType[1]
	if(!is.element(func_returnType,c("letter","three","name"))){
		stop(paste("Did not understand what kind of output was desired from defineAA: ", func_returnType,sep=""))	
	}
	# We now update the codon information to substitute any "u" values as "t"
	func_thisCodon <- gsub("u","t",tolower(as.vector(func_thisCodon, mode="character")))
	# This matrix is what will be used for returning the results, where the codons are the
	# rownames and the amino acid information is the column information
	func_refMatrix <- matrix(c(rep(c("A","Ala","Alanine"),4),
								rep(c("C","Cys","Cysteine"),2),
								rep(c("D","Asp","Aspartic Acid"),2),
								rep(c("E","Glu","Glutamic Acid"),2),
								rep(c("F","Phe","Phenylalanine"),2),
								rep(c("G","Gly","Glycine"),4),
								rep(c("H","His","Histidine"),2),
								rep(c("I","Ile","Isoleucine"),3),
								rep(c("K","Lys","Lysine"),2),
								rep(c("L","Leu","Leucine"),6),
								rep(c("M","Met","Methionine"),1),
								rep(c("N","Asn","Asparagine"),2),
								rep(c("P","Pro","Proline"),4),
								rep(c("Q","Gln","Glutamine"),2),
								rep(c("R","Arg","Arginine"),6),
								rep(c("S","Ser","Serine"),6),
								rep(c("T","Thr","Threonine"),4),
								rep(c("V","Val","Valine"),4),
								rep(c("W","Trp","Trptophan"),1),
								rep(c("Y","Tyr","Tyrosine"),2),
								rep(c("-","Stp","STOP"),3)),
							byrow=TRUE, ncol=3, 
							dimnames=list(c("gca","gcc","gcg","gct",
											"tgc","tgt",
											"gac","gat",
											"gaa","gag",
											"ttc","ttt",
											"gga","ggc","ggg","ggt",
											"cac","cat",
											"ata","atc","att",
											"aaa","aag",
											"tta","ttg","cta","ctc","ctg","ctt",
											"atg",
											"aac","aat",
											"cca","ccc","ccg","cct",
											"caa","cag",
											"aga","agg","cga","cgc","cgg","cgt",
											"agc","agt","tca","tcc","tcg","tct",
											"aca","acc","acg","act",
											"gta","gtc","gtg","gtt",
											"tgg",
											"tac","tat",
											"taa","tag","tga"),c("letter","three","name")))
	
	return( as.vector(func_refMatrix[func_thisCodon, func_returnType]) )
		
}

### This is a function from mages' blog for colouring things nicely
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

###################################################################################################
######################################## END OF FUNCTIONS #########################################
###################################################################################################




###################################################################################################
################################## BEGIN OF DATA LOADING ##########################################
###################################################################################################

# This first is an object with all the known gene interactions for P aeruginosa,
load(paste(baseDir,"PARuns/ppInteractions/ppInt_actions_links.RData",sep=""))
# We now save space by reducing this object to only those with our minimum score for gene interactions, then removing
# all other elements
geneInteractions <- ppInt_actions[which(ppInt_actions$score >= min_intScore),]
rm(ppInt_actions)
geneInteractions_detailed <- ppInt_links.detailed[which(ppInt_links.detailed$combined_score >= min_intScore),]
rm(ppInt_links.detailed)
# the other two are reference strain genome information
load(paste(baseDir,"PARuns/Sequences/PA01DBase/pa01_genicDB.RData",sep=""))
load(paste(baseDir,"PARuns/Sequences/PA14DBase/pa14_genicDB.RData",sep=""))

# We will now load our known gene interactions from previous runs, provided that one exists
if(file.exists(paste(baseDir,"PARuns/ppInteractions/",knownInteractions,"_minScore_", min_intScore,".RData",sep=""))){
	load(paste(baseDir,"PARuns/ppInteractions/",knownInteractions,"_minScore_", min_intScore,".RData",sep=""))
} else {
	# Otherwise we will need to create the knownInteractions data frame object, and pa14to01 connection object for later use
	assign("knownInteractions",
			data.frame(cbind(matrix(NA,ncol=ncol(geneInteractions_detailed),nrow=1,dimnames=list(NULL,colnames(geneInteractions_detailed))),"genePair_ID"=NA)),
			pos=".GlobalEnv")
	assign("known_notInteractions", vector(mode="character",length=0), pos=".GlobalEnv")
	assign("pa14to01", data.frame("loci_PA14"=NA,"loci_PA01"=NA), pos=".GlobalEnv")
}

# We now need to acquire the traitMat and polyAln colnames information to perform this search
if(length(list.files(path=workDir,pattern=paste(alnName,"_traitMat_", percDiffs *100,"perc_",sep=""),full.names=TRUE)) != 1){
	print("There was an error trying to load the traitMat object, please review")
	q(save="no")
} else {
	load(list.files(path=workDir,pattern=paste(alnName,"_traitMat_", percDiffs *100,"perc_",sep=""),full.names=TRUE))
}
if(length(list.files(path=workDir,pattern=paste(alnName,"_reduced.RData",sep=""),full.names=TRUE)) != 1){
	print("There was an error trying to load the traitMat object, please review")
	q(save="no")
} else {
	load(list.files(path=workDir,pattern=paste(alnName,"_reduced.RData",sep=""),full.names=TRUE))
}
posNames <- colnames(eval(as.name(alnName)))
rm(list=alnName)

# We will also need the tree object
if(!file.exists(paste(workDir, treeFile,sep=""))){
	print("There was a problem with loading a unqiue .nexus tree for use in file, please review")
	q(save="no")
} else {
	tTree <- read.nexus(paste(workDir, treeFile,sep=""))
}

# this loads the information from the Kos2015 paper, which contains a bacFrame object
if(!file.exists(paste(baseDir,"PARuns/Kos2015/Kos2015_table.RData",sep=""))){
	print("There was a problem with loading the Kos2015 table, please review")
	q(save="no")
} else {
	load(paste(baseDir,"PARuns/Kos2015/Kos2015_table.RData",sep=""))
}

# This loads the susceptibility information
if(!file.exists(paste(workDir,"susceptStrains.RData",sep=""))){
	print("There was a problem with loading the susceptStrain object, please review")
	q(save="no")
} else {
	load(paste(workDir,"susceptStrains.RData",sep=""))
}

###################################################################################################
################################### END OF DATA LOADING ###########################################
###################################################################################################



###################################################################################################
####################################### BEGIN OF BODY #############################################
###################################################################################################
# This is an assessement of the amount of polymorphism in genes and/or sites which have been reported to be related to the fluoroquinolone resistance I'm studying
# This uses the list_geneCodons and list_lossGenes objects. Other objects of interest susceptStrains and pa14_gennicDB
# We load our previous postAnalysis session information which will have 3 objects: sigPairs, unique_sigPairs and sitesInfo
if(file.exists(paste(postDir,sub(ourAccept,'0.01',file_postAnal),sep=""))){
	load(paste(postDir,sub(ourAccept,'0.01',file_postAnal),sep=""))
} else {
	print(" Could not find a postAnal results file, please review ")
	q(save="no")
}

# For each gene we will load the alignment we've built, then find the assocaited codon position, per strain (by removing gap characters and couting characters)
# Then reporting what is the AA of that codon.  For loss of function mutations this will be infered when there are stop codons that are introduced
aaDifferences <- sapply(names(list_geneCodons),function(thisGene){
					# This is the index in our reference database
					tmpIndex <- which(pa14_genicDB$geneID == thisGene)
					# We find out the strand sense of the alignment to know if we want to reverse it
					tmpForward <- pa14_genicDB$Strand[tmpIndex] == 1
					# We read the alignment and instantly remove all gap characters
					tmpAln <- sapply(read.dna(paste(seqDir,pa14_genicDB$Locus_Tag[tmpIndex],"_seqAln.fas",sep=""),as.character=TRUE, format="fasta", as.matrix=FALSE),function(x){ 
									tmpTrans <- x[which(!is.element(x,c("n","-")))]
									if(length(tmpTrans) > 0){
										return( translate(tmpTrans,sens = if(tmpForward){ "F" } else { "R" }) ) 
									} else {
										NULL
									}
								})
					# We remove any NULL returns
					tmpAln <- tmpAln[which(!sapply(tmpAln,is.null))]
					
					# Now for each of the codon positions we return what is the amino acid in question
					tmpCompare <- t(sapply(tmpAln,function(thisStrain){
									# We want to return the value for all 
									thisStrain[list_geneCodons[[thisGene]]]
									# If an amino acid was not in this alignment an NA will be returned, so we update this 
								}))
					colnames(tmpCompare) <- paste("AA",list_geneCodons[[thisGene]],sep="_")
					# We keep only the elements which are part of the phylogenetic tree we'd used for analysis
					tmpCompare <- tmpCompare[intersect(tTree$tip.label,rownames(tmpCompare)),]
					# So now per AA, we want to count the number of strains which are different from PA14, 
					# AND of these how many are not susceptStrains (ignoring PA14 and PA01)
					return( sapply(colnames(tmpCompare),function(thisSite){
									return( c("numDiff"=length(which(tmpCompare[,thisSite] != tmpCompare["PA14",thisSite])),
											 "diffResistants"=length(which(tmpCompare[setdiff(rownames(tmpCompare),c("PA14","PA01",susceptStrains)),thisSite] != tmpCompare["PA14",thisSite]))) )
							}) )
					
				},simplify=FALSE)
# Now the number of strains with a stop codon introduced prior to the last 10% of the gene
aaStops <- sapply(list_lossGenes,function(thisGene){
					# This is the index in our reference database
					tmpIndex <- which(pa14_genicDB$geneID == thisGene)
					# We find out the strand sense of the alignment to know if we want to reverse it
					tmpForward <- pa14_genicDB$Strand[tmpIndex] == 1
					# We read the alignment and instantly remove all gap characters
					tmpAln <- sapply(read.dna(paste(seqDir,pa14_genicDB$Locus_Tag[tmpIndex],"_seqAln.fas",sep=""),as.character=TRUE, format="fasta", as.matrix=FALSE),function(x){ 
									tmpTrans <- x[which(!is.element(x,c("n","-")))]
									if(length(tmpTrans) > 0){
										return( translate(tmpTrans,sens = if(tmpForward){ "F" } else { "R" }) ) 
									} else {
										NULL
									}
								})
					# We remove any NULL returns
					tmpAln <- tmpAln[which(!sapply(tmpAln,is.null))]
					# Now we look at which strains include a stop codon prior to the last 10% of it's alignment
					tmpCompare <- t(sapply(rownames(tmpAln),function(thisStrain){
									# We want to return the value for all 
									is.element("*",tmpAln[[thisStrain]][1:ceiling(length(tmpAln[[thisStrain]])*0.9)])
								}))
					})
					
# Then I also report on the the p-Value of any pair which includes a polymorphic site with diffResistants > 5
# Since I've found that there are no instances where nfxB or orfN contain premature stops in my alignment
# NOTE: These are only for sites with p-Value less than 10^-4
otherSites <- NULL
for(thisGene in names(aaDifferences)){
	# We need to know the strand sense of this gene
	tmpForward <- pa14_genicDB$Strand[which(pa14_genicDB$geneID == thisGene)] == 1
	# We identify the length of the reference gene, in terms of codons, we add 1 since the stored AA seq does not include the * amino acid
	tmpLength <- nchar(pa14_genicDB$Seq_aa[which(pa14_genicDB$geneID == thisGene)]) + 1
	# Thus any sites of interest will be given by positions with values > 30 in the 
	tmpSites <- colnames(aaDifferences[[thisGene]])[which(aaDifferences[[thisGene]]["diffResistants",] > 5)]
	tmpSites <- sub("AA",thisGene, tmpSites)
	# Now since these values are codons we look for any nucleotide positions within this codon
	tmpSites <- unlist(lapply(tmpSites,function(thisAA){
					tmp_numAA <- as.numeric(strsplit(thisAA,"_")[[1]][2])
					# this gene is read in the forward sense.
					if(tmpForward){
						paste(strsplit(thisAA,"_")[[1]][1], c((3* tmp_numAA -2):(3*tmp_numAA)),sep="_")
					# Otherwise this is a reverse sense gene	
					} else {
						paste(strsplit(thisAA,"_")[[1]][1], c((3* (tmpLength - tmp_numAA +1) -2):(3*(tmpLength - tmp_numAA +1))),sep="_")
					}
				}))
	# Now we look through our sigPairs and return any which contain a mutation in at least one site of interest
	otherSites  <- rbind(otherSites, sigPairs[which(apply(sigPairs[,siteRef_colSet],MARGIN=1,function(thisPair){ any(is.element(thisPair, tmpSites)) })), c(siteRef_colSet,"pValue_FDR")])
}
# I update my othersites information to be locus_tag's
for(thisCol in siteRef_colSet){
	if(any(grepl("GI:",otherSites[,thisCol]))){
		for(thisSite in which(grepl("GI:",otherSites[,thisCol]))){
			otherSites[thisSite,thisCol] <- pa14_genicDB$Locus_Tag[which(pa14_genicDB$geneID == strsplit(otherSites[thisSite,thisCol],"_")[[1]][1])]
		}
	}
}
	




# We load our previous postAnalysis session information which will have 3 objects: sigPairs, unique_sigPairs and sitesInfo
if(file.exists(paste(postDir,file_postAnal,sep=""))){
	load(paste(postDir,file_postAnal,sep=""))
} else {
	print(" Could not find a postAnal results file, please review ")
	q(save="no")
}

# As a step of general comparisson I'll look at the counts for number of pairs where the fixed site (siteRef2) is a gene
# from our expected vs. non-expected genes.  Using a contingency table we see if there is independence between gene type and
# the bin of q-value.  I'll also do the same for the amount of polymorphism that existed in the gene itself.  That information
# will be found in the polyFrame object
polymorphism_perGene <- if(file.exists(paste(workDir,alnName,"_polyFrame.RData",sep=""))){
							load(paste(workDir,alnName,"_polyFrame.RData",sep=""))
							# We simply tabulate the number of polymorphic sites per gene
							table(polyFrame$Gene)
						} else {
							NULL
						}
sigSites_perGene <- data.frame("Gene"= unlist(lapply(sigPairs[,c("siteRef_2")],function(x){ strsplit(x,"_")[[1]][1] })),
								"qValue"=sigPairs[,c("pValue_FDR")],
								stringsAsFactors=FALSE)
# I now make bins of the qValues by setting them to either small (meaning ourAccept -> < 10-7, or smaller)
sigSites_perGene$qBins <- cut(sigSites_perGene[,"qValue"], breaks = c(ourAccept,1e-7,min(sigSites_perGene[,"qValue"])))
sigSites_perGene$geneBins <- is.element(sigSites_perGene[,"Gene"], testedGenes$expect)
sigSites_perGene$geneBins <- factor(unlist(lapply(sigSites_perGene$geneBins,function(x){ if(x){"expect"}else{"other"} })))
# I now build the contingency table and then test for independence of the p-value bin and the gene type
sigContingency <- table(sigSites_perGene[,c("qBins","geneBins")])
sigContingency
chisq.test(sigContingency)
# I find that this fails a chisq.test, and from review of the table this shows that there are fewer
# lowly significant pairs among non-expected genes, than would be expected by chance (qualitatively).

# Next is to test the "normalised" number of significant pairs, given gene type
normTable <- NULL
if(!is.null(polymorphism_perGene)){
	normTable <- table(sigSites_perGene$Gene)
	for(thisElement in names(normTable)){
		normTable[thisElement] <- normTable[thisElement]/if(is.element(thisElement,names(polymorphism_perGene))){
																polymorphism_perGene[thisElement]
														} else {
															1
														}
	}	
}
# We perform a simply t.test asking if the mean number of normalised sigPairs differs between types that are 
# of expected genes or not
diffSets <- list("expect" = normTable[which(unlist(lapply(names(normTable),function(x){ is.element(x, testedGenes$expect) })))],
				"noExpect" = normTable[which(unlist(lapply(names(normTable),function(x){ is.element(x, testedGenes$noExpect) })))])
diffTest <- t.test(x = diffSets$expect,
					y = diffSets$noExpect)
# From a Welch Two-Sample t.test, I find no evidence that there is a difference in the normalised number of sig pairs 
# (ie total divided by the number of polymorphic sites in a gene)
diffSets
diffTest

# This test is to see if there is a correlation between the number of polymorphic sites in the "fixed" gene and the number of 
# correlated pairs it forms.
corrFrame <- matrix(c(table(polyFrame$Gene),table(sigSites_perGene$Gene)[names(table(polyFrame$Gene))]),ncol=2,
					dimnames = list(names(table(polyFrame$Gene)),c("polySites","sigPairs")))
cor.test(x = corrFrame[,"sigPairs"], y = corrFrame[,"polySites"], method = "pearson")


# Now I want to report on the number of pairwise comparisons that were effectively calculated in the analysis, this does not consider 
# the siteID length, but instead considers that when two sites share the same pattern, while we technically only perform on calculation,
# the number of pairwise comparisons calculated does not actually change... only the number of calculations performed.
print("Number of total pairwise comparisons (if we consider that identical siteID are still tested but not calculated separately)")
nrow(polyFrame) * (ncol(traitMat) - 1)

# Step one will be to review what is the breakdown of the number for unique_sigPairs with decreasing log10 pValue_FDR
sigLevels <- list()
tmpAccept <- ourAccept
while(length(which(unique_sigPairs$pValue_FDR <= tmpAccept)) > 0 && !all(unique_sigPairs$pValue_FDR[which(unique_sigPairs$pValue_FDR <= tmpAccept)] == 0) ){
	sigLevels[[as.character(tmpAccept)]] <- length(which(unique_sigPairs$pValue_FDR <= tmpAccept))
	tmpAccept <- tmpAccept/10
}
# This saves our original unique_sigPairs list, this is used quite downstream but helps with visualisation of qValue distributions
all_sigPairs <- unique_sigPairs

# Here I'm going to build an object that identifies if nucleotide positions represent synonymous changes
if(!file.exists(paste(workDir,"site_isSyn_nonInsert_polySites.RData",sep=""))){
	# This last is used by the identification of synonymous and non-synonymous mutations
	if(!file.exists(sub("rawAln/","wholeAln_index.RData",seqDir))){
		print("There was a problem with loading the wholeAln_index.RData file, please review")
		q(save="no")
	} else {
		load(sub("rawAln/","wholeAln_index.RData",seqDir))
	}
	# Using the positions in our polyAln we calculate the number of which are synonymous mutations
	site_isSyn <- foreach(thisSite = posNames[which(!(grepl("Insert",posNames)))], .combine="c")%dopar%{ isSynonymous(thisSite) }
	save(site_isSyn, file=paste(workDir,"site_isSyn_nonInsert_polySites.RData",sep=""))
	rm(wholeAln_index)
} else {
	load(paste(workDir,"site_isSyn_nonInsert_polySites.RData",sep=""))
}

# this is a test if the number of substitutions in a gene is correlated to the gene length, but to sub-divide this question for the
# expected genes, non-expected but focal genes, and all others.
numSubs_frame <- table(unlist(lapply(strsplit(posNames,"_"),function(x){ x[1] })))
numSubs_frame <- data.frame("gene"=names(numSubs_frame),
							"count"= as.vector(numSubs_frame),
							"type"=factor(unlist(lapply(names(numSubs_frame),function(x){ 
																if(is.element(x,testedGenes$expect)){
																	"expect"	
																} else if(is.element(x,testedGenes$noExpect)){
																	"noExpect"
																} else {
																	"other"
																} })),levels=c("expect","noExpect","other")))
numSubs_frame$length <- unlist(lapply(as.character(numSubs_frame$gene),function(x){
									tmpRows <- ifelse(is.element(x,pa14_genicDB$geneID),which(pa14_genicDB$geneID == x),which(grepl(x,pa14_geniDB$Molecule_Xrefs))) 
									pa14_genicDB$End[tmpRows] - pa14_genicDB$Start[tmpRows] +1 
							}))
numSubs_frame$normalised <- numSubs_frame$count/numSubs_frame$length

numSubs_compare <- sapply(levels(numSubs_frame$type),function(x){ 
								tmpRows <- which(as.character(numSubs_frame$type) == x)
								return( c("Mean"=mean(numSubs_frame$normalised[tmpRows]),"SD"=sd(numSubs_frame$normalised[tmpRows])) ) 
					})
# This is a simple ANCOVA 															
numSubs_lm <- lm(normalised ~ type,data=numSubs_frame)
anova(numSubs_lm)
summary(numSubs_lm)
library(car)
leveneTest(numSubs_lm)
hist(rstandard(numSubs_lm),breaks=100)
# I find no evidence that the type of gene, expected, focal not expected, or other has any difference in the number of substitutions found
# once we normalise for the length of the gene.  I can't reject Levene's test for homoscedacity, and that rstandard are ~ normal with only some outliers existing
# which are affecting the other's set




# Now I also care to look at if there is a correlation between the parC deleterious and beneficial mutations, using Phylogenetic least squares, 
parC_PGLS <- list()
# We define the mutations and find their traitMat location using the posNames object
parC_muts <- unique(unlist(all_sigPairs[which(all_sigPairs$pValue_FDR < 1e-11),siteRef_colSet]))
parC_muts <- parC_muts[which(grepl("parC", parC_muts))]
if(is.element("parC_2006",parC_muts)){ parC_muts <- parC_muts[-which(parC_muts == "parC_2006")] }
parC_mutTraits <- data.frame(traitMat[,sapply(parC_muts,function(x){ which(posNames == x) })])
colnames(parC_mutTraits) <- parC_muts
parC_mutTraits <- parC_mutTraits[tTree$tip.label,]
rownames(parC_mutTraits) <- tTree$tip.label
# Now we ask if any of the other three sites predict the deleterious mutation which is found to be parC_712
summary(gls(parC_712 ~ ., correlation = corBrownian(phy= tTree), data= parC_mutTraits))
# This shows that there is correlation, using this model, with 2/3 of the parC mutations, but not the same 2/3
# as predicted by BayesTraits.

# I now just want a simple frequency at which the parC_712 mutation occurs along with all other strong paired parC mutations
tmpStrains <- which(parC_mutTraits[,"parC_712"] == 1)
print(paste("The percent of times that the parC_712 mutation occurs with all 3 other strong paired parC mutations is: ",
		length(which(apply(parC_mutTraits[tmpStrains,],MARGIN=1,function(x){ sum(unlist(x)) }) == 4))/length(tmpStrains),
		" out of a total number of strains: ",
		length(tmpStrains),
		sep=""))
	
# To compare with correlation that may be because of hitchhiking we count the number of strains which poses the parC_2006 mutation
# and then all other 4 mutations
parC_muts <- unique(unlist(all_sigPairs[which(all_sigPairs$pValue_FDR < 1e-11),siteRef_colSet]))
parC_muts <- parC_muts[which(grepl("parC", parC_muts))]
parC_mutTraits <- data.frame(traitMat[,sapply(parC_muts,function(x){ which(posNames == x) })])
colnames(parC_mutTraits) <- parC_muts
tmpStrains <- which(parC_mutTraits[,"parC_2006"] == 1)
print(paste("The percent of times that the parC_2006 mutation occurs with all 4 other strong paired parC mutations is: ",
		length(which(apply(parC_mutTraits[tmpStrains,],MARGIN=1,function(x){ sum(unlist(x)) }) == 5))/length(tmpStrains),
		" out of a total number of strains: ",
		length(tmpStrains),
		sep=""))
tmpStrains <- which(parC_mutTraits[,"parC_712"] == 1)
print(paste("The percent of times that the parC_2006 and parC_712 mutation co-occurs is: ",
		length(which(apply(parC_mutTraits[tmpStrains,c("parC_712","parC_2006")],MARGIN=1,function(x){ sum(unlist(x)) }) == 2))/length(tmpStrains),
		" out of a total number of strains: ",
		length(tmpStrains),
		sep=""))





# I want to create a plot that looks at the p-value distribution of pairs which fit differing criteria:
# The basic subsetting is the distribution of pair-state which will be plotted by side, in bins.
# We create different plots for no-subsetting; only inter-genic pairs; only when at least one expected member
plottingPairs <- all_sigPairs
plottingPairs[,"siteRef_1"] <- unlist(lapply(plottingPairs[,"siteRef_1"],function(x){ strsplit(x,"::")[[1]][1] }))
plottingPairs[,"siteRef_2"] <- unlist(lapply(plottingPairs[,"siteRef_2"],function(x){ strsplit(x,"::")[[1]][1] }))

# I now create a data.frame that is a combination of all the peices of information which I intend to plot
rootString <- c("rootP_","_dep")
# we now define the columns 
rootCols <- intersect(which(grepl(rootString[1],colnames(plottingPairs))),which(grepl(rootString[2],colnames(plottingPairs))))
# I create a single large data.frame which can be subset for plotting and stats work.
if(!file.exists(paste(workDir,"qValue_histInfo_correlationType.RData",sep=""))){
	# I need to load the alignment of polymorphisms, so that I can extract the "unique-ness" of character states 
	# between pairs, prior to transforming to binary state
	load(paste(postDir,"polyAln_reduced.RData",sep=""))
	
	
	qValue_histInfo <- data.frame("siteRef_1"=plottingPairs[,"siteRef_1"],
								"siteRef_2"=plottingPairs[,"siteRef_2"],
								"state_site1"= site_isSyn[plottingPairs[,"siteRef_1"]],
								"state_site2"= site_isSyn[plottingPairs[,"siteRef_2"]],
								"expect_site1"= unlist(lapply(strsplit(unlist(lapply(plottingPairs[,"siteRef_1"],function(y){ strsplit(y,"::")[[1]][1] })),"_"),function(x){ !is.element(x[1], testedGenes$expect) })),
								"expect_site2"= unlist(lapply(strsplit(unlist(lapply(plottingPairs[,"siteRef_2"],function(y){ strsplit(y,"::")[[1]][1] })),"_"),function(x){ !is.element(x[1], testedGenes$expect) })),
								"is_intragenic"= apply(plottingPairs[,c("siteRef_1","siteRef_2")],MARGIN=1,function(thisPair){ 
															tmpReturn <- strsplit(unlist(thisPair),"_")
															return( tmpReturn[[1]][1] == tmpReturn[[2]][1] )
															}),
								"correlationType"= apply(plottingPairs,MARGIN=1,function(thisRow){
															# We extract the pair of mutations associated to this row of the all_sigPairs object
															thisPair <- unlist(thisRow[siteRef_colSet])
															# In case we have sites which had non-unqijue site patterns we simply consider the first instance in this assessment (this really isn't likely to happen)....
															thisPair <- sapply(thisPair, function(x){ strsplit(x,"::")[[1]][1] })
															# Now we try to identify the correlation type for the pair of mutations, this will be done by looking at the most likely root state for the mutations
															# and comparing this to the tip state of the pairs, first we extract the root state (from the colname) associated to the maximum likelihood
															tmp_rootState <- names(thisRow[rootCols])[which.max(thisRow[rootCols])]
															tmp_rootState <- strsplit(tmp_rootState, "_")[[1]][2]
															# Now we create a table of the tip states within our alignment and find the maximum which is not the root state, the which command removes instances
															# where the pair has not changed from the root state
															tmpTable <- table(apply(cbind(traitMat[,which(posNames == thisPair[1])],traitMat[,which(posNames == thisPair[2])]), MARGIN = 1, function(x){ 
																					paste(x,collapse="") }))
															tmpTable <- tmpTable[which(names(tmpTable) != tmp_rootState)]			
															# Now we look if the tip states dominate the opposite of the root, this means positive correlation; otherwise this means that the changes
															# we negatively correlated.
															tmpAlt_tipState <- paste(1-as.numeric(strsplit(tmp_rootState,"")[[1]]),collapse="")
															return( if(names(which.max(tmpTable)) == tmpAlt_tipState){
																					"Positive"
																				} else {
																					"Negative"
																				} )
														}),
								"ntChars"= as.vector(apply(plottingPairs[,siteRef_colSet],MARGIN=1,function(thisRow){
													# This simply asks how many a,c,t,g characters were in the pairs of sites prior to binary state transformation
													sum(apply(polyAln[,unlist(lapply(as.vector(unlist(thisRow)),function(x){ which(posNames == x) }))],MARGIN=2, function(y){ 
															length(unique(y[which(is.element(y,c("a","c","t","g")))]))
													}))
											})),
								"qValue"= plottingPairs[,"pValue_FDR"])
	
	
	# This finds out if a pair is synonymous or not by summing and then we assign factors for the values to be extracted when plotting 								
	qValue_histInfo$pair_typeState <- apply(qValue_histInfo[,c("state_site1","state_site2")],MARGIN=1,sum)
	qValue_histInfo$pair_typeState <- factor(qValue_histInfo$pair_typeState, levels = c(0,1,2),labels = c("Non Synonymous","One Synonymous","Both Synonymous"))

	qValue_histInfo$pair_expectState <- apply(qValue_histInfo[,c("expect_site1","expect_site2")],MARGIN=1,sum)
	qValue_histInfo$pair_expectState <- factor(qValue_histInfo$pair_expectState, levels = c(0,1,2),labels = c("Both Expected","One Expected","Neither Expected"))
							
	# This is the sequence of powers of 10 which represent the distribution of qValues
	tmpSeq_qValue <- seq(log10(ourAccept),round(log10(min(qValue_histInfo$qValue[which(qValue_histInfo$qValue != 0)])),0),by = -1)
	qValue_histInfo$qBins <- cut(qValue_histInfo[,"qValue"], breaks = c(10^tmpSeq_qValue,-0.01))
	levels(qValue_histInfo$qBins) <- sub("(-0.01","[0",levels(qValue_histInfo$qBins), fixed=TRUE)
	
	# I also include the physical distance between mutations, as per the reference genome of PA14
	tmpLen_refGen <- ncol(read.dna(paste(baseDir,"/PARuns/Sequences/PA_wholeGenomes/allGenomes/PA14.fas",sep=""),format="fasta"))
	qValue_histInfo$physDist <- as.numeric(apply(plottingPairs[, siteRef_colSet], MARGIN = 1, function(thisPair){
												tmpGenes <- sapply(sapply(thisPair,function(y){ strsplit(y,"::")[[1]][1]}), function(x){ strsplit(x,"_")[[1]] })
												tmpPositions <- apply(tmpGenes, MARGIN = 2, function(x){ pa14_genicDB$Start[which(pa14_genicDB$geneID == x[1])[1]] - 1 + as.numeric(x[2]) })
												# Now that we have the positions, in PA14 reference genome of each of these mutations, we find the smallest of either
												# their absolute difference, OR the sum of the smallest and the difference between the largest and the length of the genome
												return( min(abs(diff(tmpPositions)),sum(min(tmpPositions), tmpLen_refGen - max(tmpPositions) + 1)) ) 
											}))
	
	# We now save the data object for future use.
	save(qValue_histInfo,file = paste(workDir,"qValue_histInfo_correlationType.RData",sep=""))
	# We remove this large object from our working space.
	rm(polyAln)
	
} else {
	load(paste(workDir,"qValue_histInfo_correlationType.RData",sep=""))
}




# These are some rough stats on my counts of pairs by types.  Knowing how many synonymous-non-synonymous, as well as possible intra-genic vs total pairs,
# and number of sites are in expected genes vs. total possible pairs, I can derive a NULL expectation for the number of these types:
# NOTE: I haven't thought of a statistical means to comapre this to our observations, since I have a single observation and lack information 
#		that I can interpret to be used to perform a chi-squared or fisher's exact test, which also as this is a lrge number of counts those methods are
# 		to my knowledge not best used for this size of counts set.
posNames_genes <- table(unlist(lapply(strsplit(names(site_isSyn),"_"),function(x){ x[1] })))
unexpectedGenes <- names(posNames_genes[which(!is.element(names(posNames_genes),testedGenes$expect))])
possiblePairs <- c("any"=length(site_isSyn) * (length(site_isSyn)-1) /2,
						
						"Both Synonymous"= length(which(site_isSyn)) * (length(which(site_isSyn))-1)/2,
						"One Synonymous"=as.numeric(length(which(site_isSyn))) * as.numeric(length(which(!site_isSyn))),
						"Non Synonymous"=length(which(!site_isSyn)) * (length(which(!site_isSyn))-1)/2,
						
						"Intragenic"= sum(sapply(posNames_genes,function(thisGene){ thisGene * (thisGene-1)/2 })),
						"Intergenic"= sum(as.numeric(sapply(names(posNames_genes)[-length(names(posNames_genes))],function(thisGene){ 
													posNames_genes[thisGene] * 
													sum(as.numeric(posNames_genes[(which(names(posNames_genes) == thisGene)+1): length(names(posNames_genes))]))
											}))),
						
						"Both Expected"=sum(as.numeric(c(sapply(testedGenes$expect[-length(testedGenes$expect)],function(thisGene){ 
													(posNames_genes[thisGene] * (posNames_genes[thisGene] -1)/2) +
													(posNames_genes[thisGene] * sum(posNames_genes[testedGenes$expect[(which(testedGenes$expect == thisGene) +1): length(testedGenes$expect)]]))
						 					}),
						 					posNames_genes[testedGenes$expect[length(testedGenes$expect)]] * (posNames_genes[testedGenes$expect[length(testedGenes$expect)]]-1)/2))),
						"One Expected"=sum(as.numeric(sapply(testedGenes$expect,function(thisGene){
													posNames_genes[thisGene] * as.numeric(sum(posNames_genes[unexpectedGenes]))
											 }))),
						"Neither Expected"=sum(as.numeric(c(sapply(unexpectedGenes[-length(unexpectedGenes)], function(thisGene){
													(posNames_genes[thisGene] * (posNames_genes[thisGene] -1)/2) +
													(posNames_genes[thisGene] * sum(posNames_genes[unexpectedGenes[(which(unexpectedGenes == thisGene) +1): length(unexpectedGenes)]]))
											}),
											 posNames_genes[unexpectedGenes[length(unexpectedGenes)]] * (posNames_genes[unexpectedGenes[length(unexpectedGenes)]] -1)/2))) )
# This is a list for visual convenience about the NULL expectation of the types of pairs we'd expect.
typeExpectations <- possiblePairs[-which(names(possiblePairs) == "any")]/possiblePairs["any"]
list_typeExpectations <- list("paired_mutType"= typeExpectations[1:3],
						"paired_mutPositions"=typeExpectations[4:5],
						"paired_mutExpected"=typeExpectations[6:8])
					
# This object is a sanity measure incase I later change the ordering of the above items
reportOrder <- c("Both Synonymous","One Synonymous","Non Synonymous","Intragenic","Intergenic","Both Expected","One Expected","Neither Expected")	
# This is a matrix of the observed and expected proportions of my types, I break it down further by the bin of q-value and total 
array_frequencyCompare <- matrix(cbind(sapply(c("All",levels(qValue_histInfo$qBins)),function(thisSet){
								tmpRows <- if(thisSet == "All"){ 1:nrow(qValue_histInfo) } else { which(qValue_histInfo$qBins == thisSet) }
								return( c(length(which(qValue_histInfo$pair_typeState[tmpRows] == "Both Synonymous"))/length(tmpRows),
												length(which(qValue_histInfo$pair_typeState[tmpRows] == "One Synonymous"))/length(tmpRows),
												length(which(qValue_histInfo$pair_typeState[tmpRows] == "Non Synonymous"))/length(tmpRows),
												length(which(qValue_histInfo$is_intragenic[tmpRows]))/length(tmpRows),
												length(which(!qValue_histInfo$is_intragenic[tmpRows]))/length(tmpRows),
												length(which(qValue_histInfo$pair_expectState[tmpRows] == "Both Expected"))/length(tmpRows),
												length(which(qValue_histInfo$pair_expectState[tmpRows] == "One Expected"))/length(tmpRows),
												length(which(qValue_histInfo$pair_expectState[tmpRows] == "Neither Expected"))/length(tmpRows)) )
							}), typeExpectations),byrow=TRUE,nrow=length(c("All",levels(qValue_histInfo$qBins)))+1,
								dimnames=list(c("All",levels(qValue_histInfo$qBins),"Expected"), reportOrder))



# This is a means to visualise some information concerning my pairs of mutations, their distribution of q-values, and some null expectations
pdf(paste(outDir,"pairState_barplot.pdf",sep=""),width=6*2,height=4*4)
par(mfrow=c(4,2), mar = c(4,4,3,1) + 0.4)

####### PLOT Pos-Neg Assocation #######
# This is for if the pair show a signal of Positive or Negative association (this is determined by looking at the most
# probable root state and examining what is the most common tip state
plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Positive"=log2(length(which(as.character(qValue_histInfo$correlationType[tmpRows]) == "Positive"))),
								"Negative"=log2(length(which(as.character(qValue_histInfo$correlationType[tmpRows]) == "Negative")))) )
				})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Type of Correlation", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)

####### PLOT Info-Loss #######
# This is for if the pair's number of unique a,c,g,t nucleotides were part of the alignment prior to binarisation
plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Binary"=log2(length(which(as.character(qValue_histInfo$ntChars[tmpRows]) == 4))),
								"One Trinary"=log2(length(which(as.character(qValue_histInfo$ntChars[tmpRows]) == 5))),
								"Many States" = log2(length(which(as.character(qValue_histInfo$ntChars[tmpRows]) > 5)))) )
				})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","purple1","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Nucleotide Sates", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)
		
####### PLOT Expected Genes #######
# This is for when we consider the distribution of p-values grouped by if a gene was expected to have been found
# Now to create my barplot I'll make a table of this information
plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Both Expected"=log2(length(which(as.character(qValue_histInfo$pair_expectState[tmpRows]) == "Both Expected"))),
								"One Expected"=log2(length(which(as.character(qValue_histInfo$pair_expectState[tmpRows]) == "One Expected"))),
								"Neither Expected"=log2(length(which(as.character(qValue_histInfo$pair_expectState[tmpRows]) == "Neither Expected"))) ))
				})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","purple1","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Was Mutation in an Expected Gene: Observed", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)

plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Both Expected"=log2(length(tmpRows) * array_frequencyCompare["Expected","Both Expected"]),
								"One Expected"=log2(length(tmpRows) * array_frequencyCompare["Expected","One Expected"]),
								"Neither Expected"=log2(length(tmpRows) * array_frequencyCompare["Expected","Neither Expected"]) ))
				})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","purple1","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Was Mutation in an Expected Gene: Expected", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)		

####### PLOT Non-Syn vs Syn #######
# Now to create my barplot I'll make a table of this information
plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Non Synonymous"=log2(length(which(as.character(qValue_histInfo$pair_typeState[tmpRows]) == "Non Synonymous"))),
								"One Synonymous"=log2(length(which(as.character(qValue_histInfo$pair_typeState[tmpRows]) == "One Synonymous"))),
								"Both Synonymous"=log2(length(which(as.character(qValue_histInfo$pair_typeState[tmpRows]) == "Both Synonymous"))) ))
					})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","purple1","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Type of Mutation: Observed", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)

plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Non Synonymous"=log2(length(tmpRows) * array_frequencyCompare["Expected","Non Synonymous"]),
								"One Synonymous"=log2(length(tmpRows) * array_frequencyCompare["Expected","One Synonymous"]),
								"Both Synonymous"=log2(length(tmpRows) * array_frequencyCompare["Expected","Both Synonymous"]) ))
					})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","purple1","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Type of Mutation: Expected", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)

####### PLOT Inter vs Intra genic #######
# This is a plot for if the pair is intra-genic or not
# Now to create my barplot I'll make a table of this information
plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Intragenic"=log2(length(which(qValue_histInfo[tmpRows,"is_intragenic"]))),
								"Intergenic"=log2(length(which(!qValue_histInfo[tmpRows,"is_intragenic"])))) )
				})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Is the Pair Intragenic: Observed", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)

plotTable <- sapply(levels(qValue_histInfo$qBins),function(x){
					tmpRows <- which(qValue_histInfo$qBins == x)
					return( c("Intragenic"=log2(length(tmpRows) * array_frequencyCompare["Expected","Intragenic"]),
								"Intergenic"=log2(length(tmpRows) * array_frequencyCompare["Expected","Intergenic"])) )
				})
barplot(plotTable, beside = TRUE, ylim = c(-2, 13), xlab = "p-Value (FDR Corrected)",
		ylab = "Count Log2", col= c("red3","royalblue1"), names.arg = colnames(plotTable),
		legend.text = rownames(plotTable), args.legend= list(x="topleft", bty="n", title="Pair State"), xpd = FALSE,
		main = "Is the Pair Intragenic: Expected", cex = 1.8, cex.names= 0.6, cex.axis=0.9, cex.lab = 1.4)
		
dev.off()





###### Genetic-Distance vs. -log10 P value ####
pdf(paste(outDir,"distance_vs_pValue_sigPairs.pdf",sep=""), height = 4, width = 6)
par(mfrow=c(1,1))
# Since we'll be plotting log10 of our p-Value's we need to adjust zero values to be smaller
tmpY <- qValue_histInfo$qValue
tmpY[which(tmpY == 0)] <- min(tmpY[which(tmpY != 0)])/10
tmpY <- -log10(tmpY)
# We can now plot
plot(x = NULL, y = NULL, xpd = FALSE, cex.axis = 0.8,
	xlim = c(min(qValue_histInfo$physDist),max(qValue_histInfo$physDist)), 
	ylim = c(min(tmpY),max(tmpY)),
	xlab = "Distance Between Mutations (nt positions)", ylab = "-log10 p-Value (FDR Corrected)")
# This density line will help visualise where, in physical distance space, the bulk of results lay
tmpDensity <- density(qValue_histInfo$physDist, from = min(qValue_histInfo$physDist), to = max(qValue_histInfo$physDist))
# We now scale the density values to have the same y-axis range as our plot
tmpDensity$y <- tmpDensity$y * max(tmpY)/max(tmpDensity$y)
lines(tmpDensity,col=add.alpha("red3",1))
points(x = qValue_histInfo$physDist, y = tmpY, cex = 0.1, pch = 20)
dev.off()


# This is to test if p-value is a function of physical distance, this may suggest a hitchhiking effect
tmpLen_refGen <- ncol(read.dna(paste(baseDir,"/PARuns/Sequences/PA_wholeGenomes/allGenomes/PA14.fas",sep=""),format="fasta"))
tmpData <- as.numeric(apply(unique_sigPairs[, siteRef_colSet], MARGIN = 1, function(thisPair){
												tmpGenes <- sapply(sapply(thisPair,function(y){ strsplit(y,"::")[[1]][1]}), function(x){ strsplit(x,"_")[[1]] })
												tmpPositions <- apply(tmpGenes, MARGIN = 2, function(x){ pa14_genicDB$Start[which(pa14_genicDB$geneID == x[1])[1]] - 1 + as.numeric(x[2]) })
												# Now that we have the positions, in PA14 reference genome of each of these mutations, we find the smallest of either
												# their absolute difference, OR the sum of the smallest and the difference between the largest and the length of the genome
												return( min(abs(diff(tmpPositions)),sum(min(tmpPositions), tmpLen_refGen - max(tmpPositions) + 1)) ) 
											}))
tmp_subData <- data.frame("physDist"= tmpData,
						"qValue"=-log10(unique_sigPairs$pValue_FDR + min(unique_sigPairs$pValue_FDR[which(unique_sigPairs$pValue_FDR != 0)])/1e4))						
distModel <- glm(qValue ~ physDist, data = tmp_subData, family = "Gamma")
## BUT after a bit of work I can't fit a good model to this data, most likely there is no strong correlation between qValue and distance alone
# A simple correlation test?
cor.test(tmp_subData$physDist, tmp_subData$qValue, method="pearson")



# Recently I've found that if we use the p-value as our primary filter then we find our strongest 
# pairs of mutations are found to be almost entirely intragenic, though they are within genes we might
# expect to bear mutations under selection for antibiotic resistance.
smallP_sigPairs <- all_sigPairs[which(all_sigPairs$pValue_FDR <= 1e-11),]
smallP_sigPairs[,c("siteRef_1","siteRef_2","pValue_FDR")]
# This creates our processed information concerning mutations of interest.
if(!file.exists(paste(outDir,"more_mutInfo.csv",sep=""))){
	load(sub("rawAln/","wholeAln_index.RData",seqDir))
	# It's been asked that I dive deeper into the specific nature of each change, so I want to extract if the change is
	# synonymous or not, the actual strand sense value of the nucleotide positions, and other information as appropriate.
	all_mutInfo <- matrix(unlist(sapply(unique(c(unlist(smallP_sigPairs[,"siteRef_1"]),unlist(smallP_sigPairs[,"siteRef_2"]))),function(thisMut){
									more_mutInfo(paste(thisMut,collapse="_")) },simplify="matrix")),ncol=6,byrow=TRUE,
									dimnames=list(unique(c(unlist(smallP_sigPairs[,"siteRef_1"]),
													unlist(smallP_sigPairs[,"siteRef_2"]))),c("ntPos","codon","codonPos","ntChange","aaChange","isSyn")) )
	# We add on the partner(s) for each mutation
	all_mutInfo <- cbind(all_mutInfo, unlist(lapply(unique(c(unlist(smallP_sigPairs[,"siteRef_1"]),unlist(smallP_sigPairs[,"siteRef_2"]))),function(thisMut){
																paste(unique(unlist(c(smallP_sigPairs[which(smallP_sigPairs[,"siteRef_2"] == thisMut),"siteRef_1"],
																						smallP_sigPairs[which(smallP_sigPairs[,"siteRef_1"] == thisMut),"siteRef_2"]))),collapse="__:__") })))
	colnames(all_mutInfo)[ncol(all_mutInfo)] <- "partner"
	write.csv(all_mutInfo,file=paste(outDir,"more_mutInfo.csv",sep=""),row.names=TRUE)
	rm(wholeAln_index)
} else {
	all_mutInfo <- read.csv(paste(outDir,"more_mutInfo.csv",sep=""),header=TRUE, row.names=1)	
}


# Now I'm going to actually build a CAI index for each and every one of my strains.  This will be done for each strain in the 
# tree object (tTree$tip.label).
if(!file.exists(paste(outDir,"allStrains_CAI.RData",sep=""))){
	# I create the table in which I'll be summing all counts of CAI values
	tableCAI <- matrix(0,ncol=4^3,nrow=length(tTree$tip.label),
						dimnames=list(tTree$tip.label,apply(expand.grid(c("a","c","g","t"),c("a","c","g","t"),c("a","c","g","t")),MARGIN=1,function(x){ paste(x,collapse="") })))
						
						### USING EVERY GENE IS NOT HOW CAI IS CALCULATED NORMALLY....
						## Now I report on all the sequence alignment files I have that are part of a reference set of genes.  
						## As I have no particular reference set in mind, or easily graspable, I'll just use ALL genes.
						#allAlignments <- list.files(path=seqDir,pattern='_seqAln.fas')
						
						## I've followed more traditional methods of calculating CAI by taking the ribosomal RNA and protein genes, as well as transcription factors.
						## Similar to Sharp 1987, we use ribosomal proteins, elongation factors, but also include polymerases.
						#caiGenes <- c(which(grepl("ribosomal protein",pa14_genicDB$Product)),
						#				which(grepl("polymerase",pa14_genicDB$Product)),
						#				which(grepl("transcription elongation",pa14_genicDB$Product)))
						## We remove transferase and modification proteins as well as those which are sigma factors, depolymerases and sub-family types.
						## From experience revieweing the returned sets this will give what ought to be the most highly expressed genes.
						#for(thisRemove in c("transferase","modification","subfamily","depolymerase","sigma")){
						#	caiGenes <- caiGenes[which(!grepl(thisRemove,pa14_genicDB$Product[caiGenes]))]
						#}
	
	# Now to compare my calculations of CAI with that of DAMBE following the highly expressed genes reported in Hilterbrand 2012
	caiGenes <- as.vector(unlist(read.table(paste(baseDir,"PARuns/Sequences/Pseudomonas_highly_expressed_genes.ptt",sep=""),sep="\t",row.names=NULL,header=TRUE,stringsAsFactors=FALSE)))
	# I now go and find which locus_tag elements are referred to by the gene name, but since geneID in this table is from PA01, we have a two step process
	caiGenes <- sapply(caiGenes,function(x){ which(pa14_genicDB$geneID == x) })
	
	allAlignments <- unique(sapply(pa14_genicDB$Locus_Tag[caiGenes],function(x){ list.files(path=seqDir, pattern=x) }))
					
	# Now for each and every one of these files I'll open them, check their strand sense as per the PA14 reference (against which they've been aligned)
	# remove any insert sites, identified as being a gap character in the sequence, and then cutting into codons.  I note this will be a slightly imperfect
	# means of calculating CAI for strains which have many gap characters.... but I'm not fussed right now...
	for(thisFile in allAlignments){
		# We define the locus tag, then strand sense of the alignment
		tmpLocus <- strsplit(thisFile,"_seqAln")[[1]][1]
		tmpSense <- pa14_genicDB$Strand[which(pa14_genicDB$Locus_Tag == tmpLocus)]
		# we load the file
		tmpAln <- read.dna(paste(seqDir,thisFile,sep=""),as.character=TRUE, format= "fasta")
		# If the strand is the reverse, then we update our alignment
		if(tmpSense == -1){
			tmpAln <- t(apply(tmpAln,MARGIN=1,function(x){ rev(comp(x)) }))
		}
		# This is part of my personal approach to calculating CAI
		# Now we go through each sequence, removing gaps characters, and then splitting into groups of three then tabulating
		tmpCodons <- sapply(rownames(tmpAln),function(thisSequence){
							# I put the sequence's nucleotide characters into a matrix so they can be pasted into codons
							thisSequence <- tolower(tmpAln[thisSequence,which(is.element(tolower(tmpAln[thisSequence,]),c("a","c","g","t")))])
							theseCodons <- matrix(thisSequence[1:(3*floor(length(thisSequence)/3))],byrow=TRUE,ncol=3)
							if(length(theseCodons) > 6){
								return( table(apply(theseCodons,MARGIN=1,function(x){ paste(x,collapse="") })) )
							# But if there aren't at least two codons then I'll treat this as a null result
							} else {
								NULL	
							}
						},simplify=FALSE)
		# Ok, so for each of these elements we add them to the tableCAI
		for(thisStrain in names(tmpCodons)){
			if(is.element(thisStrain,rownames(tableCAI))){
				# Using the names of the table of codons for a strain, we add it's codon values to the tableCAI object
				tableCAI[thisStrain,names(tmpCodons[[thisStrain]])] <- tableCAI[thisStrain,names(tmpCodons[[thisStrain]])] + tmpCodons[[thisStrain]]
			}
		}
		
	}
	
	# Now to make a CAI table we'll turn the counts into proportion for each amino acid respectively.
	# I can use my defineAA function to find which aa is represented by which codon.
	colnamesAA <- defineAA(colnames(tableCAI),"three")
	# To calculate the CAI values we generate the RSCU - relative codon usage index, which is just the counts for a codon
	# divided by the average expected if all codons coding for an a.a were equally likely
	arrayCAI <- array(c(tableCAI,t(apply(tableCAI,MARGIN=1,function(thisStrain){  
						# Ok to calculate the RSCU we need to for through each unique amino acid, and for each of those
						# we then calculate it's relative frequency compared to an unbiased 
						for(thisAA in unique(colnamesAA)){
							thisStrain[which(colnamesAA == thisAA)] <- thisStrain[which(colnamesAA == thisAA)] * 
																		length(which(colnamesAA == thisAA))/
																		sum(thisStrain[which(colnamesAA == thisAA)])
						}
						return(thisStrain)
					}))), 
					dim=c(nrow(tableCAI),ncol(tableCAI),2), dimnames=list(rownames(tableCAI),colnames(tableCAI),c("Counts","RSCU")))
	# Now the weights for each AA are calculated by taking the difference in proportional RSCU and counts
	arrayCAI <- array(c(arrayCAI,t(apply(arrayCAI,MARGIN=1,function(thisStrain){  
									## HERE ## I will now have a matrix of the counts and RSCU for each codon
									##			and the rownames should match the positions of the colnamesAA object.
									tmpReturn <- vector(mode="numeric",length=nrow(thisStrain))
									for(thisAA in unique(colnamesAA)){
										# I identify the indexes which reflect codons for thisAA
										tmpIndexes <- which(colnamesAA == thisAA)
										# Now the weights is simply the ratio between relative (to maximum) RSCU OR counts, I use counts...
										tmpReturn[tmpIndexes] <- thisStrain[tmpIndexes,"Counts"]/max(thisStrain[tmpIndexes,"Counts"])
									}
									return(tmpReturn)
								}))),
					dim=c(dim(arrayCAI)[1],dim(arrayCAI)[2],3), dimnames=list(dimnames(arrayCAI)[[1]],dimnames(arrayCAI)[[2]],c(dimnames(arrayCAI)[[3]],"Weights")))
	# I now save this object which can be later used.
	save(arrayCAI, caiGenes,file=paste(outDir,"allStrains_CAI.RData",sep=""))
} else {
	load(paste(outDir,"allStrains_CAI.RData",sep=""))	
}

# This is a directory for the DAMBE information
dambeDir <- paste(outDir,"infoDAMBE/",sep="")
if(!dir.exists(dambeDir)){ system(paste("mkdir -p ",dambeDir,sep="")) }
		
# Now I have a CAI value which has been calculated but I also want to create information for DAMBE, what is required is
# a .ITE file, which while following a template of PA01's provided version, I will calculate the highly expressed and lowly expressed genes (which are all non-highly expressed ones
if(!file.exists(paste(outDir,"allStrains_DAMBE_infoITE.RData",sep=""))){
	# I create the table in which I'll be summing all counts of CAI values, and we store the values for the highly expressed genes and lowly expressed ones.
	tableCAI_DAMBE <- array(0,dim=c(length(tTree$tip.label),4^3,2),
						dimnames=list(tTree$tip.label,
										apply(expand.grid(c("a","c","g","t"),c("a","c","g","t"),c("a","c","g","t")),MARGIN=1,function(x){ paste(x,collapse="") }),
										c("HEG","LEG")) )
						
						### USING EVERY GENE IS NOT HOW CAI IS CALCULATED NORMALLY....
						## Now I report on all the sequence alignment files I have that are part of a reference set of genes.  
						## As I have no particular reference set in mind, or easily graspable, I'll just use ALL genes.
						#allAlignments <- list.files(path=seqDir,pattern='_seqAln.fas')
						
						## I've followed more traditional methods of calculating CAI by taking the ribosomal RNA and protein genes, as well as transcription factors.
						## Similar to Sharp 1987, we use ribosomal proteins, elongation factors, but also include polymerases.
						#caiGenes <- c(which(grepl("ribosomal protein",pa14_genicDB$Product)),
						#				which(grepl("polymerase",pa14_genicDB$Product)),
						#				which(grepl("transcription elongation",pa14_genicDB$Product)))
						## We remove transferase and modification proteins as well as those which are sigma factors, depolymerases and sub-family types.
						## From experience revieweing the returned sets this will give what ought to be the most highly expressed genes.
						#for(thisRemove in c("transferase","modification","subfamily","depolymerase","sigma")){
						#	caiGenes <- caiGenes[which(!grepl(thisRemove,pa14_genicDB$Product[caiGenes]))]
						#}
	
	# Now to compare my calculations of CAI with that of DAMBE following the highly expressed genes reported in Hilterbrand 2012
	caiGenes <- as.vector(unlist(read.table(paste(baseDir,"PARuns/Sequences/Pseudomonas_highly_expressed_genes.ptt",sep=""),sep="\t",row.names=NULL,header=TRUE,stringsAsFactors=FALSE)))
	# I now go and find which locus_tag elements are referred to by the gene name, but since geneID in this table is from PA01, we have a two step process
	caiGenes_index <- sapply(caiGenes,function(x){ which(pa14_genicDB$geneID == x) })
	
	allAlignments <- list("HEG"=unique(sapply(pa14_genicDB$Locus_Tag[caiGenes_index],function(x){ list.files(path=seqDir, pattern=x) })))
	allAlignments$LEG <- setdiff(list.files(path=seqDir, pattern="_seqAln.fas"), allAlignments$HEG)
	# Now from the LEG, since this is codon frequency for translated proteins I should remove the tRNA or other enzymatic RNA genes
	allAlignments$LEG <- allAlignments$LEG[which(sapply(allAlignments$LEG,function(x){ 
								# Step one is to remove the file information and extract the locus tag value from this filename
								tmpLocus <- sub("_seqAln.fas","",x)
								# Now we ask if the first 5 letters of the product of this gene are tRNA-
								# I should also ask if some other strings - related to enzymatic RNA - exist, but I know of none...
								# We return FALSE for the ones we don't want....
								tmpReturn <- if(substr(pa14_genicDB$Product[which(pa14_genicDB$Locus_Tag == tmpLocus)],1,5) == "tRNA-"){
													FALSE	
												} else {
													TRUE
												}
								return ( tmpReturn )
		 }))]
	
	# I want to track the number of HEG genes that I can compute for each strain
	strainHEGs <- rep(length(caiGenes),nrow(tableCAI_DAMBE))
	names(strainHEGs) <- rownames(tableCAI_DAMBE)
					
	# Now for each and every one of these files I'll open them, check their strand sense as per the PA14 reference (against which they've been aligned)
	# remove any insert sites, identified as being a gap character in the sequence, and then cutting into codons.  I note this will be a slightly imperfect
	# means of calculating CAI for strains which have many gap characters.... but I'm not fussed right now...
	for(thisSet in names(allAlignments)){
		for(thisFile in allAlignments[[thisSet]]){
			# We define the locus tag, then strand sense of the alignment
			tmpLocus <- strsplit(thisFile,"_seqAln")[[1]][1]
			tmpSense <- pa14_genicDB$Strand[which(pa14_genicDB$Locus_Tag == tmpLocus)]
			# we load the file
			tmpAln <- read.dna(paste(seqDir,thisFile,sep=""),as.character=TRUE, format= "fasta")
			# If the strand is the reverse, then we update our alignment
			if(tmpSense == -1){
				tmpAln <- t(apply(tmpAln,MARGIN=1,function(x){ rev(comp(x)) }))
			}
			# This is part of my personal approach to calculating CAI
			# Now we go through each sequence, removing gaps characters, and then splitting into groups of three then tabulating
			tmpCodons <- sapply(rownames(tableCAI_DAMBE),function(thisSequence){
								# I put the sequence's nucleotide characters into a matrix so they can be pasted into codons
								thisSequence <- tolower(tmpAln[thisSequence,which(is.element(tolower(tmpAln[thisSequence,]),c("a","c","g","t")))])
								theseCodons <- matrix(thisSequence[1:(3*floor(length(thisSequence)/3))],byrow=TRUE,ncol=3)
								if(length(theseCodons) > 6){
									return( table(apply(theseCodons,MARGIN=1,function(x){ paste(x,collapse="") })) )
								# But if there aren't at least two codons then I'll treat this as a null result
								} else {
									NULL	
								}
							},simplify=FALSE)
			# Here is where we check if a strain has a copy of this gene by asking that at least 90% of the aligned columns aren't gaps
			if(thisSet == "HEG"){
				for(thisSequence in rownames(tableCAI_DAMBE)){
					if(length(which(is.element(tolower(tmpAln[thisSequence,]),c("a","c","g","t")))) < ceiling(ncol(tmpAln) * 0.9)){
						strainHEGs[thisSequence] <- strainHEGs[thisSequence] - 1
					}
				}
			}
			
			# Ok, so for each of these elements we add them to the tableCAI_DAMBE
			for(thisStrain in names(tmpCodons)){
				if(is.element(thisStrain,rownames(tableCAI_DAMBE))){
					# Using the names of the table of codons for a strain, we add it's codon values to the tableCAI_DAMBE object
					tableCAI_DAMBE[thisStrain,names(tmpCodons[[thisStrain]]),thisSet] <- tableCAI_DAMBE[thisStrain,names(tmpCodons[[thisStrain]]),thisSet] + tmpCodons[[thisStrain]]
				}
			}
		}
	}
	
	# Now to make a CAI table we'll turn the counts into proportion for each amino acid respectively.
	# I can use my defineAA function to find which aa is represented by which codon.
	colnamesAA <- cbind(defineAA(colnames(tableCAI_DAMBE),"letter"),
				 		defineAA(colnames(tableCAI_DAMBE),"three"))
	# DAMBE expects stop codons to have the one letter code of an asterisk.. so I update with this
	colnamesAA[which(colnamesAA == "-"),1] <- "*"
	# We now re-order the colnames, and table, to suit the ordering of DAMBE, this works as I've tested that asterisks go before letters...
	tmpOrder <- order(colnamesAA[,1])
	tableCAI_DAMBE <- tableCAI_DAMBE[, tmpOrder,]
	colnamesAA <- colnamesAA[tmpOrder,]
	
	# I can now write out the frequency found for each codon, I put this into a format which DAMBE expects for ITE files
	for(thisStrain in rownames(tableCAI_DAMBE)){
		tmpFile <- paste(dambeDir,"Pseudomonas_aeruginosa_",thisStrain,".ITE",sep="")
		cat(paste("#",c(paste("Pseudomonas aeruginosa_",thisStrain,sep=""),
						paste(" 40 HEG - ",
								paste(as.vector(unlist(read.table(paste(baseDir,"PARuns/Sequences/Pseudomonas_highly_expressed_genes.ptt",sep=""),sep="\t",row.names=NULL,header=TRUE,stringsAsFactors=FALSE))),collapse=", "),
								sep=""),
						" CF_LEG - codon frequency from all non-HEG protein coding genes",
						" All protein coding genes considered are limited to those that could be found through BLASTn of this strains' published genome",
						" be that draft or complete, by querrying PA14 reference genome sequence infromation.",
						" PA14 reference sequence, gene ID, and gene product information were obtained using the database obtained through Pseudomonas.com",
						" File: UCBPP-PA14.csv, Location: www.pseudomonas.com/downloads/pseudomonas/pgd_r_16_2/Pseudomonas/complete/gtf-complete.tar.gz",						
						paste("AA","Codon",paste("CF_",strainHEGs[thisStrain],"HEGs",sep=""),"CF_LEGs",collapse="\t")),sep=""),
			sep="\n",file=tmpFile,append=FALSE)
		for(thisCodon in 1:ncol(tableCAI_DAMBE)){
			cat(paste(paste(as.vector(colnamesAA[thisCodon,1]),toupper(gsub("t","u",colnames(tableCAI_DAMBE)[thisCodon])),sep="\t"),
						paste(as.vector(tableCAI_DAMBE[thisStrain,thisCodon,]),collapse="\t")
						,sep="\t"),sep="\n",file=tmpFile,append=TRUE)
		}
		cat("",sep="\n",file=tmpFile,append=TRUE)
	}	
	save(tableCAI_DAMBE, strainHEGs, caiGenes,file=paste(outDir,"allStrains_DAMBE_infoITE.RData",sep=""))
} else {
	load(paste(outDir,"allStrains_DAMBE_infoITE.RData",sep=""))	
}

# This section is for creating the series of .fas files that will be associated with the .ITE files for passing to DAMBE

colnamesAA <- cbind(defineAA(colnames(tableCAI_DAMBE),"letter"),
				 		defineAA(colnames(tableCAI_DAMBE),"three"))
# DAMBE expects stop codons to have the one letter code of an asterisk.. so I update with this
colnamesAA[which(colnamesAA == "-"),1] <- "*"
# We now re-order the colnames, and table, to suit the ordering of DAMBE, this works as I've tested that asterisks go before letters...
tmpOrder <- order(colnamesAA[,1])
colnamesAA <- colnamesAA[tmpOrder,]
# This is a sanity check that my colnamesAA is ordered the same as the DAMBE table
if(!all(defineAA(colnames(tableCAI_DAMBE),"three") == colnamesAA[,2])){ 
	stop("There was a problem comparing the DAMBE CAI table and our reference colnamesAA object, please review") 
}

# This creates some vectors and arrays with which to work and then format our output for analysis by DAMBE,
tmp_nonSyn_pairs <- which(apply(smallP_sigPairs[,c("siteRef_1","siteRef_2")],MARGIN=1,function(x){ all(as.logical(all_mutInfo[unlist(x),"isSyn"])) }))
tmp_indep_pairedMuts <- unique(unlist(smallP_sigPairs[tmp_nonSyn_pairs,c("siteRef_1","siteRef_2")]))
tmp_indep_pairedMuts <- tmp_indep_pairedMuts[order(tmp_indep_pairedMuts)]
# This array will be key as it will be used for writting out the files to pass to DAMBE as well as then helping us parse the resultant information.
tmp_pairedDAMBE <- array("",dim=c(nrow(tableCAI_DAMBE),length(tmp_indep_pairedMuts),2),
						dimnames=list(rownames(tableCAI_DAMBE), tmp_indep_pairedMuts,c("Codon","State")))
	
for(thisPair in tmp_nonSyn_pairs){
	# For this pair of mutations we want to load the alignment of the gene(s) being considered.
	tmpInfo <- matrix(unlist(strsplit(unlist(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")]),"_")),nrow=2, dimnames=list(c("Gene","ntPos"),NULL))
	# I also want to know the strand sense of the gene(s)
	tmpSense <- sapply(unique(tmpInfo["Gene",]),function(x){ 
						pa14_genicDB$Strand[which(pa14_genicDB$geneID == x)]
				},simplify=FALSE)
	# We read the alignment(s) related to our gene(s)
	tmpAln <- sapply(unique(tmpInfo["Gene",]),function(x){ 
						tmpReturn <- read.dna(paste(seqDir,pa14_genicDB$Locus_Tag[which(pa14_genicDB$geneID == x)],"_seqAln.fas",sep=""),as.character=TRUE,format="fasta") 
						# Now we remove any inserts which are defined by gap characters in the reference strain PA14
						if(is.element("PA14",rownames(tmpReturn))){
							tmpReturn <- tmpReturn[,which(tmpReturn["PA14",] != "-")]	
						}
						
						# If the strand sense if negative we return the reverse compliment
						if(tmpSense[[x]] == -1){
							return( t(apply(tmpReturn,MARGIN=1,function(x){ rev(comp(x)) })) )
						} else {
							return( tmpReturn )
						}
				},simplify=FALSE)
				
	tmpAdd <- array(t(sapply(rownames(tableCAI_DAMBE),function(thisStrain){
													# For this strain we review that it's within the alignment (it should be....)
													if(!all(sapply(tmpAln,function(x){ is.element(thisStrain,rownames(x)) }))){
														# If the strain doesn't exist we return NULL - this ought not to happen!
														return(tmpReturn)
													} else {
														# We go through both mutations and return what their codon is.
														tmp_mutCodon <- apply(tmpInfo,MARGIN=2,function(thisMutation){
																			# We find the columns in our alignemnt which represent the mutation's coodn
																			tmp_codonCols <- as.numeric(all_mutInfo[paste(unlist(thisMutation),collapse="_"),c("ntPos","codonPos")])
																			tmp_codonCols <- (tmp_codonCols[1] - (tmp_codonCols[2]-1)):(tmp_codonCols[1] - (tmp_codonCols[2]-3))
																			# Now I find if this codon is mutant or not
																			return( c(toupper(gsub("t","u",paste(tmpAln[[thisMutation["Gene"]]][thisStrain, tmp_codonCols],collapse=""))),
																						if(paste(tmpAln[[thisMutation["Gene"]]][thisStrain, tmp_codonCols],collapse="") != paste(tmpAln[[thisMutation["Gene"]]]["PA14", tmp_codonCols],collapse="")){
																							"Mut"
																						} else {
																							"WT"
																						})) 
																		})
														return( tmp_mutCodon )
													} })),dim=c(nrow(tableCAI_DAMBE),2,2), 
													dimnames=list(NULL,NULL,apply(tmpInfo,MARGIN=2,function(thisMutation){ paste(thisMutation,collapse="_") })))
	# Now we go and add this information to our larger paired matrix
	for(thisMut in dimnames(tmpAdd)[[3]]){
		tmp_pairedDAMBE[,thisMut,] <- tmpAdd[,,thisMut]
	}
}			

# Now for each mutation we'll create a new folder and .fas file so we can take advantage of the DAMBE multi-file options
# but we calculate the ITE for both positions, that that strain, as a WT, mut1, mut2, or mut12 sequence
# This way when we compare the values of ITE, it is for the strains which became double mutatns and we're then best reporting what adaptive evidence exists
for(thisMut in as.vector(apply(smallP_sigPairs[nonSyn_pairs,c("siteRef_1","siteRef_2")],MARGIN=1,function(x){ paste(unlist(x),collapse="__w__") }))){
	# We define what are the codons for these mutations
	mutPairs <- sapply(strsplit(thisMut,"__w__")[[1]],function(x){ 
						c("WT"=unique(tmp_pairedDAMBE[which(tmp_pairedDAMBE[,x,"State"] == "WT"),x,"Codon"]),
							"Mut"=unique(tmp_pairedDAMBE[which(tmp_pairedDAMBE[,x,"State"] == "Mut"),x,"Codon"]))
				})
	# Now we create a .fas file for each strain which has a non-NA value for this mutation, and run this loop four times , one for each combination of
	# WT-WT, mut1-WT, WT-mut2, mut1-mut2
	for(thisSet in 1:4){
		# This is a function for returning which elements of the mutPairs a set should reference
		tmp_mutVector <- function(funcSet){ 
			if(funcSet == 1){
				c(1,3)
			} else if(funcSet == 2){
				c(2,3)
			} else if(funcSet == 3){
				c(1,4)
			} else if(funcSet == 4){
				c(2,4)
			} 
		}
		# This is the directory for this set of mutations, we create it if it does not exist
		tmpDir <- paste(dambeDir, thisMut,"_set",thisSet,"/",sep="")
		if(!dir.exists(tmpDir)){ system(paste("mkdir -p ",tmpDir,sep="")) }
		# Now for each strains we go ahead and write out a sequence alignment to be passed to DAMBE.
		for(thisStrain in rownames(tmp_pairedDAMBE)){
			cat(paste(">",thisMut,sep=""), 
				paste(mutPairs[tmp_mutVector(thisSet)],collapse=""),
				sep="\n",file=paste(tmpDir,"Pseudomonas_aeruginosa_",thisStrain,".fas",sep=""),append=FALSE)	
		}
		# then I also copy a version of all .ITE files in the dambeDir into this same folder
		system(paste("cd ",dambeDir,"; cp *.ITE ",tmpDir,sep=""))
	}
	
}

###########
## RIGHT now I need to go an manually peform the analysis for each folder, in DAMBE, and then I can parse that output downstream of here....
##########

## This is where the parsing starts but again it assumes that a DAMBE analysis has been performed and exists in the dambeDir sub-folders.
# I can now go through my mutant pairs 

# this builds my parsed dat to follow the format expected for my simple plotting method....
frameDAMBE <- list()
for(thisPair in tmp_nonSyn_pairs){
	# I want to extract the individual mutation identities
	tmpMuts <- as.vector(unlist(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")]))
	#this is the pair of mutations' string notation
	tmpPair <- paste(tmpMuts,collapse="__:__")
	# We now find which strains were double mutants
	tmpStrains <- rownames(tmp_pairedDAMBE)[which(apply(tmp_pairedDAMBE[,tmpMuts,"State"],MARGIN=1,function(x){ all(x == "Mut") }))]
	# Now we remove from consideration any strains where they have gap characters in their codons.
	tmpStrains <- tmpStrains[which(!sapply(tmpStrains,function(thisStrain){ any(grepl("-",tmp_pairedDAMBE[thisStrain,tmpMuts,"Codon"])) }))]
	
	# Now we go through our four sets and collect the ITE value for the pair of sites
	tmpReturn <- matrix(t(sapply(tmpStrains,function(thisStrain){ 
						# Now we go and open the .txt file for this strain, and each mutation
						tmp_Strain <- sapply(1:4,function(thisSet){  
											tmpLines <- readLines(paste(dambeDir,sub("__:__","__w__",tmpPair),"_set",thisSet,"/Pseudomonas_aeruginosa_",thisStrain,"_Out.txt",sep=""))
											return( as.numeric(strsplit(tmpLines[which(grepl("ITE = ",tmpLines))]," = ")[[1]][2]) )
										})
						return( tmp_Strain )
				 })),ncol=4,dimnames=list(tmpStrains,c("WT","mut1","mut2","mut12")))
				 
	# We now add this to our growing list of information
	frameDAMBE[[tmpPair]] <- tmpReturn
}
	
# Right now it's time to visualise this information. To do a barplot, I want to create a table that is the mean and SD of the mut1, mut2 and mut12 values
# For all those strains that are not wild type
plotTable_DAMBE <- array(sapply(names(frameDAMBE),function(x){ 
						return( apply(frameDAMBE[[x]],MARGIN=2,function(y){ return( c("Mean"=mean(y),"SD"=sd(y))) }) )
					}),dim=c(2,4,length(names(frameDAMBE))),
						dimnames=list(c("Mean","SD"),colnames(frameDAMBE[[1]]),names(frameDAMBE)))


# With the specific information for each synonymous mutation (such as nucleotide, codon position and mutation type) I can review the CAI
# influence of my mutation pairs.  I want to report these as the independent effects, and the geometric mean of the pair.
nonSyn_pairs <- which(apply(smallP_sigPairs[,c("siteRef_1","siteRef_2")],MARGIN=1,function(x){ all(as.logical(all_mutInfo[unlist(x),"isSyn"])) }))
pairedCAI <- list()
for(thisPair in tmp_nonSyn_pairs){
	# I want to extract the individual mutation identities
	tmpMuts <- as.vector(unlist(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")]))
	#this is the pair of mutations' string notation
	tmpPair <- paste(tmpMuts,collapse="__:__")
	# We now define what are the WT and mutant codons for these mutations
	tmp_codonTypes <- matrix(gsub("u","t",tolower(sapply(tmpMuts,function(y){
									# We remove any codons which are proposed to have gaps, this should ensure we're left with 2 codons
									# otherwise there will be an error call
									tmpReturn <- c(unique(tmp_pairedDAMBE[which(tmp_pairedDAMBE[,y,"State"] == "WT"),y,"Codon"]),
									 				unique(tmp_pairedDAMBE[which(tmp_pairedDAMBE[,y,"State"] == "Mut"),y,"Codon"]))
									tmpReturn <- tmpReturn[which(!grepl("-",tmpReturn))]
									return( tmpReturn ) 				
								}))),
								ncol=length(tmpMuts),
								dimnames=list(c("WT","Mut"),tmpMuts))
	# We now find which strains were double mutants, the DAMEB related object has already identified which strains are mutant or WT
	# and we've ensured previously that the rownames are identical... but we do it again here for sanity's sake
	if(!all(rownames(arrayCAI) == rownames(tmp_pairedDAMBE))){ stop("The tmp_pairedDAMBE and arrayCAI objects are not the same order, please review") }
	tmpStrains <- rownames(tmp_pairedDAMBE)[which(apply(tmp_pairedDAMBE[,tmpMuts,"State"],MARGIN=1,function(x){ all(x == "Mut") }))]
	# Now we go through the arrayCAI weights and extract the weights for the codons of the state types, for the strains which became double mutants.
	for(thisSet in 1:4){
		# This is a function for returning which elements of the tmp_codonTypes a set should reference
		tmp_mutVector <- function(funcSet){ 
			if(funcSet == 1){
				c(1,3)
			} else if(funcSet == 2){
				c(2,3)
			} else if(funcSet == 3){
				c(1,4)
			} else if(funcSet == 4){
				c(2,4)
			} 
		}
		# This will be the object of our paired codons' CAI value for the WT, mut1, mut2, and mut12 paired set
		# The CAI is the quotient of the observed and maximal RSCU for each amino acid, as per Sharp 1987
		pairedCAI[[tmpPair]] <- sapply(1:4,function(thisSet){ 
									# We find what AA each mutation represents
									tmpAA <- defineAA(tmp_codonTypes[tmp_mutVector(thisSet)],"letter")
									# We define the AA names for the columns of our object
									tmpColnames <- defineAA(colnames(arrayCAI),"letter")
									# This will then go, strain by strain, and find the CAI by taking the quotient of the geometric mean of the
									# observed RSCU values and the maximum RSCU values for our AA under consideration
									return( sapply(tmpStrains,function(thisStrain){
												geoMean(arrayCAI[thisStrain,tmp_codonTypes[tmp_mutVector(thisSet)],"RSCU"])/
												geoMean(arrayCAI[thisStrain,sapply(tmpAA,function(thisAA){ 
																					# We define which columns relate to the amino acid in question
																					tmpCols <- which(tmpColnames == thisAA)
																					return( tmpCols[which.max(arrayCAI[thisStrain,tmpCols,"Counts"])] ) 
																			}),"RSCU"])
											 }) )
								})
	}
}
					
# Right now it's time to visualise this information. To do a barplot, I want to create a table that is the mean and SD of the mut1, mut2 and mut12 values
# For all those strains that are not wild type
plotTable <- array(sapply(names(pairedCAI),function(x){
						return( apply(pairedCAI[[x]],MARGIN=2,function(y){ return( c("Mean"=mean(y),"SD"=sd(y))) }) )
					}),dim=c(2,4,length(names(pairedCAI))),
						dimnames=list(c("Mean","SD"),colnames(pairedCAI[[1]]),names(pairedCAI)))
		
# I now plot this information to be visualised.
pdf(paste(outDir,"paired_CAIvITE.pdf",sep=""),width=6,height=8)
par(mfrow=c(2,1))
plotSets <- c("CAI"="plotTable","ITE"="plotTable_DAMBE")
for(thisTable in names(plotSets)){
	tmpPlot <- barplot(eval(as.name(plotSets[thisTable]))["Mean",,],beside=TRUE,col=c("grey60","red3","royalblue1","purple1"),
						ylab = paste(thisTable,sep=""), xlab = "Pair of Mutations", ylim= c(0,1),
						names.arg = gsub("__:__","\n",names(pairedCAI)),
						legend.text = dimnames(eval(as.name(plotSets[thisTable])))[[2]], args.legend = list(x="topleft",bty="n", cex = 0.7, title="Genotype"), 
						space = c(0,2), cex.names=0.5)
	# Standard deviation bars can be added using the center of bars which gets stored in the tmpPlot object
	suppressWarnings(arrows(tmpPlot,y0=eval(as.name(plotSets[thisTable]))["Mean",,]-eval(as.name(plotSets[thisTable]))["SD",,],
									y1=eval(as.name(plotSets[thisTable]))["Mean",,]+eval(as.name(plotSets[thisTable]))["SD",,],
							lwd=1, code=3,angle=90,length=0.03, col = "grey10"))
}
dev.off()



# Now I can calculate if there is evidence of epistasis between my mutations with respect to the CAI and ITE values calculated for single and paired mutations.
# Here is use all the strains which are mutants as my replicate measures and it is from this that I can then get estimates of error.
epiMeasures <- NULL
# We measure, using error propagation, as per Trindade 2009, for epistasis, and have evidence for epistasis when 
# the multiplicative model for calculating the epsilon value lies outside of error propagation estimate of standard error
#### NOTE: If the data has not been normalised then this is likely to fail due to inter-day variance.
for(thisPair in nonSyn_pairs){
	# sanity check to ensure that our downstream will work
	if(!all(rownames(tableCAI_DAMBE) == rownames(arrayCAI))){ 
		stop("There was a problem with the arrayCAI and tableCAI_DAMBE objects' rownames, please review") 
	}
	# We find the mutations which are part of this pair
	tmp_singMuts <- apply(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")],MARGIN=2,function(x){ strsplit(x,"_")[[1]] })
	tmp_pairString <- paste(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")],collapse="__:__")
	
	# So I can now get the mean and variance "Weights" values for the CAI of our single and double mutants, again double being the geometric mean of the pair
	tmpValues_CAI <- list("WT"= pairedCAI[[tmp_pairString]][,1],
						"mut1"= pairedCAI[[tmp_pairString]][,2],
						"mut2"= pairedCAI[[tmp_pairString]][,3],
						"mut12"=pairedCAI[[tmp_pairString]][,4])
						
	# Now the epsilon and error values will be calculated with the mean and variance of the geometric mean (as per Sharp 1987) of the two positions. 
	tmpEpsilon_CAI <- (mean(tmpValues_CAI[["WT"]]) * mean(tmpValues_CAI[["mut12"]]))  - 
						(mean(tmpValues_CAI[["mut1"]]) * mean(tmpValues_CAI[["mut2"]]))
	# We now calculate the propagated error, knowing that out variance values are returned as zero when there is only one predicted folding.
	tmpError_CAI <- sqrt( (mean(tmpValues_CAI[["WT"]])^2 * var(tmpValues_CAI[["mut12"]])) +
						(mean(tmpValues_CAI[["mut12"]])^2 *var(tmpValues_CAI[["WT"]])) +
						(mean(tmpValues_CAI[["mut1"]])^2 * var(tmpValues_CAI[["mut2"]])) +
						(mean(tmpValues_CAI[["mut2"]])^2 * var(tmpValues_CAI[["mut1"]])) )
	
	# Now we perform much the same for the DAMBE information
	tmpValues_ITE <- list("WT"= frameDAMBE[[tmp_pairString]][,1],
						"mut1"= frameDAMBE[[tmp_pairString]][,2],
						"mut2"= frameDAMBE[[tmp_pairString]][,3],
						"mut12"=frameDAMBE[[tmp_pairString]][,4])
						
	# Now the epsilon and error values will be calculated with the mean and variance of the geometric mean (as per Sharp 1987) of the two positions. 
	tmpEpsilon_ITE <- (mean(tmpValues_ITE[["WT"]]) * mean(tmpValues_ITE[["mut12"]]))  - 
						(mean(tmpValues_ITE[["mut1"]]) * mean(tmpValues_ITE[["mut2"]]))
	# We now calculate the propagated error, knowing that out variance values are returned as zero when there is only one predicted folding.
	tmpError_ITE <- sqrt( (mean(tmpValues_ITE[["WT"]])^2 * var(tmpValues_ITE[["mut12"]])) +
						(mean(tmpValues_ITE[["mut12"]])^2 *var(tmpValues_ITE[["WT"]])) +
						(mean(tmpValues_ITE[["mut1"]])^2 * var(tmpValues_ITE[["mut2"]])) +
						(mean(tmpValues_ITE[["mut2"]])^2 * var(tmpValues_ITE[["mut1"]])) )
	
	# We return values about our pair, giving the absolute epsilon and error estimates
	# as well as testing what confidence we have in them by reporting on how many of the values
	# in our bootstrapped distributions are greater or smaller than our value
	epiMeasures <- rbind(epiMeasures,data.frame("Pair"= tmp_pairString,
												"CAI_Epsilon"= tmpEpsilon_CAI,
												"CAI_estError"= tmpError_CAI,
												"CAI_Epi"=tmpError_CAI < abs(tmpEpsilon_CAI),
												"ITE_Epsilon"= tmpEpsilon_ITE,
												"ITE_estError"= tmpError_ITE,
												"ITE_Epi"=tmpError_CAI < abs(tmpEpsilon_ITE)))		
	
}
print(epiMeasures)
write.csv(epiMeasures,file=paste(outDir,"epiMeasures_CAI_ITE.csv",sep=""), row.names=FALSE)




















# Now there is a concern that (other than the gyrA and parC pair) these mutations are found to be correlated for 
# no reason other than hitchiking.  In a means to suggest this is not the case I'll build the concatenated gene tree
# for my pairs, and also generate an matrix of information which shows the country of strain origin, resistance phenotype state,
# and the nucleotide characters at the sites in question.


# Let's show the number of pairs which exist at each levels of significance
tmpPlot <- unique_sigPairs$pValue_FDR
tmpPlot[which(tmpPlot == 0)] <- min(tmpPlot[which(tmpPlot != 0)])/10
pdf(paste(outDir,"hist_all_pValue_FDR.pdf",sep=""),height=8,width=6)
par(mfcol=c(2,1))
hist(log10(tmpPlot),main="Histogram of all p-Values",xlab="log10 p-Value")
hist(log10(tmpPlot[which(tmpPlot <= 1e-7)]),main="Histogram of small p-Values",xlab="log10 p-Value")
text(x=-11.5,y=25, labels="These points",cex=1, adj=c(0,0.5))
polygon(x = c(-11,-11.5,-11.5,-11), y = c(0,0,nrow(smallP_sigPairs),nrow(smallP_sigPairs)), col=alpha("red",0.5))
dev.off()




# We'll need the wholeAln_index information to run this well
# This last is used by the identification of synonymous and non-synonymous mutations
if(!file.exists(sub("rawAln/","wholeAln_index.RData",seqDir))){
	print("There was a problem with loading the wholeAln_index.RData file, please review")
	q(save="no")
} else {
	load(sub("rawAln/","wholeAln_index.RData",seqDir))
}


# This starts our frame with the country and sensitivity to levofloxacin information, as well we add if that strain bears a mutation at gyrA_T83I
# based on it being like the PA14 character value
aln_gyrA <- read.dna(list.files(path=seqDir,pattern=pa14_genicDB$Locus_Tag[which(pa14_genicDB$geneID == "gyrA")],full.names=TRUE),as.character=TRUE,format="fasta")
reportFrame <- data.frame("Strain"=tTree$tip.label, 
							"isSensitive"=unlist(lapply(tTree$tip.label,function(x){ is.element(x,c("PA01","PA14",susceptStrains)) })),
							"levoMIC"=unlist(lapply(tTree$tip.label,function(x){
															# We extract the levofloxacin MIC values tested
															tmpReturn <- as.character(bacFrame$Levofloxacin)[which(as.character(bacFrame$Isolate) == x)]
															# Now, if there is a ">" symbol, we'll strip that and return 3x the numeric value
															return( if(length(tmpReturn) == 0){
																		NA
																	} else if(grepl(">",tmpReturn)){
																		as.numeric(substr(tmpReturn, 2,nchar(tmpReturn))) * 2
																	# else if the value is ND we return NA
																	} else if (tmpReturn == "ND"){
																		NA
																	# else we simply return the numeric
																	} else {
																		as.numeric(tmpReturn)
																	} ) 
												})),
							"Country"=unlist(lapply(tTree$tip.label,function(x){ 
															tmpReturn <- as.character(bacFrame$Country)[which(as.character(bacFrame$Isolate) == x)] 
															if(length(tmpReturn) == 0){
																"Unknown"
															} else {
																tmpReturn
															}
												})),
							"gyrA_mut"=unlist(lapply(tTree$tip.label,function(x){
															as.vector(aln_gyrA[which(rownames(aln_gyrA) == x),248]) != 
															as.vector(aln_gyrA[which(rownames(aln_gyrA) == "PA14"),248])
												})) )
# Step 1 will be the building of concatenated gene trees
for(thisPair in 1:nrow(smallP_sigPairs)){
	# we find the pair's identity
	tmp_parsedInfo <- strsplit(unlist(as.character(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")])),"_")
	# Now we go and find the sequences for the gene(s) and create a tree, this requires that I find the locus_tag
	tmpLoci <- unique(unlist(lapply(tmp_parsedInfo,function(x){ pa14_genicDB$Locus_Tag[which(pa14_genicDB$geneID == x[1])] } )))
	tmpFiles <- unlist(lapply(tmpLoci,function(x){ list.files(path=seqDir,pattern=x,full.names=FALSE) }))
	# We go and copy the file(s) into out outDir, concatenate them, and then create a tree
	for(thisFile in paste("cp ",seqDir, tmpFiles," ",outDir,sep="")){ system(thisFile) }
	system(paste("cp ",catfasta_loc," ",outDir,sep=""))
	system(paste("cd ",outDir,"; ./catfasta2phyml.pl ",paste(tmpFiles,collapse=" "), "> tmpAln.fas",sep=""))
	# If the tree exists we don't make it again
	if(!file.exists(paste(outDir,paste(unique(lapply(tmp_parsedInfo,function(x){ x[1] })),collapse="_"),".tree",sep=""))){
		system(paste("cd ",outDir,"; ",fastTree_loc," -nt -gtr tmpAln.fas > ",paste(unique(lapply(tmp_parsedInfo,function(x){ x[1] })),collapse="_"),".tree",sep=""))
	}
	# Now this tree will need to be rooted using the outgroup we've used in previous steps, 
	tmp_geneTree <- read.tree(paste(outDir,paste(unique(lapply(tmp_parsedInfo,function(x){ x[1] })),collapse="_"),".tree",sep=""))
	removeSeq <- c("PA_7","AZPAE14941","AZPAE14901","AZPAE15042")
	if (length(intersect(removeSeq, tmp_geneTree $tip.label)) > 0){   
		tmp_geneTree <- root(tmp_geneTree,removeSeq)
		tmp_geneTree <- drop.tip(tmp_geneTree, unlist(sapply(removeSeq,USE.NAMES=FALSE,function(x){ which(grepl(x, tmp_geneTree $tip.label)==TRUE) })) )
	}
	write.tree(tmp_geneTree, file=paste(outDir,paste(unique(lapply(tmp_parsedInfo,function(x){ x[1] })),collapse="_"),".tree",sep=""))
	
	# Now we extract the binary state and nucleotide state information and correct the "order" of that information so it can be
	# cbound to our reportFrame data.
	tmpFrame <- NULL
	for(thisGene in tmpLoci){
		tmpId <- pa14_genicDB$geneID[which(pa14_genicDB$Locus_Tag == thisGene)]
		tmp_ntPos <- unlist(lapply(tmp_parsedInfo,function(x){ if(x[1]== tmpId){as.numeric(x[2])}else{NULL} }))
		# We load in the nucleotide information, but need to reorder it based on the tTree$tip.label object
		tmpGene <- read.dna(paste(outDir,tmpFiles[which(tmpLoci == thisGene)],sep=""),format="fasta",as.character=TRUE)
		# We tidy the tmpGene file by removing any insert positions
		# We will identify any positions in this gene which our index defines as being insert positions, these are removed
		tmpInserts <- which(!is.na(wholeAln_index$Insert_ntPos[which(wholeAln_index$Locus_Tag == thisGene)]))
		# We remove any columns that relate to insert positions
		tmpGene <- tmpGene[,setdiff(1:ncol(tmpGene),tmpInserts)]

		# We gather the nucleotide and then binary state information
		tmp_posIndex <- unlist(lapply(paste(tmpId,tmp_ntPos,sep="_"),function(x){ which(posNames == x) }))
		tmpAdd <- cbind(tmpGene[unlist(lapply(tTree$tip.label,function(x){ which(rownames(tmpGene) == x) })), tmp_ntPos],
						builtMat[corrOrder,unlist(lapply(tmp_posIndex, function(y){ 
												which(unlist(lapply(siteID,function(x){ is.element(y,x) })))
											}))])
		colnames(tmpAdd) <- c(paste("ntInfo", tmpId, tmp_ntPos,sep="_"),
								paste("binaryInfo", tmpId, tmp_ntPos,sep="_"))
		# Now I update my reportFrame object
		reportFrame <- data.frame(reportFrame,tmpAdd) 	
	}
}	

# Right now we've got a data.frame with the information of our nucleotide and binary state changes which occured
# related to the country in which the strain was isolated and the CIP resistance phenotype of the strain.
write.csv(reportFrame,file=paste(outDir,"reportFrame.csv",sep=""),row.names=FALSE)							

PGLSmodels <- list()
stats_sisClade <- list()
# This is to work out the PGLS based method of our state pair affecting MIC, then I'll alos perform the sister comparison
for(thisPair in 1:nrow(smallP_sigPairs)){
	# We define which pair of mutations are on this row
	tmp_parsedInfo <- matrix(unlist(lapply(unlist(as.character(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")])),function(x){
									strsplit(x,"_")[[1]]
								})),ncol=2,byrow=TRUE,dimnames=list(NULL,c("gene","ntPos")))
	# We generate the pairString as a convenience
	tmp_pairString <- paste(apply(tmp_parsedInfo, MARGIN=1,function(x){ paste(x,collapse="_") }),collapse="_")
	
	# We find what pairs exist and for which of the sequences
	tmpPairs <- apply(reportFrame[,apply(tmp_parsedInfo, MARGIN=1,function(x){
										which(grepl(paste("binaryInfo",paste(x,collapse="_"),sep="_"),colnames(reportFrame)))[1]
									})],MARGIN=1,function(x){ paste(x,collapse="") })
	
	# Now we can compute the PGLS for this tree
	tmpData <- data.frame("Strain"=reportFrame$Strain,"levoMIC"=log2(reportFrame$levoMIC),"pairedState"=tmpPairs,"gyrA_mut"=as.factor(reportFrame$gyrA_mut))
	# We now remove all our tmpData which has na values for levoMIC
	tmpData <- tmpData[which(!is.na(tmpData$levoMIC)),]
	# we find the tree file that relates to the phylogenetic tree of the concatenated gene(s) for this pair
	tmpTree <- list.files(path=outDir,paste(paste(unique(tmp_parsedInfo[,1]),collapse="_"),".tree",sep=""))
	tmpTree <- ladderize(read.tree(paste(outDir, tmpTree,sep="")),right=FALSE)
	# We drop the tips of sequences for which we do not have levofloxacin resistance information
	tmpTree <- drop.tip(tmpTree,setdiff(rownames(reportFrame),rownames(tmpData)))
	# Now we fit our basic model with Brownian motion, we can't use the Ornstein-Ullenbeck correlation from experience of errors
	# I also can't use the full model when the pair is the par_c_gyrA pair since the statePair and gyrA_mut predictors are correlated
	tmp_returnList <- list("bm_simple"= summary(gls(levoMIC ~ pairedState, correlation = corBrownian(phy= tmpTree), data=tmpData))[c("tTable","corBeta","sigma","dims")],
										#"ou_simple"= gls(levoMIC ~ pairedState, correlation = corMartins(1,phy= tmpTree), data=tmpData[which(!is.na(tmpData[,"levoMIC"])),]),
										"bm_full"= if(grepl("gyrA", tmp_pairString)){
														NULL
													} else {
														summary(gls(levoMIC ~ pairedState + gyrA_mut, correlation = corBrownian(phy= tmpTree), data=tmpData))[c("tTable","corBeta","sigma","dims")]
													},
										#"ou_full"= gls(levoMIC ~ pairedState + gyrA_mut, correlation = corMartins(1,phy= tmpTree), data=tmpData[which(!is.na(tmpData[,"levoMIC"])),]),	
										"bm_mut"= summary(gls(levoMIC ~ gyrA_mut, correlation = corBrownian(phy= tmpTree), data=tmpData))[c("tTable","corBeta","sigma","dims")])
										#"ou_mut"= gls(levoMIC ~ gyrA_mut, correlation = corMartins(1,phy= tmpTree), data=tmpData[which(!is.na(tmpData[,"levoMIC"])),]))
	for(thisList in 1:length(tmp_returnList)){ 
		# We handle that some objects do not get returned as anything other than NULL
		if(!is.null(tmp_returnList[[thisList]])){
			names(tmp_returnList[[thisList]]) <- c("tTable","Correlation","Residual_R^2","df") 
			# Now we update the df list items
			tmp_returnList[[thisList]][["df"]] <- c("Total"=tmp_returnList[[thisList]][["df"]]$N,
													"Residual"=tmp_returnList[[thisList]][["df"]]$N - tmp_returnList[[thisList]][["df"]]$p)
		}
	}
	
	# We now simply update our information a bit so the names and content are more legible
	PGLSmodels[[tmp_pairString]] <- tmp_returnList
	
	# In order to perform sister clade comparisons we need to first define all pairs of "11" clades and their sister(s) 
	# (should be singular unless multifurcation which is rare).  To achieve this we'll create a version of the gene(s) tree
	# where the tips are the paired states
	tmp_stateTree <- tmpTree
	# We take the paired states from here as it is part of the "removed tips" updated information
	tmp_tip_df_map <- unlist(lapply(tmp_stateTree$tip.label,function(x){ which(tmpData$Strain == x) }))
	tmp_stateTree$tip.label <- as.character(tmpData$pairedState[tmp_tip_df_map])
	# Now we look for all nodes which have all leaves being in the pairedState of "11" (which is now represented by their labels)
	tmp_focalNodes <- lapply(1: (nrow(tmp_stateTree$edge) +1),function(x){
									# If all the leaves of this node are in the pairedState then we want to extrac the sister clade(s) are.
									if(all(as.character(tips(tmp_stateTree, x)) == "11")){
										return(list("node"=x,"leaves"=tips(tmpTree,x)))
									} else {
										return( NULL )
									} })
	
	# Now we reduce this to the nodes relating to unique subsets, this is found by using the original tip labels and looking for unique largest sets
	tmp_nodes_ofInterest <- which(unlist(lapply(tmp_focalNodes,function(x){ !any(is.null(x)) })))
	# So that my set of logcails works downstream, I need to go backward through the largest sets to the smallest
	tmp_nodes_ofInterest <- tmp_nodes_ofInterest[order(sapply(tmp_nodes_ofInterest,function(x){length(tmp_focalNodes[[x]]$leaves) }),decreasing=TRUE)]
	# We initialise with the single largest set
	tmp_cladeHeads <- list(tmp_focalNodes[[tmp_nodes_ofInterest[1]]]$leaves)
	names(tmp_cladeHeads) <- as.character(tmp_nodes_ofInterest[1])
	for(thisNode in tmp_nodes_ofInterest[-1]){
		# We ask if all the tips labels of this node are already present in an element of tmp_cladeHeads and vice versa
		# We also ask if any of the set is part of the clade
		tmp_cladeCompares <- sapply(tmp_cladeHeads, function(x){ c(all(is.element(tmp_focalNodes[[thisNode]]$leaves, x)),
																	any(is.element(tmp_focalNodes[[thisNode]]$leaves, x)),
																	all(is.element(x,tmp_focalNodes[[thisNode]]$leaves))) })
		# If this set shares no tips with any clade, then this is a new set
		if(all(! tmp_cladeCompares[2,])){
			# We only create a new clade provided that this set is larger than one
			if(length(tmp_focalNodes[[thisNode]]$leaves) > 1){
				tmp_cladeHeads[[length(tmp_cladeHeads) + 1]] <- tmp_focalNodes[[thisNode]]$leaves
				names(tmp_cladeHeads)[length(tmp_cladeHeads)] <- as.character(thisNode)
			}
		# If there is an overlap with one other clade this means it is either larger or smaller, and this is found through 
		# the first and last calls asking if the clade is within or greater than the set
		} else {
			if(length(which(tmp_cladeCompares[2,])) > 1){ stop(paste("Problem when trying to find the clade for node ",thisNode,sep="")) }
			# If the first return is true then all the clade is within the 
			tmp_focalClade <- which(tmp_cladeCompares[2,])
			if(tmp_cladeCompares[3, tmp_focalClade]){
				tmp_cladeHeads[[tmp_focalClade]] <- tmp_focalNodes[[tmp_focalClade]]$leaves
				names(tmp_cladeHeads)[tmp_focalClade] <- as.character(tmp_focalClade)
			# Otherwise this means that the set is within the clade, but we check that this is the case...
			} else if (!tmp_cladeCompares[1, tmp_focalClade]){
				stop(paste("Problem when trying to find the clade for node ",thisNode,sep=""))
			}
			
		}
	}
	# We create an object to store the stats for our clades
	cladeStats <- sapply(as.numeric(names(tmp_cladeHeads)), function(thisClade){
	# Now we go through the head nodes of our clades and find the crown node which adjoins the sister clade
		# We now need to look for the node which has depth greater than this one so we first extract
		# the node upstream, which will be that one which is connected to this
		tmp_crownNode <- as.vector(unlist(tmp_stateTree$edge[which(tmp_stateTree$edge[,2] == thisClade),1]))
		# Now as a sanity I'll complain if there is internal multifurcation, meaning this is multiple values
		if(length(tmp_crownNode) != 1){ return(paste("Found more than one crown node for node ", thisClade,sep="")) }
		# Now we find the alternate connection(s) for the crown node
		tmp_sisterClades <- setdiff(as.vector(unlist(tmp_stateTree$edge[which(tmp_stateTree$edge[,1] == tmp_crownNode),2])), thisClade)
		# Ok, now we have the sister clade(s) defined we create a small data.frame of information so we can run stats
		tmp_subData <- data.frame("Strain"= tmp_cladeHeads[[as.character(thisClade)]],
									"levoMIC"=log2(reportFrame$levoMIC[unlist(lapply(tmp_cladeHeads[[as.character(thisClade)]],function(x){ which(reportFrame$Strain == x)}))]),
									"clade"="11",stringsAsFactors=FALSE)
		for(thisSis in tmp_sisterClades){
			# We only add this sister clade if there are more than one sistering tips
			if(length(tips(tmpTree,thisSis)) > 1){
				tmp_subData <- rbind(tmp_subData, data.frame("Strain"= tips(tmpTree, thisSis),
															"levoMIC"=log2(reportFrame$levoMIC[unlist(lapply(tips(tmpTree, thisSis),function(x){ which(reportFrame$Strain == x)}))]),
															"clade"=paste("Sister",thisSis,sep="_"),stringsAsFactors=FALSE))
			}
		}
		# We only proceed if there were any sister clades identified
		if(length(unique(tmp_subData$clade)) > 1){
			rownames(tmp_subData) <- tmp_subData$Strain
			# We update certain columns to be factors
			for(thisCol in colnames(tmp_subData)[which(!grepl("levoMIC",colnames(tmp_subData)))]){ tmp_subData[,thisCol] <- as.factor(tmp_subData[,thisCol]) }
			
			# Great, now we perform a simple t.test on the means of these which is likely not smart... but whatever at this point
			# and then a PGLS to see if clade predicts levoMIC
			tmp_crownTree <- drop.tip(tmpTree,setdiff(tmpTree$tip.label,as.character(tmp_subData$Strain)))
			# For the PGLS analysis we take the p.value for each clade being different from our "11" clade, which will inherit as being the intercept
			# and then we calculate the log likelihood from a chi2 distribuiton with df = number of comparisons
			# However it's possible that there is no variation in levoMIC across the clades, in which case 
			
			if(length(unique(tmp_subData$levoMIC)) == 1){
				return( c("tTest"=1, "PGLS"=1) )
			} else {
				# We extract the p value and then perform the log likelihood equation calculation as suggested by Slowinski 1993
				tmp_PGLS <- summary(gls(levoMIC ~ clade, data = tmp_subData, correlation = corBrownian(phy = tmp_crownTree)))$tTable[,"p-value"]
				
				tmp_PGLS <- pchisq(q = -2*log( sum(tmp_PGLS[which(!grepl("Intercept",names(tmp_PGLS)))]) ),
						df = length(tmp_PGLS[which(!grepl("Intercept",names(tmp_PGLS)))]) )
				
				return( c("tTest"=t.test(levoMIC ~ clade, data = tmp_subData)$p.value,
						"PGLS"= tmp_PGLS
						)
					)
			}
			
		} else {
			# This means we return null information
			return( c("tTest"=NA,
						"PGLS"=NA) )
		}
	})
	colnames(cladeStats) <- names(tmp_cladeHeads)
	stats_sisClade[[length(stats_sisClade) + 1]] <- cladeStats
	names(stats_sisClade)[length(stats_sisClade)] <- tmp_pairString
	
}
save(PGLSmodels,stats_sisClade,file=paste(outDir,"statsModels_sisterComapre.RData",sep=""))




# We now load the rest of the libraries and then continue analysis
library(seqinr)
library(igraph)
#library(ADGofTest) # This allows me to perform the Anderson-Darling test for distribution fits.
library(genetics) # This allows us to get LD, D, D' and r values calculated quickly.
library(plot3D)
library(mlogit) # this allows multinomial regression
library(randomForest)
library(Hmisc) # This allows us to call rcorr which gives correlation coefficient with P value
library(party) # this is an alternative randomForest method which can be plotted
library(colorRamps)
# This initialises and object for performing randomForest calls checking if either Country of Resistance state can 
# predict the state pair
importanceList <- list()
chisqMatrix <- data.frame(matrix(,ncol=3,dimnames=list(NULL,c("mutPair","Country","isSensitive"))))
for(thisPair in 1:nrow(smallP_sigPairs)){
	# We define which pair of mutations are on this row
	tmp_parsedInfo <- matrix(unlist(lapply(unlist(as.character(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")])),function(x){
									strsplit(x,"_")[[1]]
								})),ncol=2,byrow=TRUE,dimnames=list(NULL,c("gene","ntPos")))
	# We generate the pairString as a convenience
	tmp_pairString <- paste(apply(tmp_parsedInfo, MARGIN=1,function(x){ paste(x,collapse="_") }),collapse="_")
	
	# We find what pairs exist and for which of the sequences
	tmpPairs <- apply(reportFrame[,apply(tmp_parsedInfo, MARGIN=1,function(x){
										which(grepl(paste("binaryInfo",paste(x,collapse="_"),sep="_"),colnames(reportFrame)))[1]
									})],MARGIN=1,function(x){ paste(x,collapse="") })
	# Now we record if either Country or isSensitive can predict the statePair
	importanceList[[tmp_pairString]] <- randomForest(statePair ~ Country + isSensitive, data = data.frame("statePair"=as.factor(tmpPairs),
																					"Country"=as.factor(as.character(reportFrame$Country)),
																					"isSensitive"=as.factor(reportFrame$isSensitive)),
																importance=TRUE)$importance
	# This gives us a plottable classification tree which if singleton means that neither sensitivity nor country predicted the state pair
	tmpTrees <- ctree(statePair ~ Country + isSensitive, data = data.frame("statePair"=as.factor(tmpPairs),
																					"Country"=as.factor(as.character(reportFrame$Country)),
																					"isSensitive"=as.factor(reportFrame$isSensitive)))
	pdf(paste(outDir,"classTree_",tmp_pairString,".pdf",sep=""),height=6,width=6)
	plot(tmpTrees, type="simple")
	dev.off()
	
	# We now calculate the correlation of nomial values by doing a chisq.test on the factors
	# I suppress warnings as I know there are class imbalances.... I don't really care
	chisqMatrix <- suppressWarnings(rbind(chisqMatrix, c(tmp_pairString ,
														chisq.test(x=as.factor(tmpPairs),y=reportFrame$Country)$p.value,
														chisq.test(x=as.factor(tmpPairs),y=as.factor(reportFrame$isSensitive))$p.value)))
}
importanceList
chisqMatrix



# These are colouring objects for helping visualise information on our tree.
sensitivity_plotCols <- c("TRUE"="blue","FALSE"="red")
paired_plotCols <- c("00"="blue","11"="red","01"="purple","10"="orange")
country_plotCols <- c(brewer.pal(11,"RdYlBu")[-6],brewer.pal(11,"PiYG")[-6],"black")
names(country_plotCols) <- levels(bacFrame$Country)
ntPairs_plotCols <- c(brewer.pal(11,"RdYlBu")[-c(6:7)],brewer.pal(11,"PiYG")[-c(6:7)],brewer.pal(11,"BrBG")[-c(5:8)])
names(ntPairs_plotCols) <- unlist(lapply(c("A","T","C","G","-"),function(x){ paste(x,c("A","T","C","G","-"),sep="") }))
micCol_plotCols <- colorRamps::blue2red(length(unique(reportFrame$levoMIC[which(!is.na(reportFrame$levoMIC))])))
names(micCol_plotCols) <- as.character(unique(reportFrame$levoMIC[which(!is.na(reportFrame$levoMIC))])[order(unique(reportFrame$levoMIC[which(!is.na(reportFrame$levoMIC))]))])

# This is a plot of all the binary and nt state pairs put together as one plot.
plotDims <- c(nrow(smallP_sigPairs),5)
# This is ascaling adjustment to be applied to the edge branch length so we can position our
# plotted characters within the frame of the image but so they can be reasonably distinguished.
spaceAdjust <- c(1,2,3)/10

pdf(paste(outDir,"all_pairedInfo_trees.pdf",sep=""),height=9*plotDims[1],width=6*plotDims[2])
par(mfrow= plotDims)
for(thisPair in 1:nrow(smallP_sigPairs)){
	# We define which pair of mutations are on this row
	tmp_parsedInfo <- matrix(unlist(lapply(unlist(as.character(smallP_sigPairs[thisPair,c("siteRef_1","siteRef_2")])),function(x){
									strsplit(x,"_")[[1]]
								})),ncol=2,byrow=TRUE,dimnames=list(NULL,c("gene","ntPos")))
	# We generate the pairString as a convenience
	tmp_pairString <- paste(apply(tmp_parsedInfo, MARGIN=1,function(x){ paste(x,collapse="_") }),collapse="_")
	
	
	# We find what pairs exist and for which of the sequences
	tmp_binaryPairs <- apply(reportFrame[,apply(tmp_parsedInfo, MARGIN=1,function(x){
										which(grepl(paste("binaryInfo",paste(x,collapse="_"),sep="_"),colnames(reportFrame)))[1]
									})],MARGIN=1,function(x){ paste(x,collapse="") })
									
	tmp_ntPairs <- toupper(apply(reportFrame[,apply(tmp_parsedInfo, MARGIN=1,function(x){
							which(grepl(paste("ntInfo",paste(x,collapse="_"),sep="_"),colnames(reportFrame)))[1]
						})],MARGIN=1,function(x){ paste(x,collapse="") }))
	
	# we find the tree file(s) for this pair
	tmpTree <- list.files(path=outDir,paste(paste(unique(tmp_parsedInfo[,1]),collapse="_"),".tree",sep=""))
	# We read the tree and then extract the mean edge length
	tmpTree <- ladderize(read.tree(paste(outDir, tmpTree,sep="")),right=FALSE)
	adjFactor <- mean(tmpTree$edge.length) * spaceAdjust
	# This creates a vector of the colours which should be assicated to the tip labels, this is based on
	# this types of pairs being considered and our reference colour vector
	tmp_binCols <- unlist(lapply(tmpTree$tip.label,function(thisTip){ paired_plotCols[tmp_binaryPairs[which(names(tmp_binaryPairs) == thisTip)]] }))
	tmp_ntCols <- unlist(lapply(tmpTree$tip.label,function(thisTip){ ntPairs_plotCols[tmp_ntPairs[which(names(tmp_ntPairs) == thisTip)]] }))
	tmp_countryCols <- unlist(lapply(tmpTree$tip.label,function(thisTip){ 
								country_plotCols[which(names(country_plotCols) == as.character(reportFrame$Country)[which(as.character(reportFrame$Strain) == thisTip)])] 
							}))
	tmp_sensitivityCols <- unlist(lapply(tmpTree$tip.label,function(thisTip){ 
								tmpReturn <- reportFrame$isSensitive[which(as.character(reportFrame$Strain) == thisTip)]
								if(is.na(tmpReturn)){
									"white"
								} else {
									sensitivity_plotCols[which(names(sensitivity_plotCols) == as.character(tmpReturn))] 
								}
							}))
	tmp_propMIC <- unlist(lapply(tmpTree$tip.label,function(thisTip){
								tmpReturn <- reportFrame$levoMIC[which(reportFrame$Strain == thisTip)]
								if(is.na(tmpReturn)){
									"white"
								} else {
									micCol_plotCols[which(names(micCol_plotCols) == as.character(tmpReturn))] 
								}
							}))
								
	# We now strip the actual labels from the tree so they don't tmeselves plot
	tmpTree$tip.label <- rep("",length(tmpTree$tip.label))
	# Now we plot, and then add labels for our information if interest, first the binary state, then nucleotide state, then country, then resistance phenotype
	plot(tmpTree, main= tmp_pairString, edge.width = 2)
	tiplabels(pch=16, col = tmp_binCols,  frame="none",cex=2)
	legend("topleft",unique(names(tmp_binCols)),title = "Binary Pair",col= unique(tmp_binCols),bty="n",cex=3, pch =16)
	
	plot(tmpTree, main= tmp_pairString, edge.width = 2)
	tiplabels(pch=15, col = tmp_ntCols, frame="none",cex=2)
	legend("topleft",unique(names(tmp_ntCols)), title = "Nucleotide Pair",col= unique(tmp_ntCols),bty="n",cex=3,ncol=1, pch = 15)
	
	plot(tmpTree, main= tmp_pairString, edge.width = 2)
	tiplabels(pch=17, col = tmp_countryCols, frame="none",cex=2)
	legend("topleft",unique(names(tmp_countryCols)), title = "Country of Isolation",col= unique(tmp_countryCols),bty="n",cex=1.25,ncol=2, pch = 17)
	
	plot(tmpTree, main= tmp_pairString, edge.width = 2)
	tiplabels(pch=18, col = tmp_sensitivityCols, frame="none",cex=2)
	legend("topleft",unique(names(tmp_sensitivityCols)), title = "Levofloxacin Sensitive",col= unique(tmp_sensitivityCols),bty="n",cex=3,ncol=1, pch = 18)
	#########
	plot(tmpTree, main= tmp_pairString, edge.width = 2)
	tiplabels(pch=19, col = tmp_propMIC , frame="none",cex=2)
	legend("topleft",unique(names(tmp_propMIC)[order(as.numeric(names(tmp_propMIC)))]), title = "Levofloxacin MIC",
			col= unique(tmp_propMIC[order(as.numeric(names(tmp_propMIC)))]),bty="n",cex=2,ncol=2, pch = 19)

}
dev.off()


# Now let's plot one of those classic netwrok diagrams
# Here we visualise networks of mutations, we have not yet made any specification about expected genes...
# We will also now include a plot of the network of correlated evolution that is proposed by our keepPairs$ppInt_vs_sngOverlap
networkCols <- brewer.pal(4, name = "Set1")
# We start by creating the data.frame of info we want to plot
keepPairs_network <- smallP_sigPairs[,c(siteRef_colSet,"pValue_FDR")]													
# Now we adjust our pValue_FDR and "fix" any values which are 0, since adjacency matrix will not handle these normally
# We will set all values of 0 to be 10% smaller than the smalles value, we also handle if the value of NA exists and we ignore it.
if (length(intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))) > 0 && is.element(0,keepPairs_network$pValue_FDR) ){
	min_qVal <- min(keepPairs_network$pValue_FDR[intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))])
	keepPairs_network$pValue_FDR[which(keepPairs_network$pValue_FDR == 0)] <- as.numeric(min_qVal) * 0.1
} else if (all(keepPairs_network$pValue_FDR == 0, na.rm = TRUE)) {
	# This means that our qValues are all equal to zero so we assign some bening small value
	keepPairs_network$pValue_FDR[which(!is.na(keepPairs_network$pValue_FDR))] <- 2e-16
}
# In order to use the iGraph adjacency matrix function with weighted edges we need to create this matrix
# We will extract the siteRef information from rows which have qValues (so are not NA)
adjVertices <- unique(as.vector(sapply(which(grepl("siteRef_",colnames(keepPairs_network))),function(j){
						sapply(which(!is.na(keepPairs_network$pValue_FDR)),function(i){ 
							return( keepPairs_network[i,j] ) 
						}) })))						
adjMatrix <- matrix(0,ncol=length(adjVertices),nrow=length(adjVertices), dimnames=list(adjVertices, adjVertices))
# Now we go through and extract qValue information from rows which are not NA (meaning this data had no instance)
for (thisRow in which(!is.na(keepPairs_network$pValue_FDR))){
	tmpVal <- as.numeric(keepPairs_network$pValue_FDR[thisRow])
	adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]]] <- tmpVal
	adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]]] <- tmpVal
}
aPlot <- graph.adjacency(adjMatrix, mode="undirected", weighted=TRUE)
# In order to view these edges properly we will exxagerate the edge weights
# We use the iGraph function of E()
E(aPlot)$weight <- -log10(E(aPlot)$weight)
pdf(file=paste(outDir,"network_smallP_sigPairs.pdf",sep=""),width=10,height=10)
# This layout method of plotting is from Stephane's example of the iGraph graph.adjacency script see: "iGraph_network_Stephane.R"
# Many of these plotting options are graphical parameters of igraph's plot.igraph

# Here we create a vector to identify how we will colour our vertices, this uses information about if the gene holding a mutation
# is one our our expected genes
vertCols <- sapply(adjVertices,function(thisVertex){
	thisReturn <- NULL
	# Is the siteRef gene information one of the cannonical resistance genes?  Also is this a synonymous or non-synonymous site mutation?
	if (is.element(substr(thisVertex,1,gregexpr("_",thisVertex)[[1]][length(gregexpr("_",thisVertex)[[1]])]-1), testedGenes$expect)){
		if(site_isSyn[thisVertex]){
			return( networkCols[1] )
		} else {
			return( networkCols[2] )
		} 	
	} else {
		if(site_isSyn[thisVertex]){
			return( networkCols[3] )
		} else {
			return( networkCols[4] )
		}
	} })
		
# Now we will handle the edge.widths by creating a binned set of sizes, we will make this set with 6 bins
# Since this will create a vector of factors we use the factor level to represent the widths
thickenLine <- function(inNums){ 3* seq(2,11,length.out = 10)[as.factor(inNums)] }
edgeBin <- cut(E(aPlot)$weight, breaks=5); edge_plotThickness <- thickenLine(edgeBin) 
# NOTE: this is a number of graphical parameters for the plot.igraph option as found through: "http://igraph.org/r/doc/plot.common.html"
#plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= round(E(aPlot)$weight,3), vertex.size = 20, edge.label.cex=3, edge.width=edgeBin, 
#	vertex.label.cex=3, vertex.color= vertCols, asp = 0.66, rescale = .5)
plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= NULL, vertex.size = 25, edge.label.cex=1, edge.width= edge_plotThickness , 
	vertex.label.cex=2, vertex.color= vertCols, asp = 1)

legend("bottomleft",legend=levels(edgeBin),lwd= thickenLine(1:length(edgeBin)) ,title="Edge -log10(qValue)",bty="n",cex=1)
legend("topleft",legend=c("Cannonical_Gene_synMut", "Cannonical_Gene_nonsynMut","Novel_Gene_synMut","Novel_Gene_nonsynMut"),
		fill= networkCols,title="Mutation Type",bty="n",cex=1)
dev.off()
	








 

# We end our script here, anything below is old stuff kept about so it can be quickly accesed

q(save="no")





















































































































############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
################################### BELOW HERE IS OLD CONTENT KEPT FOR QUICK ACCESS ########################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################


# Now if we are limiting ourselves to only pairs with nonSynonymous mutations we act here to subset our work
if(use_nonSyn_only){
	if(!file.exists(paste(postDir,"nonSyn_sigPairs_",ourAccept,".RData",sep=""))){
		# Now we look at all pairs which involve a non-synonymous mutation in a gene of interest
		isPair_nonSyn <- apply(unique_sigPairs[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
									tmpPair <- sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] })
									# We check that both of these are within site_isSyn (which drops certain instances....)
									if(!all(is.element(tmpPair,names(site_isSyn)))){
										return( NA )
									# We ask if any of these are nonSynonymous
									} else if(all(!site_isSyn[tmpPair])){
										return( TRUE )
									} else {
										return( FALSE )
									}
								})
		nonSyn_sigPairs <- unique_sigPairs[which(isPair_nonSyn),]
		save(nonSyn_sigPairs, isPair_nonSyn, file=paste(postDir,"nonSyn_sigPairs_",ourAccept,".RData",sep=""))
	} else { 
		load(paste(postDir,"nonSyn_sigPairs_",ourAccept,".RData",sep=""))
	}
	# We then actually make certain this object is unique_sigPairs only....
	unique_sigPairs <- nonSyn_sigPairs
	duplicateRows <- find_duplicateRows(unique_sigPairs,"siteRef_1","siteRef_2")
	removeRows <- if (length(duplicateRows) == 0){ 
						NULL 
					} else { 
						unlist(sapply(1:length(duplicateRows),function(x){
							# We look if this list entry is NULL, if not we need to review which to keep
							if (length(duplicateRows[[x]]) > 0){
								tmpSet <- c(as.numeric(names(duplicateRows)[x]),duplicateRows[[x]])
								# now we return the non-minimum pValue_FDR rows as reported by sigPairs row x or those duplicates.
								return( tmpSet[-which(unique_sigPairs$pValue_FDR[tmpSet] == min(unique_sigPairs$pValue_FDR[tmpSet]))[1]] )
							} else {
								return( NULL )
							} })) 
					}
	assign("unique_sigPairs", if (length(removeRows) > 0){ unique_sigPairs[-removeRows,] } else { unique_sigPairs }, pos=".GlobalEnv")
}


if(!file.exists(paste(postDir,"nonSyn_expectPairs_minInt_",min_intScore,"_ourAccept",ourAccept,".RData",sep=""))){
	# Now we look at all pairs which involve a non-synonymous mutation in a gene of interest
	nonSyn_expecCheck <- apply(unique_sigPairs[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
								tmp_pairInfo <- sapply(sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] }), function(y){ strsplit(y,"_")[[1]] })
								# If none of these are in expected genes we can stop here
								tmp_isExpect <- is.element(tmp_pairInfo[1,], testedGenes$expect)
								if(!any(tmp_isExpect)){ 
									return( FALSE )
								# Now we check if either of these are nonSynonymous sites
								} else {
									tmpPair <- apply(tmp_pairInfo, MARGIN = 2,function(x){ paste(x,collapse="_") })
									# We check that both of these are within site_isSyn (which drops certain instances....)
									if(!all(is.element(tmpPair,names(site_isSyn)))){
										return( NA )
									}
									# We ask if any of these are nonSynonymous
									if(any(!site_isSyn[tmpPair])){
										return( TRUE )
									} else {
										return( FALSE )
									}
								}
							})
	# This simply checks if pairs are non synonymous, nothing else						
	nonSyn_check <- apply(unique_sigPairs[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
		tmp_pairInfo <- sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] }) 
		return(all(!site_isSyn[tmp_pairInfo]))
		})
	nonSyn_both_unique_sigPairs <- unique_sigPairs[which(nonSyn_check),]
	# We can now subset our data to the pairs which have at least one non-synonymous mutation in an expected gene
	nonSyn_expectPairs <- unique_sigPairs[which(nonSyn_expecCheck),]
	nonSyn_expectPairs <- nonSyn_expectPairs[order(nonSyn_expectPairs$pValue_FDR,decreasing = FALSE),]
	# We can now subset our data further to the pairs which have both non-synonymous mutations in an expected gene
	nonSyn_expectBoth <- nonSyn_expectPairs[which(apply(nonSyn_expectPairs[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
								tmp_pairInfo <- sapply(sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] }), function(y){ strsplit(y,"_")[[1]] })
								# If none of these are in expected genes we can stop here
								tmp_isExpect <- is.element(tmp_pairInfo[1,], testedGenes$expect)
								if(!all(tmp_isExpect)){ 
									return( FALSE )
								# Now we check if either of these are nonSynonymous sites
								} else {
									tmpPair <- apply(tmp_pairInfo, MARGIN = 2,function(x){ paste(x,collapse="_") })
									# We check that both of these are within site_isSyn (which drops certain instances....)
									if(!all(is.element(tmpPair,names(site_isSyn)))){
										return( NA )
									}
									# We ask if any of these are nonSynonymous
									if(all(!site_isSyn[tmpPair])){
										return( TRUE )
									} else {
										return( FALSE )
									}
								}
							})), ]
	# I remove duplicate rows
	# I now run a re-setting of the unique_sigPairs in case there are duplicates the crept in from this step.... (found in past)
	# Now we define the unique significant pairs, first we get a list, then for each list entry we find which row has the smallest pValue_FDR
	duplicateRows <- find_duplicateRows(nonSyn_expectPairs,"siteRef_1","siteRef_2")
	removeRows <- if (length(duplicateRows) == 0){ 
						NULL 
					} else { 
						unlist(sapply(1:length(duplicateRows),function(x){
							# We look if this list entry is NULL, if not we need to review which to keep
							if (length(duplicateRows[[x]]) > 0){
								tmpSet <- c(x,duplicateRows[[x]])
								# now we return the non-minimum pValue_FDR rows as reported by sigPairs row x or those duplicates.
								return( tmpSet[-which(nonSyn_expectPairs$pValue_FDR[tmpSet] == min(nonSyn_expectPairs$pValue_FDR[tmpSet]))[1]] )
							} else {
								return( NULL )
							} })) 
					}
	nonSyn_expectPairs <- if (length(removeRows) > 0){ nonSyn_expectPairs[-removeRows,] } else { nonSyn_expectPairs }

	# Now I go about removing any intragenic pairs with a significant LD
	tmpData <- vector(mode="list",length=nrow(nonSyn_expectPairs))
	for(x in 1:nrow(nonSyn_expectPairs)){
		thisPair <- unlist(nonSyn_expectPairs[x, siteRef_colSet])
		# In case we have sites which had non-unqijue site patterns we simply consider the first instance in this assessment
		thisPair <- sapply(thisPair, function(x){ strsplit(x,"::")[[1]][1] })
		# We extract the sitePattern for these two positions
		tmp_sitePatterns <- cbind(traitMat[,which(posNames == thisPair[1])],traitMat[,which(posNames == thisPair[2])])
		# We remove any gap characters, which are stored as 2.
		tmp_sitePatterns <- tmp_sitePatterns[intersect(which(tmp_sitePatterns[,1] != 2),which(tmp_sitePatterns[,2] != 2)),]
		# We now find our which of the site pairs is sensitive and resistant
		# This will create genotype objects which are the "changes" based on the root state 
		tmpData[[x]] <- LD(genotype(paste(tmp_sitePatterns[,1],tmp_sitePatterns[,1],sep=""),sep=""),genotype(paste(tmp_sitePatterns[,2],tmp_sitePatterns[,2],sep=""),sep=""))
	}
	# Now I get the distance data for my pairs, since this is a circular genome we must know the maximum length of the reference genome 
	# and look at the smallest value of the linear distance between them, and the sum of the smallest position and the diff of the largest and genome length
	# This will require us to quickly look at the length of the reference genome
	tmpLen_refGen <- ncol(read.dna(paste(baseDir,"/PARuns/Sequences/PA_wholeGenomes/allGenomes/PA14.fas",sep=""),format="fasta"))
	tmpDist <- as.numeric(apply(nonSyn_expectPairs[, siteRef_colSet], MARGIN = 1, function(thisPair){
						tmpGenes <- sapply(sapply(thisPair,function(y){ strsplit(y,"::")[[1]][1]}), function(x){ strsplit(x,"_")[[1]] })
						tmpPositions <- apply(tmpGenes, MARGIN = 2, function(x){ pa14_genicDB$Start[which(pa14_genicDB$geneID == x[1])[1]] - 1 + as.numeric(x[2]) })
						# Now that we have the positions, in PA14 reference genome of each of these mutations, we find the smallest of either
						# their absolute difference, OR the sum of the smallest and the difference between the largest and the length of the genome
						return( min(abs(diff(tmpPositions)),sum(min(tmpPositions), tmpLen_refGen - max(tmpPositions) + 1)) ) }))
	tmpReview <- cbind(nonSyn_expectPairs[, c(siteRef_colSet,"pValue_FDR")], tmpDist, sapply(tmpData,function(x){ return(x[["P-value"]])}))
	colnames(tmpReview)[4:5] <- c("Dist","LD_pValue")
	
	is_intraGenic <- apply(nonSyn_expectPairs[, siteRef_colSet], MARGIN = 1, function(x){ strsplit(x[1],"_")[[1]][1] == strsplit(x[2],"_")[[1]][1] })
		
	pdf(paste(postDir,"plots_distLD_nonSyn_expectPairs.pdf",sep=""), height = 12, width = 12)
	par(mfrow=c(2,2))
	plot(tmpReview[,"LD_pValue"],nonSyn_expectPairs$pValue_FDR, xlab = "LD P-Value", ylab = "Corr. Evol. qValue", log="y")
	lines(density(tmpReview[,"LD_pValue"]),col="red")
	plot(tmpDist,nonSyn_expectPairs$pValue_FDR, xlab = "Physical Distance (nt positions)", ylab = "Corr. Evol. qValue",log = "y")
	plot(tmpDist[which(is_intraGenic)],tmpReview[which(is_intraGenic),"LD_pValue"],ylab = "LD P-Value", xlab = "Physical Distance (nt positions)")
	#plot(tmpDist[which(is_intraGenic)],nonSyn_expectPairs$pValue_FDR[which(is_intraGenic)], xlab = "Physical Distance (nt positions)", ylab="Corr. Evol. qValue",log="y")
	# We now visualise our data plotted in 3D space, we colour based on gene types being in the set.
	tmp_qVals <- -log10(nonSyn_expectPairs$pValue_FDR); tmpMax <- max(tmp_qVals[which(tmp_qVals != "Inf")])*1.2
	for(thisVal in which(tmp_qVals == Inf)){ tmp_qVals[thisVal] <- tmpMax }
	scatter3D(tmpDist, tmpReview[,"LD_pValue"], tmp_qVals, xlab="Physical Distance (nt positions)", 
				ylab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 300, phi = 30, pch = 19, cex.lab = 1)
	dev.off()
	
	# This is a scatterplot of the points based on the LD values, signal of corr evol, and the maximum likelihood probability of the root state
	# We use the siteRef_# information to reference the unique sig pairs data frame and find the maximum likelihood root state recorded.
	# We do this for both the dependent and independent root states
	tmpReview <- cbind(tmpReview, t(apply(tmpReview[,siteRef_colSet], MARGIN = 1, function(x){ 
						tmpRow <- intersect(which(nonSyn_expectPairs$siteRef_1 == x[1]),which(nonSyn_expectPairs$siteRef_2 == x[2]))
						return( c("mlIndep" = max(nonSyn_expectPairs[tmpRow,intersect(which(grepl("rootP_",colnames(nonSyn_expectPairs))),which(grepl("_indep",colnames(nonSyn_expectPairs))))]), 
									"mlDep" = max(nonSyn_expectPairs[tmpRow,intersect(which(grepl("rootP_",colnames(nonSyn_expectPairs))),which(grepl("_dep",colnames(nonSyn_expectPairs))))])) ) 
						})))
	pdf(paste(postDir,"plots_rootPvsLD_nonSyn_expectPairs.pdf",sep=""), height = 6, width = 18)
	par(mfrow=c(1,3))
	tmp_qVals <- -log10(nonSyn_expectPairs$pValue_FDR); tmpMax <- max(tmp_qVals[which(tmp_qVals != "Inf")])*1.2
	for(thisVal in which(tmp_qVals == Inf)){ tmp_qVals[thisVal] <- tmpMax }
	scatter3D(tmpReview[,"LD_pValue"], tmpReview[,"mlIndep"], tmp_qVals, ylab="Root State ML (indep model)", 
				xlab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)
	scatter3D(tmpReview[,"LD_pValue"], tmpReview[,"mlDep"], tmp_qVals, ylab="Root State ML (dep model)", 
				xlab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)
	# Here we plot the colours as the qValue
	scatter3D(tmpReview[,"mlIndep"], tmpReview[,"mlDep"], tmpReview[,"LD_pValue"], ylab="Root State ML (dep model)", 
				xlab = "Root State ML (indep model)", zlab = "LD P-Value", ticktype = "detailed", colvar = tmp_qVals,
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)			
	dev.off()
	
	# From having a look at this we will filter by removing those sites which show significant LD AND are intra-genic
	nonSyn_expectPairs_intraLD <- nonSyn_expectPairs[union(which(!is_intraGenic),intersect(which(is_intraGenic),which(tmpReview[,"LD_pValue"] > 0.05))),]
	# We will now extract all the pairs in which both mutations are non-synonymous
	nonSyn_bothExpect_intraLD <- nonSyn_expectPairs_intraLD[which(apply(nonSyn_expectPairs_intraLD[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
								tmp_pairInfo <- sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] })
								# We ask if both pairs are non-synonymous mutations or not.
								if(!all(is.element(tmp_pairInfo,names(site_isSyn)))){
									return( NA )
								}
								# We ask if any of these are nonSynonymous
								if(all(!site_isSyn[tmp_pairInfo])){
									return( TRUE )
								} else {
									return( FALSE )
								}
							})), ]
	# Now I go about removing any intragenic pairs with a significant LD
	tmpData <- vector(mode="list",length=nrow(nonSyn_bothExpect_intraLD))
	for(x in 1:nrow(nonSyn_bothExpect_intraLD)){
		thisPair <- unlist(nonSyn_bothExpect_intraLD[x, siteRef_colSet])
		# In case we have sites which had non-unqijue site patterns we simply consider the first instance in this assessment
		thisPair <- sapply(thisPair, function(x){ strsplit(x,"::")[[1]][1] })
		# We extract the sitePattern for these two positions
		tmp_sitePatterns <- cbind(traitMat[,which(posNames == thisPair[1])],traitMat[,which(posNames == thisPair[2])])
		# We remove any gap characters, which are stored as 2.
		tmp_sitePatterns <- tmp_sitePatterns[intersect(which(tmp_sitePatterns[,1] != 2),which(tmp_sitePatterns[,2] != 2)),]
		# We now find our which of the site pairs is sensitive and resistant
		# This will create genotype objects which are the "changes" based on the root state 
		tmpData[[x]] <- LD(genotype(paste(tmp_sitePatterns[,1],tmp_sitePatterns[,1],sep=""),sep=""),genotype(paste(tmp_sitePatterns[,2],tmp_sitePatterns[,2],sep=""),sep=""))
	}
	# Now I get the distance data for my pairs, since this is a circular genome we must know the maximum length of the reference genome 
	# and look at the smallest value of the linear distance between them, and the sum of the smallest position and the diff of the largest and genome length
	# This will require us to quickly look at the length of the reference genome
	tmpLen_refGen <- ncol(read.dna(paste(baseDir,"/PARuns/Sequences/PA_wholeGenomes/allGenomes/PA14.fas",sep=""),format="fasta"))
	tmpDist <- as.numeric(apply(nonSyn_bothExpect_intraLD[, siteRef_colSet], MARGIN = 1, function(thisPair){
						tmpGenes <- sapply(sapply(thisPair,function(y){ strsplit(y,"::")[[1]][1]}), function(x){ strsplit(x,"_")[[1]] })
						tmpPositions <- apply(tmpGenes, MARGIN = 2, function(x){ pa14_genicDB$Start[which(pa14_genicDB$geneID == x[1])[1]] - 1 + as.numeric(x[2]) })
						# Now that we have the positions, in PA14 reference genome of each of these mutations, we find the smallest of either
						# their absolute difference, OR the sum of the smallest and the difference between the largest and the length of the genome
						return( min(abs(diff(tmpPositions)),sum(min(tmpPositions), tmpLen_refGen - max(tmpPositions) + 1)) ) }))
	tmpReview <- cbind(nonSyn_bothExpect_intraLD[, c(siteRef_colSet,"pValue_FDR")], tmpDist, sapply(tmpData,function(x){ return(x[["P-value"]])}))
	colnames(tmpReview)[4:5] <- c("Dist","LD_pValue")
	
	is_intraGenic <- apply(nonSyn_bothExpect_intraLD[, siteRef_colSet], MARGIN = 1, function(x){ strsplit(x[1],"_")[[1]][1] == strsplit(x[2],"_")[[1]][1] })
		
	pdf(paste(postDir,"plots_distLD_nonSyn_bothExpect_intraLD.pdf",sep=""), height = 12, width = 12)
	par(mfrow=c(2,3))
	plot(tmpReview[,"LD_pValue"],nonSyn_bothExpect_intraLD$pValue_FDR, xlab = "LD P-Value", ylab = "Corr. Evol. qValue", log="y")
	lines(density(tmpReview[,"LD_pValue"]),col="red")
	plot(tmpDist,nonSyn_bothExpect_intraLD$pValue_FDR, xlab = "Physical Distance (nt positions)", ylab = "Corr. Evol. qValue",log = "y")
	# We now visualise our data plotted in 3D space, we colour based on gene types being in the set.
	tmp_qVals <- -log10(nonSyn_bothExpect_intraLD$pValue_FDR); tmpMax <- max(tmp_qVals[which(tmp_qVals != "Inf")])*1.2
	for(thisVal in which(tmp_qVals == Inf)){ tmp_qVals[thisVal] <- tmpMax }
	scatter3D(tmpDist, tmpReview[,"LD_pValue"], tmp_qVals, xlab="Physical Distance (nt positions)", 
				ylab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 300, phi = 30, pch = 19, cex.lab = 1)
	# This is a scatterplot of the points based on the LD values, signal of corr evol, and the maximum likelihood probability of the root state
	# We use the siteRef_# information to reference the unique sig pairs data frame and find the maximum likelihood root state recorded.
	# We do this for both the dependent and independent root states
	tmpReview <- cbind(tmpReview, t(apply(tmpReview[,siteRef_colSet], MARGIN = 1, function(x){ 
						tmpRow <- intersect(which(nonSyn_bothExpect_intraLD$siteRef_1 == x[1]),which(nonSyn_bothExpect_intraLD$siteRef_2 == x[2]))
						return( c("mlIndep" = max(nonSyn_bothExpect_intraLD[tmpRow,intersect(which(grepl("rootP_",colnames(nonSyn_bothExpect_intraLD))),which(grepl("_indep",colnames(nonSyn_bothExpect_intraLD))))]), 
									"mlDep" = max(nonSyn_bothExpect_intraLD[tmpRow,intersect(which(grepl("rootP_",colnames(nonSyn_bothExpect_intraLD))),which(grepl("_dep",colnames(nonSyn_bothExpect_intraLD))))])) ) 
						})))
	tmp_qVals <- -log10(nonSyn_bothExpect_intraLD$pValue_FDR); tmpMax <- max(tmp_qVals[which(tmp_qVals != "Inf")])*1.2
	for(thisVal in which(tmp_qVals == Inf)){ tmp_qVals[thisVal] <- tmpMax }
	scatter3D(tmpReview[,"LD_pValue"], tmpReview[,"mlIndep"], tmp_qVals, ylab="Root State ML (indep model)", 
				xlab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)
	scatter3D(tmpReview[,"LD_pValue"], tmpReview[,"mlDep"], tmp_qVals, ylab="Root State ML (dep model)", 
				xlab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)
	# Here we plot the colours as the qValue
	scatter3D(tmpReview[,"mlIndep"], tmpReview[,"mlDep"], tmpReview[,"LD_pValue"], ylab="Root State ML (dep model)", 
				xlab = "Root State ML (indep model)", zlab = "LD P-Value", ticktype = "detailed", colvar = tmp_qVals,
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)			
	dev.off()
	# We now arrange our nonSynonymous both pairs, and we remove instances when there are non-unique site-patterns
	nonSyn_bothExpect_intraLD <- nonSyn_bothExpect_intraLD[order(nonSyn_bothExpect_intraLD$pValue_FDR),]
	nonSyn_bothExpect_intraLD_uniquePattern <- nonSyn_bothExpect_intraLD[intersect(which(nonSyn_bothExpect_intraLD$len_siteID_1 <= 1),which(nonSyn_bothExpect_intraLD$len_siteID_2 <= 1)),]
	
	# will be to review what is the breakdown of the number for unique_sigPairs with decreasing log10 pValue_FDR
	nonSyn_sigLevels <- list()
	tmpAccept <- ourAccept
	while(length(which(nonSyn_bothExpect_intraLD_uniquePattern$pValue_FDR <= tmpAccept)) > 0 && !all(nonSyn_bothExpect_intraLD_uniquePattern$pValue_FDR[which(nonSyn_bothExpect_intraLD_uniquePattern$pValue_FDR <= tmpAccept)] == 0) ){
		nonSyn_sigLevels[[as.character(tmpAccept)]] <- length(which(nonSyn_bothExpect_intraLD_uniquePattern$pValue_FDR <= tmpAccept))
		tmpAccept <- tmpAccept/10
	}
	
	# Here we visualise networks of mutations, we have not yet made any specification about expected genes...
	# We will also now include a plot of the network of correlated evolution that is proposed by our keepPairs$ppInt_vs_sngOverlap
	networkCols <- brewer.pal(4, name = "Set1")
	# We start by creating the data.frame of info we want to plot
	keepPairs_network <- nonSyn_bothExpect_intraLD_uniquePattern[,c(siteRef_colSet,"pValue_FDR")]													
	# Now we adjust our pValue_FDR and "fix" any values which are 0, since adjacency matrix will not handle these normally
	# We will set all values of 0 to be 10% smaller than the smalles value, we also handle if the value of NA exists and we ignore it.
	if (length(intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))) > 0 && is.element(0,keepPairs_network$pValue_FDR) ){
		min_qVal <- min(keepPairs_network$pValue_FDR[intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))])
		keepPairs_network$pValue_FDR[which(keepPairs_network$pValue_FDR == 0)] <- as.numeric(min_qVal) * 0.1
	} else if (all(keepPairs_network$pValue_FDR == 0, na.rm = TRUE)) {
		# This means that our qValues are all equal to zero so we assign some bening small value
		keepPairs_network$pValue_FDR[which(!is.na(keepPairs_network$pValue_FDR))] <- 2e-16
	}
	# In order to use the iGraph adjacency matrix function with weighted edges we need to create this matrix
	# We will extract the siteRef information from rows which have qValues (so are not NA)
	adjVertices <- unique(as.vector(sapply(which(grepl("siteRef_",colnames(keepPairs_network))),function(j){
							sapply(which(!is.na(keepPairs_network$pValue_FDR)),function(i){ 
								return( keepPairs_network[i,j] ) 
							}) })))						
	adjMatrix <- matrix(0,ncol=length(adjVertices),nrow=length(adjVertices), dimnames=list(adjVertices, adjVertices))
	# Now we go through and extract qValue information from rows which are not NA (meaning this data had no instance)
	for (thisRow in which(!is.na(keepPairs_network$pValue_FDR))){
		tmpVal <- as.numeric(keepPairs_network$pValue_FDR[thisRow])
		adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]]] <- tmpVal
		adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]]] <- tmpVal
	}
	aPlot <- graph.adjacency(adjMatrix, mode="undirected", weighted=TRUE)
	# In order to view these edges properly we will exxagerate the edge weights
	# We use the iGraph function of E()
	E(aPlot)$weight <- -log10(E(aPlot)$weight)
	pdf(file=paste(postDir,"nonSyn_bothExpect_intraLD_uniquePattern_min_intScore600_network.pdf",sep=""),width=120,height=120)
	# This layout method of plotting is from Stephane's example of the iGraph graph.adjacency script see: "iGraph_network_Stephane.R"
	# Many of these plotting options are graphical parameters of igraph's plot.igraph
	
	# Here we create a vector to identify how we will colour our vertices, this uses information about if the gene holding a mutation
	# is one our our expected genes
	vertCols <- sapply(adjVertices,function(thisVertex){
		thisReturn <- NULL
		# Is the siteRef gene information one of the cannonical resistance genes?  Also is this a synonymous or non-synonymous site mutation?
		if (is.element(substr(thisVertex,1,gregexpr("_",thisVertex)[[1]][length(gregexpr("_",thisVertex)[[1]])]-1), testedGenes$expect)){
			if(site_isSyn[thisVertex]){
				return( networkCols[1] )
			} else {
				return( networkCols[2] )
			} 	
		} else {
			if(site_isSyn[thisVertex]){
				return( networkCols[3] )
			} else {
				return( networkCols[4] )
			}
		} })
			
	# Now we will handle the edge.widths by creating a binned set of sizes, we will make this set with 6 bins
	# Since this will create a vector of factors we use the factor level to represent the widths
	thickenLine <- function(inNums){ 3* seq(2,11,length.out = 10)[as.factor(inNums)] }
	edgeBin <- cut(E(aPlot)$weight, breaks=10); edge_plotThickness <- thickenLine(edgeBin) 
	# NOTE: this is a number of graphical parameters for the plot.igraph option as found through: "http://igraph.org/r/doc/plot.common.html"
	#plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= round(E(aPlot)$weight,3), vertex.size = 20, edge.label.cex=3, edge.width=edgeBin, 
	#	vertex.label.cex=3, vertex.color= vertCols, asp = 0.66, rescale = .5)
	plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= NULL, vertex.size = 5, edge.label.cex=1, edge.width= edge_plotThickness , 
		vertex.label.cex=3, vertex.color= vertCols, asp = 1)
	
	legend("bottomleft",legend=levels(edgeBin),lwd= thickenLine(1:length(edgeBin)) ,title="Edge -log10(qValue)",bty="n",cex=10)
	legend("topleft",legend=c("Cannonical_Gene_synMut", "Cannonical_Gene_nonsynMut","Novel_Gene_synMut","Novel_Gene_nonsynMut"),
			fill= networkCols,title="Mutation Type",bty="n",cex=10)
	dev.off()
	
	# To be noted: From my working dataset there is only one pair of non-synonymous mutations both in expected genes... guess which (or just review).
	save(nonSyn_both_unique_sigPairs, nonSyn_expectPairs_intraLD, nonSyn_expectPairs, nonSyn_expectBoth, nonSyn_bothExpect_intraLD, 
			nonSyn_bothExpect_intraLD_uniquePattern, file = paste(postDir,"nonSyn_allPairs_minInt_",min_intScore,"_ourAccept",ourAccept,".RData",sep=""))
} else {
	load(paste(postDir,"nonSyn_expectPairs_minInt_",min_intScore,"_ourAccept",ourAccept,".RData",sep=""))
}


# This is a first pass qualitative assessment of whether our analytical method found more significantly correlated pairs
# when at least one site is within our six expected genes, the alternate being any other pairs.  We will break down our groups 
# into those pairs with, 0,1,2 sites with an expect gene.  We need ALL our data for this....
if(!file.exists(paste(postDir, fileName_expGene_numSigs,sep=""))){
	if(file.exists(paste(postDir,"compiledRun_data.RData",sep=""))){
		# We now load our previously made object which compiled all our runs of BayesTraits.
		load(paste(postDir,"compiledRun_data.RData",sep=""))
		# This is a matrix which tells us if a siteRef_<column> contains one of the expected genes
		expectGene_refMatrix <- cbind(apply(sapply(testedGenes$expect,function(x){ grepl(x, compiledRun_data$siteRef_1) }), MARGIN = 1, any),
										apply(sapply(testedGenes$expect,function(x){ grepl(x, compiledRun_data$siteRef_2) }), MARGIN = 1, any))
		colnames(expectGene_refMatrix) <- c("siteRef_1","siteRef_2")
		# We use the fact that logicals can be summed 
		tmpSum <- apply(expectGene_refMatrix, MARGIN = 1, FUN = sum)
		pairs_byType <- list("none"=which(tmpSum == 0),
						"one"=which(tmpSum == 1),
						"two"=which(tmpSum == 2)) 
		# Let's look at the varying levels of significance we could tolerate and see if at those levels we find alternate 
		# results for our analysis of correlated evolution showing more significant pairs with our expected genes.
		analysisSig <- list()
		for(tmpAccept in as.numeric(names(sigLevels))){
			# Now the goal is to see how many of these are significant or not
			sig_byType <- sapply(pairs_byType, function(x){ return( c(length(which(compiledRun_data$pValue_FDR[x] <= tmpAccept)), 
																	length(which(compiledRun_data$pValue_FDR[x] > tmpAccept))) ) })
			rownames(sig_byType) <- c("signif","notSignif")
			# Here we review the proportion of each type of analyses that are significant
			analysisSig[[as.character(tmpAccept)]] <- list("sig_byType"=sig_byType, "proportions"=apply(sig_byType,MARGIN = 2, function(x){ return(x[1]/sum(x)) }),
					"anySites"=poisson.test(sum(sig_byType["signif",c("one","two")]), r = sig_byType["signif","none"]/sum(sig_byType[,"none"]), T = sum(sig_byType[,c("one","two")]), alternative = "greater"),
					"dblSites"=poisson.test(sum(sig_byType["signif","two"]), r = sig_byType["signif","none"]/sum(sig_byType[,"none"]), T = sum(sig_byType[,"two"]), alternative = "greater"))
		}
		# We will now instead of looking at expected genes, simply look for intra vs inter geneic pairs.
		# We do have an expectation that intra should be higher.
		is_intraGenic <- apply(compiledRun_data[,c("siteRef_1","siteRef_2")], MARGIN = 1, function(x){ strsplit(x[1],"_")[[1]][1] == strsplit(x[2],"_")[[1]][1] })
		intraGenic_sigList <- list()
		for(tmpAccept in as.numeric(names(sigLevels))){
			sig_byIntra <- sapply(list("intra"=which(is_intraGenic),"inter"=which(!is_intraGenic)),function(x){
									return( c(length(which(compiledRun_data$pValue_FDR[x] <= tmpAccept)), 
												length(which(compiledRun_data$pValue_FDR[x] > tmpAccept))) ) })
			rownames(sig_byIntra) <- c("signif","notSignif")
			intraGenic_sigList[[as.character(tmpAccept)]] <- list("sig_byIntra"= sig_byIntra,
					"intraSig"=poisson.test(sig_byIntra["signif","intra"], r = sig_byIntra["signif","inter"]/sum(sig_byIntra[,"inter"]), T = sum(sig_byIntra[,"intra"]), alternative = "greater"))
		}	
		# Now we save our results												
		save(analysisSig, intraGenic_sigList, file = paste(postDir, fileName_expGene_numSigs,sep=""))
		rm(compiledRun_data)
	} else {
		print("Could not find the compiledRun_data.RData file, did not do this part")
	}
} else {
	load(paste(postDir, fileName_expGene_numSigs,sep=""))
}
## FROM ACTUALLY PERFORMING THIS WE FIND THAT OUR ANALYTICAL METHOD IS NOT STRICTLY BETTER AT FINDING CORRELATED EVOLUTION AMONG
## EXPECTED GENES, IT DOES FIND MORE INTRA GENIC CORRELATED EVOLUTION WHICH IS EXPECTED FROM THEORY ON EPISTASIS, BUT
## IS ALSO MORE PROBABLE IF HITCHIKING IS A FACTOR, OR OTHER BIOLOGICAL REASONS.

# We will now create a condensed object that is the list of unique site interactions we will be considering
unique_sigPairs_list <- extractPair(unique_sigPairs, siteRef_colSet)

interactingPairs <- NULL
############## THIS IS FOR THE STRING DATABASE INTERACTIONS INFORMATION ##############
# Here we will use the gene interactions information to help visualise which of our protein interactions are supported.
# These are the new pa01_ & pa14_genicDB and geneInteractions & _links.detailed objects that I can use.... NOTE: for GI information I need
# to work around a bit better as the GI of PA01 != GI: of PA14
# CHALLENGE IS THAT THE ppInt INFROMATION IS FOR THE REFERENCE STRAIN PA01 BUT ALL MY ANALYSIS HAS USED THE PA14 REFERENCE
# THIS MEANS WE WILL HAVE TO USE BLAST TO TRY AND FIND THE SEQUENCE ANALOGS BETWEEN GI: SEQUENCES.
if(file.exists(paste(postDir,"interacting_uniquePairs_minInt",min_intScore,".RData",sep=""))){	
	 load(paste(postDir,"interacting_uniquePairs_minInt",min_intScore,".RData",sep=""))
} else {
	# We will extract those of our proposed pairs we want to explore and look if there are ppInt that exist
	# We will use all unique sngOverlap, bldOverlap and dblSupport sites
	for(thisPair in unique_sigPairs_list){
		# We only really need the gene information not the actual site of the mutation.
		tmpPair <- unlist(lapply(strsplit(thisPair,pairString)[[1]],function(x){ substr(x,1,regexpr("_",x)-1) }))
		this_genePair <- paste(tmpPair,collapse=pairString)
		# First we verify if this_genePair is known to not interact
		if(is.element(this_genePair, known_notInteractions)){
			# In this case there is simply nothing to be done we move onto the next.
			next
		# Then we verify if this pair exists in our knownInteractions
		} else if(is.element(this_genePair,knownInteractions$genePair_ID)){
			# This means the pair is part of our knownInterctions and thus we will add it to these interacting pairs
			interactingPairs <- c(interactingPairs, which(unique_sigPairs_list == thisPair))
		} else {
			# This means that this_genePair is neither known to interact or not, thus we verify it.
			# Now we look to see if there is a previously identified reference connection for the PA14 to PA01 loci data
			if(all(sapply(tmpPair,function(x){ is.element(x, pa14to01$loci_PA14) }))){
				tmpLoci <- sapply(tmpPair,function(x){ pa14to01$loci_PA01[which(pa14to01$loci_PA14 == x)] },USE.NAMES=TRUE)
			} else {
				# We now look if there are PA###### Locus_Tag's identifiable for all elements of thisPair
				tmpLoci <- unlist(sapply(tmpPair,function(x){ 
					# We want to handle if nothing is found by returning an NA value
					tmpIndex <- which(pa01_genicDB$geneID == x)
					if(length(tmpIndex) == 0){
						return( NA )
					} else {
						return( pa01_genicDB$Locus_Tag[which(pa01_genicDB$geneID == x)] )
					} })) 
				# Now we handle circumstance where nothing was found bu geneID reference by using a blast search of the sequence
				if(is.element(NA,tmpLoci)){
					# First we check if there is a pa01_genicDB database else we make one
					if (!file.exists(paste(blastDB,"pa01_genicDB.nsq",sep=""))){
						print(paste("Missing the blastdb for pa01_genicDB trying to create it",sep=" "))
						# We need to make certain there is a Seq_nt sequence otherwise we won't include that to our blastDB
						tmpIndex_seqNT <- which(nchar(pa01_genicDB$Seq_nt) > 0)
						write.fasta(as.list(pa01_genicDB$Seq_nt[tmpIndex_seqNT]),names=pa01_genicDB$Locus_Tag[tmpIndex_seqNT],file=paste(blastDB,"pa01_genicDB.fas",sep=""))
						system(paste(blastLoc,"makeblastdb -in ",blastDB,"pa01_genicDB.fas -dbtype nucl -input_type fasta -out ",blastDB,"pa01_genicDB",sep=""))
						unlink( paste(blastDB,"pa01_genicDB.fas",sep="") )
					}
					# Now we will blast the Seq_nt for this gene from pa14 to the pa01 db
					for(thisGene in tmpPair[which(is.na(tmpLoci))]){
						# We hand the ENTREZ values as some instances of my siteRef's may not have the same caps values, this shouldn't happen for anything other than Entrez: values.
						tmpSeq <- pa14_genicDB$Seq_nt[which(toupper(pa14_genicDB$geneID)== toupper(thisGene))]
						tmpFile <- paste(blastDB,"tmp_",gsub("[[:punct:]]","_",thisGene),".fas",sep="")
						# If we found a sequence then we will perform a blast search for this
						write.fasta(as.list(tmpSeq),names=thisGene,file=tmpFile)
						system(paste(blastLoc,"blastn -query ", tmpFile," -db ",blastDB,"pa01_genicDB -out ", blastDB,"tmpQuery_blast.table -outfmt '6 sseqid sstart send qstart qend qlen evalue pident sseq'",sep=""))
						# Now we look if we found anything by searching is a _blast.table file exists and if it has content
						if(file.exists(paste(blastDB,"tmpQuery_blast.table",sep="")) && length(readLines(paste(blastDB,"tmpQuery_blast.table",sep=""))) > 0){
							# We read in the blast.table information and run a bit of quality control by verifying the percIdent value and eValue
							tmpFrame <- read.table(paste(blastDB,"tmpQuery_blast.table",sep=""),header=FALSE, stringsAsFactors = FALSE)
							# These colnames are based on my knowledge of how I am making the blast call
							colnames(tmpFrame) <- c("refID","refStart","refEnd","qStart","qEnd","qLen","eValue","PercIdent","refSeq")
							if(tmpFrame$eValue <= 1e-40 || tmpFrame$PercIdent >= 90){
								tmpLoci[thisGene] <- tmpFrame$refID[union(which(tmpFrame$eValue == min(tmpFrame$eValue)),which(tmpFrame$PercIdent == max(tmpFrame$PercIdent)))[1]]
							}
							# We no longer need this tmpQuery_blast.table information
							unlink( paste(blastDB,"tmpQuery_blast.table",sep="") )
						}
						# we also no longer need the tmp_thisGene.fas file
						unlink( tmpFile )	
					}
				} # This closes out the calls to blast and other cross-referencing resources
			} # This closes out using the cross-reference short speed pa14to01 object
				
			# Now we review if our quick reference object for PA01 and PA14 loci has these elements
			# If we've been able to identify PA01 loci for this_genePair then we will add the quick reference information
			for(thisLoci in tmpLoci){
				if(!is.na(thisLoci) && !is.element(thisLoci,pa14to01$loci_PA01)){
					pa14to01 <- rbind(pa14to01, c(names(tmpLoci)[which(tmpLoci == thisLoci)],thisLoci))
				}
			}	
			 
			#  We check that there are no more NA values in the tmpLoci, otherwise we've failed to find meaningfull paired information
			# and we will simply move on without recording anything.
			if(!is.element(NA,tmpLoci)){
				# This lists the rows of our interactions object which contain instance of each gene, the recursive call finds the intersect.
				tmpList <- lapply(tmpLoci, function(x){ union(which(grepl(x, geneInteractions_detailed$protein1)),
																which(grepl(x, geneInteractions_detailed$protein2))) })
				tmpIntersect <- recursiveIntersect(tmpList)																					
				# If we have found any tmp_intersectAction than we will return thisPair, else we do not
				if(length(tmpIntersect) > 0){	
					knownInteractions <- rbind(knownInteractions,cbind(geneInteractions_detailed[tmpIntersect, ],"genePair_ID"=paste(tmpPair,collapse=pairString)) ) 
					interactingPairs <- c(interactingPairs, which(unique_sigPairs_list == thisPair))
				} else {
					# This means that there is no intersection between the two loci and thus they do not interact
					known_notInteractions <- c(known_notInteractions, this_genePair)
				}
			}
		} # This closes out our object being in our previously verified knownInteractions or known_notInteraction objects	
	} 
	# We now save our knownInteractions object after removing any NA instances
	knownInteractions <- knownInteractions[which(!is.na(knownInteractions[,1])),]
	pa14to01 <- pa14to01[which(!is.na(pa14to01[,1])),]
	known_notInteractions <- known_notInteractions[which(!is.na(known_notInteractions))]
	save(knownInteractions, known_notInteractions, pa14to01,file=paste(baseDir,"PARuns/ppInteractions/knownInteractions_minScore_", min_intScore,".RData",sep=""))
	# Now since we have the same indexes of our unique_sigPairs_list and unique_sigPairs, we extract our subsetted sites with this
	interacting_uniquePairs <- unique_sigPairs[interactingPairs,]
	# We save this large object created as it is relatively time consuming.... and at this point deemed quite meaningful.
	save(interacting_uniquePairs,file=paste(postDir,"interacting_uniquePairs_minInt",min_intScore,".RData",sep=""))
}

############## THIS USES COG TERMS TO LOOK FOR PAIRS AMONG SIMILAR ORTHOLOGOUS GROUPS ##############
# We can look for clusetering by Orthologous Group using the pa14_genicDB$OrthologGrp level information
# NOTE: This information may be highly similar to what is found through the geneInteractions information
if(!file.exists(paste0(postDir,"reviewItems_q",ourAccept,".RData"))){
	sameOG_uniquePairs <- unlist(lapply(1:length(unique_sigPairs_list),function(x){
					thisPair <- sapply(strsplit(unique_sigPairs_list[x], pairString)[[1]],function(y){ substr(y,1,regexpr("_",y)-1) })
					# We find the POG information for each of these pairs
					pogPair <- list("first" = pa14_genicDB$OrthologGrp[which(pa14_genicDB$geneID == thisPair[1])],
									"second" = pa14_genicDB$OrthologGrp[which(pa14_genicDB$geneID == thisPair[2])])
					#if(length(pogPair[["first"]]) > 1){
					#	print(paste0("Review this ",x))
					#} else if(length(pogPair[["second"]]) > 1){ 
					#	print(paste0("Review this ",x))
					#}
					# Now because within our dataset there are some sets where there is no POG information we verify this first
					if(all(unlist(lapply(unlist(pogPair),nchar)) > 0)){
						# We look to see if the Orthologous Group tags are identical based on the geneID information
						if( any(is.element(pogPair[["first"]], pogPair[["second"]])) ){
							return( x )
						} else {
							return( NULL )
						}
					} else {
						return( NULL ) 
					} }))
	sameOG <- unique_sigPairs[sameOG_uniquePairs,]
	##### If you care to see what is novel within this set
	#unique_sigPairs[setdiff(sameOG_uniquePairs,interactingPairs),]
				
		# Then we can look if the COG terms using pa14_genicDB$COG and look for pairs that are within the same COG or category 
		# (through the pa01_cog.csv reference - which derives from Gallperin, 2015)
		# Being part of the same COG group would be very weak support, what I'd want is for the two proteins to be in a correlated pathway (hence gene interactions...)
		################## I'VE THOUGHT ABOUT THIS AND AM NOT CONVINCED IT WILL ACTUALLY YIELD ANY USEFULL INFROMATION ###############
		
		# Also we can look at sub-ceullar localisation with pa14_genicDB$Subcellular_Localization
		################# THE INFORMATION WITHIN OUR DATABASE IS SIMPLY TOO SIMPLE TO HAVE THIS BE A MEANINGFUL COMPARISON
	
		# We can also look at biological processes as defined in the pa14_genicDB$Gene_Ontology information, 
		# We can further divide this into the [biological process] and [molecular function] information
		###### CAN'T THERE SIMPLY AREN'T ENOUGH PROTEINS WITH INFORMATION HELD HERE ##########
	
	
	############## NOW WE CAN START PERFORMING SOME ANALYSIS OF OUR SIGPAIRS WHICH HAVE SHOWN BIOLOGICAL RATIONAL ##############
	combinedPairs <- unique(rbind(sameOG, interacting_uniquePairs))
	# Let's review again the distribution of q-Value information within this set
	combined_sigLevels <- list()
	tmpAccept <- ourAccept
	while(length(which(combinedPairs$pValue_FDR <= tmpAccept)) > 0 && !all(combinedPairs$pValue_FDR[which(combinedPairs$pValue_FDR <= tmpAccept)] == 0) ){
		combined_sigLevels[[as.character(tmpAccept)]] <- length(which(combinedPairs$pValue_FDR <= tmpAccept))
		tmpAccept <- tmpAccept/10
	}
	
	# This is just a call to save some of the review objects that have been created to date and does not need to persist.
	save(combinedPairs, sameOG, combined_sigLevels, sigLevels, unique_sigPairs_list, file=paste0(postDir,"reviewItems_q",ourAccept,".RData"))
} else {
	load(paste0(postDir,"reviewItems_q",ourAccept,".RData"))
}




################################################################################################################
################################################################################################################
############## HERE WE CONSIDER LOOKING AT THE P-VALUE OBTAINED AGAINST A SITE PATTERN AND #####################
################################################################################################################
################################################################################################################
# Now we look to see if the have an object saved which informs us of our interacting_unique_sigPairs
### NOTE: This overwrites our unique_sigPairs object and is therefore accepting this object's contents as a 
###			hiearchial filter for pairs being meaningfull
if(is.element("interacting_uniquePairs",ls()) && use_supportSTRING){
	assign("unique_sigPairs", interacting_uniquePairs, pos=".GlobalEnv")
}

# I now run a re-setting of the unique_sigPairs in case there are duplicates the crept in from this step.... (found in past)
# Now we define the unique significant pairs, first we get a list, then for each list entry we find which row has the smallest pValue_FDR
duplicateRows <- find_duplicateRows(unique_sigPairs,"siteRef_1","siteRef_2")
removeRows <- if (length(duplicateRows) == 0){ 
						NULL 
				} else { 
					unlist(sapply(1:length(duplicateRows),function(x){
						# We look if this list entry is NULL, if not we need to review which to keep
						if (length(duplicateRows[[x]]) > 0){
							tmpSet <- c(as.numeric(names(duplicateRows)[x]),duplicateRows[[x]])
							# now we return the non-minimum pValue_FDR rows as reported by sigPairs row x or those duplicates.
							return( tmpSet[-which(unique_sigPairs$pValue_FDR[tmpSet] == min(unique_sigPairs$pValue_FDR[tmpSet]))[1]] )
						} else {
							return( NULL )
						} })) 
				}
unique_sigPairs <- if (length(removeRows) > 0){ unique_sigPairs[-removeRows,] } else { unique_sigPairs }


# This will check if we need to recollect some farmed parts or if we will simply be taking our saved object
if(!file.exists(paste(postDir,"pairedTests_ourAccept",ourAccept,".RData",sep=""))){
	allLists <- list.files(path=hpcvl_farmOut, pattern=paste("pairedTests_ourAccept",ourAccept,sep=""),full.names=TRUE)
	pairedTests <- NULL
	# If nothing was returned there is nothing to recollect...
	if(length(allLists) > 1){
		listNames <- sapply(1:length(allLists),function(x){ paste(rep("tmp",x),collapse="_") })
		for(x in 1:length(allLists)){ 
			load(allLists[x])
			assign(listNames[x],pairedTests,pos=".GlobalEnv") 
		}
		# We assume all parts will have the same shape as the first instance
		combinedTests <- eval(as.name(listNames[1]))
		# We now look at the names and sub-names of each loaded element
		#allTests <- lapply(listNames[-1],function(x){ return(matrix(c(names(eval(as.name(x))),unlist(lapply(eval(as.name(x)),length))),ncol=2)) })
		allTests <- lapply(listNames,function(x){ return(matrix(c(names(eval(as.name(x))),unlist(lapply(eval(as.name(x)),length))),ncol=2)) })
		
		# We will now go through all unique names of our leadingSites (first column in the returned matrix) and add unique elements to our combinedTests
		for(thisSite in unique(unlist(lapply(allTests,function(x){ return(x[,1]) })))){
			# We go through each allTests element and seek if thisSite is within it
			#for(thisLoaded in listNames[-1]){
			for(thisLoaded in listNames){
			
				tmpData <- eval(as.name(thisLoaded))
				if(is.element(thisSite,names(tmpData))){
					# We check if thisSite is already part of combinedTests
					if(is.element(thisSite,names(combinedTests))){
						# We now combine new elements in tmpData to our combinedTests, but to preserve the names of elements we perform it this slower means
						# rather than simply calling unique()
						tmpNew <- lapply(names(tmpData[[thisSite]]),function(x){ 
									if(!is.element(x,names(combinedTests[[thisSite]]))){ 
										tmpReturn <- tmpData[[thisSite]][[x]]
										#names(tmpReturn) <- x
										return( tmpReturn )
									} else {
										return( NULL ) 
									} })
						names(tmpNew) <- names(tmpData[[thisSite]])			
						# Now we remove any NULL elements of tmpNew
						for(thisElement in names(tmpNew)){
							if(is.null(tmpNew[[thisElement]])){
								tmpNew[[thisElement]] <- NULL
							}
						}
						# If there is anything new to add do it, otherwise perform no action
						if(length(tmpNew) > 0){
							combinedTests[[thisSite]] <- c(combinedTests[[thisSite]], tmpNew)
						}
					} else {
						# This means we simple add a new element to combinedTests
						combinedTests[[thisSite]] <- tmpData[[thisSite]]
					}
				}
			}
		}
		pairedTests <- combinedTests
		rm(list= listNames)
		# This is the old method which can't be trusted unless each leadingSite is uniquely searched and basic list is blank.		
		# We now combine this basic pairedTests list with all other parts looking for unique elements in each to append
		#combinedTests <- lapply(1:length(combinedTests),function(x){ 
		#					# We look if this list is empty, if so we look for data in the other tmpNames elements
		#					if(length(combinedTests[[x]]) == 0){ 
		#						useList <- which(sapply(listNames[-1],function(y){ length(eval(as.name(y))[[x]]) > 0 }))
		#						if(length(useList) > 0){
		#							return( eval(as.name(listNames[useList + 1]))[[x]] )
		#						} else {
		#							return( combinedTests[[x]] )
		#						}
		#					} else { 
		#						return( combinedTests[[x]] ) 
		#					} })
		#pairedTests <- combinedTests
		#names(pairedTests) <- names(eval(as.name(listNames[1])))
	
					### OLD ####### I NEED TO MERGE MY CURRENT WORKED OUT pairedTests WITH THOSE REQUIRED BY THIS CURRENT FILTER STATUS ###
					## this should come from my old lists
					#keepLists <- pairedTests
					## this is my new list
					#pairedTests <- sapply(unique(unique_sigPairs$siteRef_2),simplify="list",FUN = function(x){ return(vector(mode="list",length=0)) })
					## Now we go through each of the new list's leading site's and look for those pairs that ought to exist, if they are in the keepLists
					## set then we import them
					#for(this_leadPair in names(pairedTests)){
					#	allPaired <- unique(unique_sigPairs$siteRef_1[which(unique_sigPairs$siteRef_2 == this_leadPair)])
					#	for(thisPaired in allPaired){
					#		if(is.element(thisPaired,names(keepLists[[this_leadPair]]))){
					#			# This means we've found a list entry that alraedy exists, so we will port it over
					#			pairedTests[[this_leadPair]][[thisPaired]] <- keepLists[[this_leadPair]][[thisPaired]]
					#		}
					#	}
					#}
 		save(pairedTests,file=paste(postDir,"pairedTests_ourAccept",ourAccept,".RData",sep=""))
	} else {
		load(allLists[x])
	}
} else {
	load(paste(postDir,"pairedTests_ourAccept",ourAccept,".RData",sep=""))
} ### HERE ####

# We start by collecting the siteRef information for all unique site patterns, RECALL that some of these are collapsed
# sitePair information, thus when we recollect the traitMat sitePattern information we'll only use the first value but assign
# the result to all related sitePairs.
all_uniqueSites <- unique(unlist(c(unique_sigPairs$siteRef_1,unique_sigPairs$siteRef_2)))
all_uniqueSites <- all_uniqueSites[order(all_uniqueSites)]



##############################

# Now we will go through each unique site pattern and then find the distribution of p-values obtained from runing it against BayesTraits with min_testReps
# randomisations of the paired sitePattern.  Provided the tree is sufficiently long, and the site pattern is not largely monomorphic, then there should be 
# many permutations of the site pattern which are unique.
if(file.exists(paste(postDir,"pairedTests_ourAccept",ourAccept,".RData",sep=""))){
	load(paste(postDir,"pairedTests_ourAccept",ourAccept,".RData",sep=""))
} else {
	pairedTests <- sapply(unique(unique_sigPairs$siteRef_2),simplify="list",FUN = function(x){ return(vector(mode="list",length=0)) }) 
}

##################################################################################################################
##################### WE START BY LOOKING FOR SIGNIFICANT LD AND REMOVING INTRAGENEIC PAIRS ######################
##################################################################################################################

if(!file.exists(paste(postDir,"filteredLD_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))){
	tmpData <- vector(mode="list",length=nrow(unique_sigPairs))
	for(x in 1:nrow(unique_sigPairs)){
		thisPair <- unlist(unique_sigPairs[x, siteRef_colSet])
		# In case we have sites which had non-unqijue site patterns we simply consider the first instance in this assessment
		thisPair <- sapply(thisPair, function(x){ strsplit(x,"::")[[1]][1] })
		# We extract the sitePattern for these two positions
		tmp_sitePatterns <- cbind(traitMat[,which(posNames == thisPair[1])[1]],traitMat[,which(posNames == thisPair[2])[1]])
		# We remove any gap characters, which are stored as 2.
		tmp_sitePatterns <- tmp_sitePatterns[intersect(which(tmp_sitePatterns[,1] != 2),which(tmp_sitePatterns[,2] != 2)),]
		# We now need to check that there is diveristy among our site patterns else we simply move on...
		if(any(apply(tmp_sitePatterns,MARGIN=2,function(x){ length(unique(x)) <= 1 }))){
			tmpData[[x]] <- list("D"=NA,"D'"= NA,"r"=NA,"R^2"=NA,"n"=NA,"X^2"=NA,"P-value"=NA)
		} else {
			# We now find our which of the site pairs is sensitive and resistant
			# This will create genotype objects which are the "changes" based on the root state 
			tmpData[[x]] <- LD(genotype(paste(tmp_sitePatterns[,1],tmp_sitePatterns[,1],sep=""),sep=""),genotype(paste(tmp_sitePatterns[,2],tmp_sitePatterns[,2],sep=""),sep=""))
		}
	}
	
	# Now I get the distance data for my pairs, since this is a circular genome we must know the maximum length of the reference genome 
	# and look at the smallest value of the linear distance between them, and the sum of the smallest position and the diff of the largest and genome length
	# This will require us to quickly look at the length of the reference genome
	tmpLen_refGen <- ncol(read.dna(paste(baseDir,"/PARuns/Sequences/PA_wholeGenomes/allGenomes/PA14.fas",sep=""),format="fasta"))
	tmpDist <- as.numeric(apply(unique_sigPairs[, siteRef_colSet], MARGIN = 1, function(thisPair){
						tmpGenes <- sapply(sapply(thisPair,function(y){ strsplit(y,"::")[[1]][1]}), function(x){ strsplit(x,"_")[[1]] })
						tmpPositions <- apply(tmpGenes, MARGIN = 2, function(x){ pa14_genicDB$Start[which(pa14_genicDB$geneID == x[1])[1]] - 1 + as.numeric(x[2]) })
						# Now that we have the positions, in PA14 reference genome of each of these mutations, we find the smallest of either
						# their absolute difference, OR the sum of the smallest and the difference between the largest and the length of the genome
						return( min(abs(diff(tmpPositions)),sum(min(tmpPositions), tmpLen_refGen - max(tmpPositions) + 1)) ) }))
	tmpReview <- cbind(unique_sigPairs[, c(siteRef_colSet,"pValue_FDR")], tmpDist, sapply(tmpData,function(x){ return(x[["P-value"]])}))
	colnames(tmpReview)[4:5] <- c("Dist","LD_pValue")
	
	is_intraGenic <- apply(unique_sigPairs[, siteRef_colSet], MARGIN = 1, function(x){ strsplit(x[1],"_")[[1]][1] == strsplit(x[2],"_")[[1]][1] })
	
	# We will test to see if the signal of correlated evolution is at all a signal of LD or distance
	# Which is likely a bad model since we assume that Ld and distance are correalted, however, there is no statistical linear correlation
	ldModel <- lm(LD_pValue ~ Dist, data = tmpReview)
	corrModel <- lm(pValue_FDR ~ Dist * LD_pValue, data = tmpReview)
	# There is evidence that some points have high leverage, the R^2 or corrModel is 0.002.. so ignored.
	
	pdf(paste(postDir,"plots_distLD.pdf",sep=""), height = 12, width = 12)
	par(mfrow=c(2,2))
	plot(tmpReview[,"LD_pValue"],unique_sigPairs$pValue_FDR, xlab = "LD P-Value", ylab = "Corr. Evol. qValue", log="y")
	lines(density(tmpReview[,"LD_pValue"]),col="red")
	plot(tmpDist,unique_sigPairs$pValue_FDR, xlab = "Physical Distance (nt positions)", ylab = "Corr. Evol. qValue",log = "y")
	plot(tmpDist[which(is_intraGenic)],tmpReview[which(is_intraGenic),"LD_pValue"],ylab = "LD P-Value", xlab = "Physical Distance (nt positions)")
	#plot(tmpDist[which(is_intraGenic)],unique_sigPairs$pValue_FDR[which(is_intraGenic)], xlab = "Physical Distance (nt positions)", ylab="Corr. Evol. qValue",log="y")
	# We now visualise our data plotted in 3D space, we colour based on gene types being in the set.
	tmp_qVals <- -log10(unique_sigPairs$pValue_FDR); tmpMax <- max(tmp_qVals[which(tmp_qVals != "Inf")])*1.2
	for(thisVal in which(tmp_qVals == Inf)){ tmp_qVals[thisVal] <- tmpMax }
	scatter3D(tmpDist, tmpReview[,"LD_pValue"], tmp_qVals, xlab="Physical Distance (nt positions)", 
				ylab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 300, phi = 30, pch = 19, cex.lab = 1)
	dev.off()
	
	# This is a scatterplot of the points based on the LD values, signal of corr evol, and the maximum likelihood probability of the root state
	# We use the siteRef_# information to reference the unique sig pairs data frame and find the maximum likelihood root state recorded.
	# We do this for both the dependent and independent root states
	tmpReview <- cbind(tmpReview, t(apply(tmpReview[,siteRef_colSet], MARGIN = 1, function(x){ 
						tmpRow <- intersect(which(unique_sigPairs$siteRef_1 == x[1]),which(unique_sigPairs$siteRef_2 == x[2]))
						return( c("mlIndep" = max(unique_sigPairs[tmpRow,intersect(which(grepl("rootP_",colnames(unique_sigPairs))),which(grepl("_indep",colnames(unique_sigPairs))))]), 
									"mlDep" = max(unique_sigPairs[tmpRow,intersect(which(grepl("rootP_",colnames(unique_sigPairs))),which(grepl("_dep",colnames(unique_sigPairs))))])) ) 
						})))
	pdf(paste(postDir,"plots_rootPvsLD.pdf",sep=""), height = 6, width = 18)
	par(mfrow=c(1,3))
	tmp_qVals <- -log10(unique_sigPairs$pValue_FDR); tmpMax <- max(tmp_qVals[which(tmp_qVals != "Inf")])*1.2
	for(thisVal in which(tmp_qVals == Inf)){ tmp_qVals[thisVal] <- tmpMax }
	scatter3D(tmpReview[,"LD_pValue"], tmpReview[,"mlIndep"], tmp_qVals, ylab="Root State ML (indep model)", 
				xlab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)
	scatter3D(tmpReview[,"LD_pValue"], tmpReview[,"mlDep"], tmp_qVals, ylab="Root State ML (dep model)", 
				xlab = "LD P-Value", zlab = "Corr. Evol. qValue (-log10)", ticktype = "detailed",
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)
	# Here we plot the colours as the qValue
	scatter3D(tmpReview[,"mlIndep"], tmpReview[,"mlDep"], tmpReview[,"LD_pValue"], ylab="Root State ML (dep model)", 
				xlab = "Root State ML (indep model)", zlab = "LD P-Value", ticktype = "detailed", colvar = tmp_qVals,
				theta = 335, phi = 30, pch = 19, cex.lab = 1, alpha = 0.8)			
	dev.off()
	
	
	# From having a look at this we will filter by removing those sites which show significant LD AND are intra-genic
	unique_sigPairs <- unique_sigPairs[union(which(!is_intraGenic),intersect(which(is_intraGenic),which(tmpReview[,"LD_pValue"] > 0.05))),]
	
	save(unique_sigPairs, ldModel, corrModel, tmpReview, file=paste(postDir,"filteredLD_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))
} else {
	load(paste(postDir,"filteredLD_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))
}
####### THIS IS CODE FROM A PREVIOUS METHOD ##########
## This vector will let us know is a strain is sensitive, resistant, or unknown
#load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")	
#kos_phenoVec <- unlist(sapply(rownames(traitMat),function(x){ 
#if (is.element(x, c("PA14","PA01"))){ 
#	return( "ND" )
#} else {
#	as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
#} } ))
#rm(bacFrame)
# This vector and list will record if a strain is resistant, or other.
#sensitivityVec <- sapply(kos_phenoVec,function(x){ is.element("(R)",x) })
#apply(unique_sigPairs[, siteRef_colSet], MARGIN = 1, function(thisPair){
#		# thisPair is sent forward as a data.frame which is just not practical we want it as a vector...
#		thisPair <- unlist(thisPair)
#		# We extract the sitePattern for these two positions
#		tmp_sitePatterns <- cbind(traitMat[,which(posNames == thisPair[1])],traitMat[,which(posNames == thisPair[2])])
#		
#		# This will create genotype objects which are the "changes" based on the root state 
#		LD(genotype(tmp_sensitivityList$Sensitive,sep=""),genotype(sample(tmp_sensitivityList$Resistant, length(tmp_sensitivityList$Sensitive), replace=FALSE),sep=""))
#		
#		# We now find our which of the site pairs is sensitive and resistant
#		tmpPairs <- apply(tmp_sitePatterns, MARGIN = 1, function(x){ paste(x,collapse = "") })
#		tmp_sensitivityList <- list("Sensitive"= tmpPairs[which(sensitivityVec)], "Resistant"= tmpPairs[which(! sensitivityVec)])
#		# Now we can calculate D, the linkage disequilibrium, and from that the R^2 value (coefficient of correlation)
#		tmpD <- lapply(tmp_sensitivityList, function(x){ 
#						tmpVals <- c("00"=length(which(x == "00")), 
#										"01"=length(which(x == "01")), 
#										"10"=length(which(x == "10")), 
#										"11"=length(which(x == "11")))/length(x)
#						return(tmpVals["11"]*tmpVals["00"] - tmpVals["10"]*tmpVals["01"]) })
#		tmp_rSqr <- lapply(names(tmp_sensitivityList), function(x){ 
#						tmpRows <- which(is.element(rownames(tmp_sitePatterns),names(tmp_sensitivityList[[x]]))) 
#						tmpVals <- c("p0"=length(which(tmp_sitePatterns[tmpRows,1] == 0)), 
#									"p1"=length(which(tmp_sitePatterns[tmpRows,1] == 1)), 
#									"q0"=length(which(tmp_sitePatterns[tmpRows,2] == 0)), 
#									"q1"=length(which(tmp_sitePatterns[tmpRows,2] == 1)))/length(tmpRows) 
#						return( tmpD[[x]]^2 /(tmpVals[1]*tmpVals[2]*tmpVals[3]*tmpVals[4]) ) })
#		# Now we know the correlation of coefficient for both the sensitive and resistant strains
#})

# We see that the sensitive type for this pair has a strong signal of linkage equilibrium, and for resistant this is moderate.

	
##################################################################################################################
###################### NOW WE ASK IF WE FARM TO HPCVL OR IF WE SIMPLY RUN FROM THIS SCRIPT #######################
##################################################################################################################
if(farmHPCVL){
	# First we find which parts of the pairedTests object need to be worked upon
	sites_needWork <- names(which(sapply(names(pairedTests),function(x){ 
		# The first asks if the length of our site lists match, the seonc asks if each instance has min_testReps
		if(length(pairedTests[[x]]) != length(unique(unique_sigPairs$siteRef_1[which(unique_sigPairs$siteRef_2 == x)]))){
			return( TRUE )
		} else if(!all(sapply(names(pairedTests[[x]]),function(y){ nrow(as.matrix(pairedTests[[x]][[y]])) == 1000 }))) { 
			return( TRUE )
		} else {
			return( FALSE )
		} })))
	# Now we divide the sites that need work into parts, we start by identifying the lead elements
	all_workingLeads <- unlist(lapply(sites_needWork, function(x){ names(pairedTests)[which(names(pairedTests) == x)] }))
	# Now for each lead we identify which paired sites need work
	all_workingPairs <- list()
	# We go through all leads and find which potential pairs have not yet been run successfully.
	for(thisLead in all_workingLeads){
		# This exsitst to allow us to re-start these runs at mid-points, we verify if this list obejct has the proper length
		# we look to see if there have been sub lists created and if they have been filled	
		if(length(pairedTests[[thisLead]]) != 0){
			# Otherwise we want to see if there has been sufficient work done, and mark if we need to run a pair
			for(thisPaired in unique(unique_sigPairs$siteRef_1[which(unique_sigPairs$siteRef_2 == thisLead)])){
				# Then we check if the list has named objects and if this name exists
				if(is.element(thisPaired, names(pairedTests[[thisLead]]))){
					# If this paired site exists, we check if it has been run to completion.
					if(nrow(pairedTests[[thisLead]][[thisPaired]]) >= min_testReps){
						next
					} else {
						# This means we need to perform some work
						all_workingPairs[[thisLead]] <- c(all_workingPairs[[thisLead]], thisPaired)
					}
				} else {
					# If this paired site does not have a name in the pairedTests list then we need it run
					all_workingPairs[[thisLead]] <- c(all_workingPairs[[thisLead]], thisPaired)
				}
			}
		} else {
			# If there is no work that has yet been done for thisLead, then we need to do all possible partners
			all_workingPairs[[thisLead]] <- unique(unique_sigPairs$siteRef_1[which(unique_sigPairs$siteRef_2 == thisLead)])
		}
	}
	# Now the number of runs requires is calculated, as well the appropriate division of labour (1st and 2nd index elements created)
	tmpDims <- c(length(unlist(all_workingPairs)),ceiling(length(unlist(all_workingPairs))/num_farmParts))
	# We know roughly hw big each part should be, now we go about finding the division of lead sites which will accomodate this division
	workSizes <- sapply(all_workingPairs, length)
	
	# Now we assign elements to the farmable parts
	aCounter <- 1
	theseDivisions <- list()
	for(thisDivision in 1: num_farmParts){
		counterStart <- aCounter 
		thisPart <- 0
		while(thisPart <= tmpDims[2] && aCounter <= length(all_workingPairs)){
			thisPart <- sum(thisPart, workSizes[aCounter])	
			aCounter <- aCounter + 1
		}
		theseDivisions[[thisDivision]] <- names(all_workingPairs)[counterStart:aCounter][which(!is.na(names(all_workingPairs)[counterStart:aCounter]))]
	}
	
	
	# Then we will modify our template script and then send requests to the server of HPCVL to work
	for(thisPart in 1: num_farmParts){
		tmpLines <- readLines(paste(postDir, hpcvlScript,sep=""))
		# This indicates which "part" - a unique save an indexing means - this set will be
		tmpLines[which(grepl("farmPart <- ",tmpLines))[1]] <- paste("farmPart <- ",thisPart,sep="")
		# This will identify which elements of the pairedTests list should be worked upon
		tmpLines[which(grepl("farmRange <- ",tmpLines))[1]] <- paste("farmRange <- ", paste(theseDivisions[thisPart],collapse=" "),sep="")
		# This sets the number of cores this job ought to request
		tmpLines[which(grepl('registerDoMC',tmpLines))[1]] <- paste("registerDoMC(cores=", num_hpcvlCores,")",sep="")
		# These are other variables we wish to reset
		tmpLines[which(grepl("ourAccept <- ",tmpLines))[1]] <- paste("ourAccept <- ",ourAccept,sep="")
		tmpLines[which(grepl("percDiffs <- ",tmpLines))[1]] <- paste("percDiffs <- ",percDiffs,sep="")
		tmpLines[which(grepl("min_testReps <- ",tmpLines))[1]] <- paste("min_testReps <- ",min_testReps,sep="")
		tmpLines[which(grepl("pairString <- ",tmpLines))[1]] <- paste('pairString <- "',pairString,'"',sep="")
		tmpLines[which(grepl("farmHPCVL <- ",tmpLines))[1]] <- paste("farmHPCVL <- ",TRUE,sep="")
		writeLines(tmpLines,con=paste(postDir, sub("_mod2.r",paste("_part", thisPart,".r",sep=""),hpcvlScript),sep=""))
		# Now we submit this script as a job to the HPCVL server
		#system(paste("qsub -q abaqus.q -pe dist.pe ", num_hpcvlCores," -N reNew_",thisPart," -cwd -b y -j y /home/hpc3058/bin/R-3.2.3/bin/R CMD BATCH ",
		#				postDir, sub("_mod.r",paste("_mod",thisPart,".r",sep=""),hpcvlScript),sep=""))
	}
	# We now exit this script as it has sent jobs to be worked that must complete before any further steps are taken
	print("Jobs ought to be submitted to the server queue, now exiting")
	q(save="no")
} else {
	# For each unqiue site we will test it against all sites for which it was found to be paired and get a distribution of the pValues for randomised
	# pairedSite site patterns.
	for(thisSite in unique(unique_sigPairs$siteRef_2)){
		for(thisPaired in unique(unique_sigPairs$siteRef_1[which(unique_sigPairs$siteRef_2 == thisSite)])){
			# This exsitst to allow us to re-start these runs at mid-points, we verify if this list obejct has the proper length
			# we look to see if there have been sub lists created and if they have been filled
			if(length(pairedTests[[thisSite]]) != 0){
				# Then we check if the list has named objects and if this name exists
				if(is.element(thisPaired, names(pairedTests[[thisSite]]))){
					if(nrow(pairedTests[[thisSite]][[thisPaired]]) >= min_testReps){
						next
					}
				}
			}
			# we collect the site patterns to use and recall when using the siteRef information to strsplit this 
			# in the event it is a collapsed string.
			thisPattern <- traitMat[,which(posNames == strsplit(thisSite, pairString)[[1]][1])]
			pairedPattern <- traitMat[,which(posNames == strsplit(thisPaired, pairString)[[1]][1])]
			# We want to create min_testReps number of unique samplings of the order
			testPatterns <- unique(sapply(1: min_testReps,function(x){ sample.int(length(tTree$tip.label)) }))
			if(!is.matrix(testPatterns)){
				print("There was a problem with trying to create a unique set of testPatterns, did not return a matrix")
				q(save="no")
			}
			# If the sampling did not produce a proper number of entirely unqiue samplings then we will reshuffle
			aCounter <- 1
			while(aCounter < 1000 && ncol(testPatterns) < min_testReps){
				testPatterns <- unique(sapply(1: min_testReps,function(x){ sample.int(length(tTree$tip.label)) }))
				aCounter <- aCounter + 1
			}
			# We verify if our method was able to produce a unique (for ordering not actual site pattern) set of sitePattern combinations.
			if(aCounter >= 1000 && ncol(testPatterns) < min_testReps){
				print("There was a problem with trying to create a unique set of testPatterns after 1000 replications, review please")
				q(save="no")
			}
			# Now we take our test patterns and generate a BayesTraits pValue for each of these
			test_pValues <- foreach(randomPattern = 1: min_testReps,.combine="rbind")%dopar%{
				# we now put in our normal BayesTraits loop and use the columns as our site patterns, the test site being a randomisation
				col1 <- thisPattern
				# Here we have traitMat sites as the second site against which we compare
				col2 <- pairedPattern[testPatterns[,randomPattern]]
				# This will handle any possible values of 2 which have persisted through the builtMat and siteCoulmn writting of previous steps.
				# NOTE: mostly a failsafe
				if (length(which(col1 == "2")) > 0){ col1[which(col1 == 2)] <- "-" }
				if (length(which(col2 == "2")) > 0){ col2[which(col2 == 2)] <- "-" } 
				col1[is.na(col1)] <- "-"; col2[is.na(col2)] <- "-"; col1[which(col1 == "NA")] <- "-"; col2[which(col2 == "NA")] <- "-"
				# These are objects to hold our siteColumns which we will remove dashes from and evalutate if the likelihood surface would be flat due to identity
				col1bis <- col1 
				col2bis <- col2
				# remove cols w/ "-" if there is an "-" in either site both sites are NAed
				col1bis[which(col2bis=="-")] <- NA
				col2bis[which(col2bis=="-")] <- NA
				col2bis[which(col1bis=="-")] <- NA
				col1bis[which(col1bis=="-")] <- NA
				# This call to identical will remove identical sites which are of this type because of initial values or after the above NA steps.
				if(!identical(col1bis,col2bis)){
					# We build a dataframe of what we want analysed and reorder the columns based on the corrected order from corrOrder
					trait_data <- data.frame(tTree$tip.label, col1[corrOrder], col2[corrOrder])
					colnames(trait_data) <- NULL
					write.table(trait_data, paste(tmpDir,"traits.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""),quote=F,row.names=F)
					
					# run BayesTraits under the dependent and indep models
					########### NOTE: The Input files InputBT_ML_indep.txt, InputBT_ML_dep.txt must be placed into each wd
					system(paste(bayesLoc, "BayesTraits ",workDir,treeFile," ",tmpDir,"traits.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt < ", workDir,"InputBT_ML_indep.txt > ",tmpDir,"traits_res_indep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					system(paste(bayesLoc, "BayesTraits ",workDir,treeFile," ",tmpDir,"traits.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt < ", workDir,"InputBT_ML_dep.txt > ",tmpDir,"traits_res_dep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					
					# extract likelihood values, by counting the # of lines in the traits_res_indep.i.j.txt file
					system(paste("wc -l ",tmpDir,"traits_res_indep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt > ",tmpDir,"wc", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					# This reads the wc file produce which should have simply the number of lines in the file
					wc <- read.table(paste(tmpDir,"wc", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))            
					# Since the last line in the "traits" file is the valuable output from ByesTraits, we skip wc[1] -1 lines, to read the last.
					res_indep <- read.table(paste(tmpDir,"traits_res_indep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""),header=FALSE,fill=TRUE,skip=as.numeric(wc[1]-wcAdj))
					# Repeat now for the dependent model
					system(paste("wc -l ",tmpDir,"traits_res_dep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt > ",tmpDir,"wc", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					wc <- read.table(paste(tmpDir,"wc", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))           
					res_dep <- read.table(paste(tmpDir,"traits_res_dep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""),header=FALSE,fill=TRUE,skip=as.numeric(wc[1]-wcAdj))
					# This adjustment is required only beacuse on the HPCVL server the traits_res_indep & _dep files contain an additional line
					# after results section leading to a table read in with 2 rows, yet we only need to 1 row's value
				
					res_dep <- res_dep[1,]
					colnames(res_dep) <- c("TreeNo","Lh","q12","q13","q21","q24","q31","q34","q42","q43","P(00)","P(01)","P(10)","P(11)")
					res_indep <- res_indep[1,]
					colnames(res_indep) <- c("TreeNo","Lh","alpha1","beta1","alpha2","beta2","P(00)","P(01)","P(10)","P(11)")
	
					# cleanup tmp files releasing them from memory
					unlink(paste(tmpDir,"traits_res_indep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					unlink(paste(tmpDir,"traits_res_dep.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					unlink(paste(tmpDir,"traits.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					unlink(paste(tmpDir,"traits.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt.log.txt",sep=""))
					unlink(paste(tmpDir,"wc", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt",sep=""))
					unlink(paste(tmpDir,"traits.", gsub("[[:punct:]]","",thisSite),".", randomPattern,".txt.log.txt.Schedule.txt",sep=""))
					
					# Since we expect our Bayes values to be negative this tests to assure us we are within expectations otherwise ignore
					# Note these likelihood values are the log-likelihoods of the models.
					if((res_indep[,2] < 0) & (res_dep[,2] < 0)){
					    # As this is a two sided test we multiply our test stat by two
					    test_stat <- 2 * abs(res_indep[,2] - res_dep[,2])
					    # this should be the last instruction of the foreach loop.  This caculatees the significance of our LRT
					    return( cbind(data.frame("pValue"=1-pchisq(test_stat,4),"testStat"= test_stat), res_indep[,3:6],res_dep[,3:10],stringsAsFactors=FALSE) )
					}
				}
			} # This closes out the inner loops of our BayesTraits analysis for a paired pattern
			# Now we will save the pValues obtained for thisSite and it's pairedSite
			pairedTests[[thisSite]][[thisPaired]] <- test_pValues
			# We will save this object as we grow it, this will allow us to extract and review mid points
			save(pairedTests, file=paste(postDir,"pairedTests_ourAccept",ourAccept,".RData",sep=""))
		}
	}
}


# These are methods for reviewing the output from this stage
# These are dimensions/parameters for our plotting
numPlots <- sum(unlist(lapply(pairedTests,function(x){ length(names(x)) })))
numCols <- min(ceiling(sqrt(numPlots)),8)
numRows <- min(ceiling(numPlots/numCols)+1,10)
numBreaks <- 100
numPoints <- length(pairedTests[[1]][[1]]$testStat)

# Now make a plot for the test stat of all values considered and tested.
pdf(paste(postDir,"pairedTests_ourAccept", ourAccept,"_testStat.pdf",sep=""),height=4 * numRows,width= 6 * numCols)
par(mfrow=c(min(10,numRows), min(8,numCols)))
# We plot some chisq distributions with the same number of breaks as our testStat data to compare
for(x in 1: numCols){
	hist(rchisq(numPoints,4),breaks= numBreaks,main=paste("NULL",x,sep=" "),cex.main = 1.5)
}
for(thisPlot1 in names(pairedTests)){
	for(thisPlot2 in names(pairedTests[[thisPlot1]])){
		hist(pairedTests[[thisPlot1]][[thisPlot2]]$testStat,breaks= numBreaks,main=paste(thisPlot1,thisPlot2,sep="  and  "),cex.main=1.5)
	}
}
dev.off()

## HOWEVER: This more complex reivew is still ideal, use the Li et al 2016 timeTree_ofLife info to test and workout a script.
## NOTE: From my work so far, the keep pairs list will be essentially zero and does not offer an exceptionally meaningfull filter method.
##		I may only want to bias parameter estimates which relate to transition from the root to observed states.  However, teasing out what the 
##		proportion of site pairs are in which state, other than a simple table of pasted vectors function, is the observed values of interest.  We will
##		have two sets, the (00, 11) and the (01, 10), so if we look for which pair is in maximal predominance will inform which parameters need to be 
##		propoerly estimated to generate the observed data with a proper test statistic.  (use the compiledRun_data.RData object to get the
##	 	root_<statePair>_<model_type> columns of each sitePair's root information, then the geneSites and sitePatterns.RData objects for observed.
## IMPORTANT: From a review of the gyrA_248 and parC_2006, there are columns of parameters which quite literally plateau as ~100% 0 or 1000.  Need to review
##				These in the context of, as above, meaningfull parameter estimates based on the root and observed site patterns.
if(!file.exists(paste(postDir,"filtered_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))){
	# We will cycle through each pair as defined in the unique_sigPairs object, then for each of these find out what the most common observed
	# sitePair types were.  Here we are using objects of a somewhat different shape but we can always work backward
	allPairs <- paste(unique_sigPairs$siteRef_1,unique_sigPairs$siteRef_2,sep=pairString)
	# We identify the colnames of unique_sigPairs which are for the root state probability of our models
	root_stateCols <- list("indep"=intersect(which(grepl("_indep",colnames(unique_sigPairs))), which(grepl("rootP_",colnames(unique_sigPairs)))) ,
							"dep"=intersect(which(grepl("_dep",colnames(unique_sigPairs))), which(grepl("rootP_",colnames(unique_sigPairs)))))
	# This matrix should help to identify which routes are required to get from a root to tip point
	pairedStates <- c("00","01","10","11")
	dep_routeMatrix <- matrix(c(NA,"q21","q31",paste("01","10",sep=pairString),
							"q12",NA,paste("00","11",sep=pairString),"q42",
							"q13",paste("00","11",sep=pairString),NA,"q43",
							paste("01","10",sep=pairString),"q24","q34",NA),nrow=4,dimnames=list(pairedStates, pairedStates))
	indep_routeMatrix <- matrix(c(NA,"beta2","beta1",paste("01","10",sep=pairString),
							"alpha2",NA,paste("00","11",sep=pairString),"beta1",
							"alpha1",paste("00","11",sep=pairString),NA,"beta2",
							paste("01","10",sep=pairString),"alpha1","alpha2",NA),nrow=4,dimnames=list(pairedStates, pairedStates))
	##### THIS MAY BE TOO CONSERVATIVE WITH RESPECT TO SITES INVOLVED IN ANTAGONISTIC CORRELATED EVOLUTION ##########
	if(use_parBounds){
				## HOLY SHIT IT WORKS, This method will find those sites that do not contain parameter estimates in parameter space belonging to 
				## Important parts of their phylogenetic trajectory from root to tip.
				## BUT From manual inspection is has not grabbed all the "best" sets, perhaps I need to relax my call to all, and look for not more than 1 or 2...?
				workOutput <- sapply(allPairs, function(x){
									# We use our tmpSites information to go find the aln colnames position to get the binary states from traitMat
									tmpSites <- strsplit(x,pairString)[[1]] 
									if(length(tmpSites) == 2){
										#print(x)
										
										tmpRow <- intersect(which(unique_sigPairs$siteRef_1 == tmpSites[1]),which(unique_sigPairs$siteRef_2 == tmpSites[2]))
										tmpPatterns <- traitMat[,c(which(posNames == strsplit(tmpSites[1],"::")[[1]][1]),which(posNames == strsplit(tmpSites[2],"::")[[1]][1]))]
										tmpStates <- table(apply(tmpPatterns,MARGIN = 1,function(i){ paste(as.character(i),collapse="") }))
										# We now ask which of the root states was most likely for each of the "indep" and "dep" models
										tmpRoots <- sapply(c("dep","indep"),function(y){
														tmpReturn <- colnames(unique_sigPairs)[root_stateCols[[y]][which.max(unique_sigPairs[tmpRow,root_stateCols[[y]]])]]
														tmpReturn <- substr(sub("rootP_","",tmpReturn),1,2)
														return( tmpReturn )})
										# We now look to see if the tmpStates are dominated by the 11  OR 01, 10 states, We ignore the 00 states
										# since a no response is the "default" situation.  SHOULD I WEIGHT THE 11 by 2???
										for(thisState in pairedStates){ 
											if(!is.element(thisState,names(tmpStates))){
												tmpVal <- 0
												names(tmpVal) <- thisState
												tmpStates <- c(tmpStates, tmpVal)
											} 
										}
										tmpFinal <- if(tmpStates[which(names(tmpStates) == "11")] * 2 > sum(tmpStates[which(names(tmpStates) == "01")], tmpStates[which(names(tmpStates) == "10")])){
														c("00","11")
													} else if(tmpStates[which(names(tmpStates) == "11")] * 2 < sum(tmpStates[which(names(tmpStates) == "01")], tmpStates[which(names(tmpStates) == "10")])) {
														c("10","01")
													} else {
														# If the two are equal then it's not distinguishable and we will not bias.
														print(paste("Can't distinguish the most important state pair for",x))
														return( NULL )
													}
										# Knowing the root and tip parameters we'll find the parameters that should not be zero values			
										importantParameters <- lapply(c("dep","indep"),function(y){
																	thisFinal <- tmpFinal[which(sapply(tmpFinal,function(x){ !is.element(x,tmpRoots[y]) }))]
																	tmpRoutes <- expand.grid(tmpRoots[y],thisFinal,stringsAsFactors=FALSE)
																	return( unique(apply(tmpRoutes,MARGIN = 1,function(i){
																			
																			tmpPath <- as.character(eval(as.name(paste(y,"_routeMatrix",sep="")))[unlist(strsplit(unlist(i[1]),pairString)),unlist(i[2])])
																			#print(str(tmpPath))
																			#print(x)
																			#print(y)
																			
																			if(any(is.element(strsplit(tmpPath,pairString)[[1]], pairedStates))){
																				tmpPath <- unlist(lapply(list(c(i[1],tmpPath),c(tmpPath,i[2])),function(xyz){
																						tmpPoints <- unlist(xyz)
																						return( eval(as.name(paste(y,"_routeMatrix",sep="")))[unlist(strsplit(unlist(tmpPoints[1]),pairString)),unlist(strsplit(unlist(tmpPoints[2]),pairString))] )
																				}))
																			}
																			return(tmpPath)
																		}))) })
																	
										# Now we review if any of the important parameters are at the boundary space
										#return( any(is.element(unique_sigPairs[tmpRow,unlist(importantParameters)],c(0,1000))) )
										# The above works if I have rate values in the standard output, otherwise we use the pairedTests information.
										return( any(sapply(unlist(importantParameters), function(thisCol){
												tmpValues <- c(mean(pairedTests[[tmpSites[2]]][[tmpSites[1]]][,thisCol]),sd(pairedTests[[tmpSites[2]]][[tmpSites[1]]][,thisCol]))
												#if(any(sapply(tmpValues,is.na))){ print( x )}
												return( tmpValues[1] - tmpValues[2] <= 0 || tmpValues[1] + tmpValues[2] >= 1000 )
												})) )
									} else { return( NULL )}
								})
				noBound_pars_keepPairs <- unlist(lapply(1:length(workOutput),function(x){  
											if(!is.na(workOutput[[x]]) && !is.null(workOutput[[x]])){
												#print(x)
												if(!workOutput[[x]]){ return(names(workOutput)[x]) } 
											} }))
	} else {
		noBound_pars_keepPairs <- allPairs
	}
	##############################################################################################
											
	## Here we further subset our results by looking at those pairs which have a distribution that follows some form of chi2 distribution with df = {1,4}
	## We will jacknife our testStatistics and then ask if there is some fit, and we will accept all distribution fit returns as realised return values
	## With this we will take the mean of the naive df (4) and the realised df, in order to create a new p-Value.
	ks_keepPairs_range <- lapply(noBound_pars_keepPairs,function(x){
						# We extract the two sites that were tested
						tmpSites <- strsplit(x,pairString)[[1]]
						# This is a null return, which if not overwritten, will be what is sent back!
						tmpReturn <- matrix(FALSE,nrow=length(dfRange),ncol=tmp_numReps,dimnames=list(paste("chi2_df", dfRange,sep=""),NULL))
						if(length(tmpSites) == 2){
							# In the event we did not complete a particular run we'll simply return some blank data
							# but not fail out the analysis... these check that an instance exists to be analysed.
							if(is.element(tmpSites[2],names(pairedTests))){
								if(is.element(tmpSites[1],names(pairedTests[[tmpSites[2]]]))){
								
									tmpInfo <- pairedTests[[tmpSites[2]]][[tmpSites[1]]]
									tmpRows <- sample(1:nrow(tmpInfo),tmp_numSamples * tmp_numReps,replace=FALSE)
									# We now perform replicate sub-analyses of our jacknifed data set
									tmpReturn <- sapply(1:tmp_numReps,function(thisRep){
										theseRows <- tmpRows[(1+(tmp_numSamples * (thisRep-1))):(tmp_numSamples * thisRep)]
										# This returns a logical indicating if we fail to reject the null hypothesis that this distribution is drawn 
										# from the pchisq distribution with df = 4.  i.e. when TRUE, this is a pair to keep because we fail to reject the null.
										return(sapply(dfRange,function(y){ return( ks.test(x = tmpInfo[theseRows,"testStat"],
												"pchisq" ,df = y )$p.value > ourDist_pAccept ) }))
									})
									rownames(tmpReturn) <- paste("chi2_df", dfRange,sep="")
									#return(tmpReturn)
								}# else {
									#return( NULL )
								#}
							}#else {
								#return( NULL )
							#}							
						} #else {
							#tmpReturn <- matrix(FALSE,nrow=length(dfRange),ncol=tmp_numReps,dimnames=list(paste("chi2_df", dfRange,sep=""),NULL))
							#return( tmpReturn )
						#}
						return( tmpReturn )
		})#, simplify = "list")
	
	ad_keepPairs_range <- lapply(noBound_pars_keepPairs,function(x){
						# We extract the two sites that were tested
						tmpSites <- strsplit(x,pairString)[[1]]
						# This is a null return, which if not overwritten, will be what is sent back!
						tmpReturn <- matrix(FALSE,nrow=length(dfRange),ncol=tmp_numReps,dimnames=list(paste("chi2_df", dfRange,sep=""),NULL))
						if(length(tmpSites) == 2){
							# In the event we did not complete a particular run we'll simply return some blank data
							# but not fail out the analysis... these check that an instance exists to be analysed.
							if(is.element(tmpSites[2],names(pairedTests))){
								if(is.element(tmpSites[1],names(pairedTests[[tmpSites[2]]]))){
								
									tmpInfo <- pairedTests[[tmpSites[2]]][[tmpSites[1]]]
									tmpRows <- sample(1:nrow(tmpInfo),tmp_numSamples * tmp_numReps,replace=FALSE)
									# We now perform replicate sub-analyses of our jacknifed data set
									tmpReturn <- sapply(1:tmp_numReps,function(thisRep){
										theseRows <- tmpRows[(1+(tmp_numSamples * (thisRep-1))):(tmp_numSamples * thisRep)]
										# This returns a logical indicating if we fail to reject the null hypothesis that this distribution is drawn 
										# from the pchisq distribution with df = 4.  i.e. when TRUE, this is a pair to keep because we fail to reject the null.
										return(sapply(dfRange,function(y){ return( ad.test(x = tmpInfo[theseRows,"testStat"],
												distr.fun = pchisq ,df = y )$p.value > ourDist_pAccept ) }))
									})
									rownames(tmpReturn) <- paste("chi2_df", dfRange,sep="")
								} #else {
									#return( NULL )
								#}
							} #else {
								#return( NULL )
							#}							
						} #else {
							#tmpReturn <- matrix(FALSE,nrow=length(dfRange),ncol=tmp_numReps,dimnames=list(paste("chi2_df", dfRange,sep=""),NULL))
							#return( tmpReturn )
						#}
						return(tmpReturn)
		})#, simplify = "list")
	# We now assign names to these lists so we can reference to them easily in the future
	names(ks_keepPairs_range) <- names(ad_keepPairs_range ) <- noBound_pars_keepPairs
	
	# Now we find which of our keepPairs has a chi2 fit and what the df associated with that fit are
	### NOTE: We will not be adjusting the p-value as our initial acceptance criteria has such a low q-value acceptance that testing on the
	### distribution with lower df will only mean a more significant p-value, and we are not interested in reporting on these values but rather
	### identifying pairs with support.
	final_keepPairs <- unlist(sapply(noBound_pars_keepPairs,function(x){
							tmp_ksInfo <- apply(ks_keepPairs_range[[x]], MARGIN = 1,function(y){
												sum(y) >= ceiling(tmp_numReps / 2) })
							tmp_adInfo <- apply(ad_keepPairs_range[[x]], MARGIN = 1,function(y){
												sum(y) >= ceiling(tmp_numReps / 2) })
							if(is.element(TRUE, tmp_ksInfo) && is.element(TRUE, tmp_adInfo)){
								tmpReturn <- intersect(which(tmp_ksInfo),which(tmp_adInfo))
								# In the event there is no intersection we return the minimum of these two piece of information
								if(length(tmpReturn) == 0){
									return( min(c(which(tmp_ksInfo),which(tmp_adInfo))) )
								} else {
									# If there is an intersection we only want to know the "best" chi2 df matching.
									return( max(tmpReturn) ) 
								}
							} else {
								return( NULL )
							} }, simplify="list"))
							# We now find if there is at least one row with a majority of TRUE values
							
							
	# Now let's plot out final_keepPairs and asses the performance of our test
	tmp_numPlots <- length(final_keepPairs)
	tmp_numCols <- ceiling(sqrt(tmp_numPlots))
	tmp_numRows <- ceiling(tmp_numPlots/tmp_numCols)+1
	tmp_numBreaks <- 100
	tmp_numPoints <- length(pairedTests[[1]][[1]]$testStat)
	all_keepPairs <- matrix(unlist(strsplit(names(final_keepPairs),pairString)),nrow=2)
	pdf(paste(postDir,"pairedTests_ourAccept", ourAccept,"_testStat_noBound_pars_keepPairs.pdf",sep=""),height=4 * tmp_numRows,width= 6 * tmp_numCols)
	par(mfrow=c(tmp_numRows, tmp_numCols))
	# We plot some chisq distributions with the same number of breaks as our testStat data to compare
	for(x in 1: tmp_numCols){
		hist(rchisq(tmp_numPoints,4),breaks= tmp_numBreaks,main=paste("NULL",x,sep=" "),cex.main = 1.5)
	}
	for(thisPlot in 1:ncol(all_keepPairs)){
		hist(pairedTests[[all_keepPairs[2, thisPlot]]][[all_keepPairs[1, thisPlot]]]$testStat,breaks= tmp_numBreaks,main=paste(all_keepPairs[2, thisPlot],all_keepPairs[1, thisPlot],sep="  and  "),cex.main=1.5)
	}
	dev.off()
	
	# From this filter we want to keep those pairs named in final_keepPairs
	unique_sigPairs <- unique_sigPairs[sapply(names(final_keepPairs),function(x){ 
												intersect(which(unique_sigPairs$siteRef_1 == strsplit(x,pairString)[[1]][1]),
												which(unique_sigPairs$siteRef_2 == strsplit(x,pairString)[[1]][2]) ) }) ,]
	save(unique_sigPairs, file=paste(postDir,"filtered_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))
} else {
	load(paste(postDir,"filtered_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))
}
# This is simply a tool for visualising the distribution of our significance.
filtered_sigLevels <- list()
tmpAccept <- ourAccept
while(length(which(unique_sigPairs$pValue_FDR <= tmpAccept)) > 0 && !all(unique_sigPairs$pValue_FDR[which(unique_sigPairs$pValue_FDR <= tmpAccept)] == 0) ){
	filtered_sigLevels[[as.character(tmpAccept)]] <- length(which(unique_sigPairs$pValue_FDR <= tmpAccept))
	tmpAccept <- tmpAccept/10
}

################################################################################################################
################################################################################################################
############### END OF CONSIDER LOOKING AT THE P-VALUE OBTAINED AGAINST A SITE PATTERN AND #####################
################################################################################################################
################################################################################################################				

# We will remove any of these pairs where the len_siteID > 2
#unique_sigPairs <- unique_sigPairs[intersect(which(unique_sigPairs$len_siteID_1 <= 2),which(unique_sigPairs$len_siteID_2 <= 2)),]



## THIS CALL FRIGHTENS ME BY WHAT IT SUGGETS, lpd3 is interacting with everything! ##
unique_sigPairs[order(unique_sigPairs$pValue_FDR),c(siteRef_colSet,"pValue_FDR")]

### DO WE SET A HIGHER STRING DATABASE REQUIREMENT FOR SIGNIFICANCE?  LETS TRY CALLING IT AGAIN BUT THE MINIMUM SET TO 600, 800??, 900??
######### HERE #############
# I have manually gone and done this, it can be perfected later.... it is for min_intScore == 600 (to allow parC 2006 and gyrA 248 to be in the list - medium score)
if(file.exists(paste(postDir, "filtered_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))){
	load(paste(postDir, "filtered_uniqueSigPairs_ourAccept_",ourAccept,".RData",sep=""))
}
#save(unique_sigPairs, file=paste(postDir, "filtered_sigPairs_mannual_min_intScore600.RData",sep=""))


############## THIS USES CANNONICAL RESISTANCE MUTATIONS AND KEEPS PAIRS THAT CAN TRACE THEMSELVES TO HAVE SOME CONNECTION ##############



#### UPDATE WITH MY METHOD FOR DETERMINING IF A MUTATION IS Synonymous or Non-synonymous then re-colour.
# Here we visualise networks of mutations, we have not yet made any specification about expected genes...
# We will also now include a plot of the network of correlated evolution that is proposed by our keepPairs$ppInt_vs_sngOverlap
networkCols <- brewer.pal(4, name = "Set1")
# We start by creating the data.frame of info we want to plot
# We can set a qValue threshold if we want
max_qValue <- 1e-4
keepPairs_network <- unique_sigPairs[which(unique_sigPairs$pValue_FDR <= max_qValue),c(siteRef_colSet,"pValue_FDR")]													
# We can also require that gene be characterised ones, to do this we remove all pairs with GI: or ENTREZ geneID information
limit_toCharacterised <- FALSE
if(limit_toCharacterised){ 
	keepPairs_network <- keepPairs_network[which(apply(keepPairs_network,MARGIN = 1,function(x){ 
													all(sapply(c("GI:","ENTREZ"), function(y){ 
																all(!grepl(y,unlist(x))) })) })),]
}
# Now we adjust our pValue_FDR and "fix" any values which are 0, since adjacency matrix will not handle these normally
# We will set all values of 0 to be 10% smaller than the smalles value, we also handle if the value of NA exists and we ignore it.
if (length(intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))) > 0 && is.element(0,keepPairs_network$pValue_FDR) ){
	min_qVal <- min(keepPairs_network$pValue_FDR[intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))])
	keepPairs_network$pValue_FDR[which(keepPairs_network$pValue_FDR == 0)] <- as.numeric(min_qVal) * 0.1
} else if (all(keepPairs_network$pValue_FDR == 0, na.rm = TRUE)) {
	# This means that our qValues are all equal to zero so we assign some bening small value
	keepPairs_network$pValue_FDR[which(!is.na(keepPairs_network$pValue_FDR))] <- 2e-16
}
# In order to use the iGraph adjacency matrix function with weighted edges we need to create this matrix
# We will extract the siteRef information from rows which have qValues (so are not NA)
adjVertices <- unique(as.vector(sapply(which(grepl("siteRef_",colnames(keepPairs_network))),function(j){
						sapply(which(!is.na(keepPairs_network$pValue_FDR)),function(i){ 
							return( keepPairs_network[i,j] ) 
						}) })))						
adjMatrix <- matrix(0,ncol=length(adjVertices),nrow=length(adjVertices), dimnames=list(adjVertices, adjVertices))
# Now we go through and extract qValue information from rows which are not NA (meaning this data had no instance)
for (thisRow in which(!is.na(keepPairs_network$pValue_FDR))){
	tmpVal <- as.numeric(keepPairs_network$pValue_FDR[thisRow])
	adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]]] <- tmpVal
	adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]]] <- tmpVal
}
aPlot <- graph.adjacency(adjMatrix, mode="undirected", weighted=TRUE)
# In order to view these edges properly we will exxagerate the edge weights
# We use the iGraph function of E()
E(aPlot)$weight <- -log10(E(aPlot)$weight)
pdf(file=paste(postDir,"filtered__sigPairs_max_qValue", max_qValue,"_limitCharacterised_", limit_toCharacterised,"_network.pdf",sep=""),width=45,height=45)
# This layout method of plotting is from Stephane's example of the iGraph graph.adjacency script see: "iGraph_network_Stephane.R"
# Many of these plotting options are graphical parameters of igraph's plot.igraph

# Here we create a vector to identify how we will colour our vertices, this uses information about if the gene holding a mutation
# is one our our expected genes
vertCols <- sapply(adjVertices,function(thisVertex){
	thisReturn <- NULL
	# Is the siteRef gene information one of the cannonical resistance genes?  Also is this a synonymous or non-synonymous site mutation?
	if (is.element(substr(thisVertex,1,gregexpr("_",thisVertex)[[1]][length(gregexpr("_",thisVertex)[[1]])]-1), testedGenes$expect)){
		if(site_isSyn[thisVertex]){
			return( networkCols[1] )
		} else {
			return( networkCols[2] )
		} 	
	} else {
		if(site_isSyn[thisVertex]){
			return( networkCols[3] )
		} else {
			return( networkCols[4] )
		}
	} })
		
# Now we will handle the edge.widths by creating a binned set of sizes, we will make this set with 6 bins
# Since this will create a vector of factors we use the factor level to represent the widths
thickenLine <- function(inNums){ 3* seq(2,11,length.out = 5)[as.factor(inNums)] }
edgeBin <- cut(E(aPlot)$weight, breaks=5); edge_plotThickness <- thickenLine(edgeBin) 
# NOTE: this is a number of graphical parameters for the plot.igraph option as found through: "http://igraph.org/r/doc/plot.common.html"
#plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= round(E(aPlot)$weight,3), vertex.size = 20, edge.label.cex=3, edge.width=edgeBin, 
#	vertex.label.cex=3, vertex.color= vertCols, asp = 0.66, rescale = .5)
plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= NULL, vertex.size = 15, edge.label.cex=1, edge.width= edge_plotThickness , 
	vertex.label.cex=3, vertex.color= vertCols, asp = 1)

legend("bottomleft",legend=levels(edgeBin),lwd= thickenLine(1:length(edgeBin)) ,title="Edge -log10(qValue)",bty="n",cex=3)
legend("topleft",legend=c("Cannonical_Gene_synMut", "Cannonical_Gene_nonsynMut","Novel_Gene_synMut","Novel_Gene_nonsynMut"),
		fill= networkCols,title="Mutation Type",bty="n",cex=3)
dev.off()


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# This is a space if we want to change the unique_sigPairs for any reason, a place to put code

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################


# This last plot is of the pairs that I will use and is just a means of visualising the alignment for the sites and site pairs being considered.
load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")	
if(length(list.files(path=workDir,pattern=paste(alnName,"_reduced.RData",sep=""),full.names=TRUE)) != 1){
	print("There was an error trying to load the traitMat object, please review")
	q(save="no")
} else {
	load(list.files(path=workDir,pattern=paste(alnName,"_reduced.RData",sep=""),full.names=TRUE))
}

kos_phenoVec <- unlist(sapply(rownames(traitMat),function(x){ 
if (is.element(x, c("PA14","PA01"))){ 
	return( "ND" )
} else {
	as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
} } ))
# We define the rows worth considering and setup some plotting information to be used downstream, we will keep the pairs that have ppInt support and double
# phnoetypic support PLUS one site pair which we are proposing to be part of a complete empirical fitness landscape we will construct - "rpoB_864__:__gyrA_259"
keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
plotSites <- unlist(strsplit(unique(unlist(unique_sigPairs[which(unique_sigPairs$pValue_FDR <= max_qValue),siteRef_colSet])),"::"))
plotSites <- plotSites[order(plotSites)]
plotLength <- length(plotSites)
numPlots <- 2
ntNames <- c("A","C","G","T","-")
phenoNames <- c("(R)","(S)")
# I'm pallate of colours from a qualitative set that by knowledge accepts the length of my set (max = 9)
# We set the phenoCols first as from visualising the set we know that blue and red are the first two choices
allCols <- brewer.pal(n= length(c(phenoNames,ntNames)),name= "Set1")
names(allCols) <- c(phenoNames,ntNames)

pdf(paste(postDir,"viewAlignments_max_qValue", max_qValue,".pdf",sep=""),width=6, height=24)
# This is the phenotype alignments for the sites of interest which I'm considering now as the "sngOverlap_vs_dblSupport" group
# This first plot is of the 
plotMat <- cbind(kos_phenoVec[keepRows],toupper(polyAln[keepRows, plotSites]))
# Now we create a blank plot space and then will simply plot coloured pch squares based on the phenotype and nt types
plot(x = NULL, y = NULL, xlim = c(1, plotLength + 1), ylim = c(1, nrow(plotMat)), ylab = "Phenotype and Nucleotide Alignment", xlab = "", axes = FALSE)	
# We use the colour sets and their names to be how we will search through this, these are not robust since the lists are manually generated.
for(thisState in names(allCols)){
	tmpPlot <- which(plotMat == thisState,arr.ind = TRUE)
	points(x = tmpPlot[,"col"], y = tmpPlot[,"row"], pch = 15, cex = 1, col = allCols[thisState])
}
# We will add text to the bottom of the plot to indicate what information is being shown
text(x = 0: plotLength + 1, y = 1, labels = c("Phenotype", plotSites) ,pos = 1, offset = 2, srt = 90)
legend("bottom",cex=1 ,names(allCols),bty="n",fill= allCols, ncol = 4)
dev.off()




# This last is used by the identification of synonymous and non-synonymous mutations
if(!file.exists(sub("rawAln/","wholeAln_index.RData",seqDir))){
	print("There was a problem with loading the wholeAln_index.RData file, please review")
	q(save="no")
} else {
	load(sub("rawAln/","wholeAln_index.RData",seqDir))
}


# This will return for use the mutations associated to the nucleotide position columns
siteMutations <- lapply(unique(unlist(unique_sigPairs[,siteRef_colSet])),whatAA)
names(siteMutations) <- unique(unlist(unique_sigPairs[,siteRef_colSet]))
save(siteMutations, file=paste(postDir,"siteMutations_max_qValue", max_qValue,".RData",sep=""))
rm(wholeAln_index)
# We save the element of unique_sigPairs we've reached here
save(unique_sigPairs,file=paste(postDir,"final_keepPairs_max_qValue", max_qValue,".RData",sep=""))



# Lastly I want to collect some information concerning my mutations so that I can present a table and 
# figure relating characteristics, we start by loading an AA reference object.
aaRef_frame <- read.csv(paste(baseDir,"Jonathan_Dench/MyResearch/Resources/AA_NT_genome/aaRef_table.csv",sep=""),header=TRUE, stringsAsFactors = FALSE)
# We'll alos nee the STRING database information
load(paste(baseDir,"PARuns/ppInteractions/ppInt_actions_links.RData",sep=""))
# These are the strings used to identify the root state columns of our unique_sigPairs matrix
# These should also be the string values which flank the actual root state
rootString <- c("rootP_","_dep")
# we now define the columns 
rootCols <- intersect(which(grepl(rootString[1],colnames(unique_sigPairs))),which(grepl(rootString[2],colnames(unique_sigPairs))))

mutSummary_frame <- NULL
for(thisRow in 1:nrow(unique_sigPairs)){
	# We extract the pair of mutations associated to this row of the unique_sigPairs object
	thisPair <- unlist(unique_sigPairs[thisRow, siteRef_colSet])
	# In case we have sites which had non-unqijue site patterns we simply consider the first instance in this assessment (this really isn't likely to happen)....
	thisPair <- sapply(thisPair, function(x){ strsplit(x,"::")[[1]][1] })
	# For each mutation in the pairs we will create a row in our data frame
	for(thisMutation in thisPair){
		# We will use the aaRef_frame object to see if the alternative AA of this mutation have different size, shape, Physiochemical group
		tmpAA <- aaRef_frame[which(is.element(aaRef_frame$X1Letter,names(siteMutations[[thisMutation]]$aa))),]
		# Sanity check that we found something in the reference frame
		if(nrow(tmpAA) <= 0){ stop(paste("Problem finding AA reference table information for", thisMutation,sep=" ")) }
		# Otherwise we define the changes of amino acids caused by these mutations
		tmpPhys <- length(unique(unlist(tmpAA$Physiochemical))) > 1
		tmpShape <- length(unique(unlist(tmpAA$Shape))) > 1
		tmpSize <- mean(unlist(lapply(1:(length(tmpAA$Size)-1),function(x){ unlist(lapply((x+1):length(tmpAA$Size),function(y){ abs(tmpAA$Size[x] - tmpAA$Size[y]) })) })))
		# Now we try to identify the correlation type for the pair of mutations, this will be done by looking at the most likely root state for the mutations
		# and comparing this to the tip state of the pairs, first we extract the root state (from the colname) associated to the maximum likelihood
		tmp_rootState <- colnames(unique_sigPairs[rootCols])[which.max(unique_sigPairs[thisRow, rootCols])]
		tmp_rootState <- strsplit(tmp_rootState, "_")[[1]][2]
		# Now we create a table of the tip states within our alignment and find the maximum which is not the root state, the which command removes instances
		# where the pair has not changed from the root state
		tmpTable <- table(apply(cbind(traitMat[,which(posNames == thisPair[1])],traitMat[,which(posNames == thisPair[2])]), MARGIN = 1, function(x){ 
								paste(x,collapse="") }))
		tmpTable <- tmpTable[which(names(tmpTable) != tmp_rootState)]			
		# Now we look if the tip states dominate the opposite of the root, this means positive correlation; otherwise this means that the changes
		# we negatively correlated.
		tmpAlt_tipState <- paste(1-as.numeric(strsplit(tmp_rootState,"")[[1]]),collapse="")
		tmp_corrType <- if(names(which.max(tmpTable)) == tmpAlt_tipState){
							"Positive"
						} else {
							"Negative"
						}
		# Now we go about getting the STRING database value for this pair of mutations				
		# We only really need the gene information not the actual site of the mutation.
		tmpGenes <- unlist(lapply(thisPair,function(x){ substr(x,1,regexpr("_",x)-1) }))
		# We now look if there are PA###### Locus_Tag's identifiable for all elements of thisPair
		tmpLoci <- unlist(sapply(tmpGenes,function(x){ 
			# We want to handle if nothing is found by returning an NA value
			tmpIndex <- which(pa01_genicDB$geneID == x)
			if(length(tmpIndex) == 0){
				return( NA )
			} else {
				return( pa01_genicDB$Locus_Tag[which(pa01_genicDB$geneID == x)] )
			} })) 		
		# Now we handle circumstance where nothing was found bu geneID reference by using a blast search of the sequence
		if(is.element(NA,tmpLoci)){
			# First we check if there is a pa01_genicDB database else we make one
			if (!file.exists(paste(blastDB,"pa01_genicDB.nsq",sep=""))){
				print(paste("Missing the blastdb for pa01_genicDB trying to create it",sep=" "))
				# We need to make certain there is a Seq_nt sequence otherwise we won't include that to our blastDB
				tmpIndex_seqNT <- which(nchar(pa01_genicDB$Seq_nt) > 0)
				write.fasta(as.list(pa01_genicDB$Seq_nt[tmpIndex_seqNT]),names=pa01_genicDB$Locus_Tag[tmpIndex_seqNT],file=paste(blastDB,"pa01_genicDB.fas",sep=""))
				system(paste(blastLoc,"makeblastdb -in ",blastDB,"pa01_genicDB.fas -dbtype nucl -input_type fasta -out ",blastDB,"pa01_genicDB",sep=""))
				unlink( paste(blastDB,"pa01_genicDB.fas",sep="") )
			}
			# Now we will blast the Seq_nt for this gene from pa14 to the pa01 db
			for(thisGene in tmpGenes[which(is.na(tmpLoci))]){
				# We hand the ENTREZ values as some instances of my siteRef's may not have the same caps values, this shouldn't happen for anything other than Entrez: values.
				tmpSeq <- pa14_genicDB$Seq_nt[which(toupper(pa14_genicDB$geneID)== toupper(thisGene))]
				tmpFile <- paste(blastDB,"tmp_",gsub("[[:punct:]]","_",thisGene),".fas",sep="")
				# If we found a sequence then we will perform a blast search for this
				write.fasta(as.list(tmpSeq),names=thisGene,file=tmpFile)
				system(paste(blastLoc,"blastn -query ", tmpFile," -db ",blastDB,"pa01_genicDB -out ", blastDB,"tmpQuery_blast.table -outfmt '6 sseqid sstart send qstart qend qlen evalue pident sseq'",sep=""))
				# Now we look if we found anything by searching is a _blast.table file exists and if it has content
				if(file.exists(paste(blastDB,"tmpQuery_blast.table",sep="")) && length(readLines(paste(blastDB,"tmpQuery_blast.table",sep=""))) > 0){
					# We read in the blast.table information and run a bit of quality control by verifying the percIdent value and eValue
					tmpFrame <- read.table(paste(blastDB,"tmpQuery_blast.table",sep=""),header=FALSE, stringsAsFactors = FALSE)
					# These colnames are based on my knowledge of how I am making the blast call
					colnames(tmpFrame) <- c("refID","refStart","refEnd","qStart","qEnd","qLen","eValue","PercIdent","refSeq")
					if(tmpFrame$eValue <= 1e-40 || tmpFrame$PercIdent >= 90){
						tmpLoci[thisGene] <- tmpFrame$refID[union(which(tmpFrame$eValue == min(tmpFrame$eValue)),which(tmpFrame$PercIdent == max(tmpFrame$PercIdent)))[1]]
					}
					# We no longer need this tmpQuery_blast.table information
					unlink( paste(blastDB,"tmpQuery_blast.table",sep="") )
				}
				# we also no longer need the tmp_thisGene.fas file
				unlink( tmpFile )	
			}
		}
		#  We check that there are no more NA values in the tmpLoci
		tmpSTRING <- if(!is.element(NA,tmpLoci)){
						# This looks for row entries in the STRING database object "ppInt_links.detailed" for entries for this pair of mutations
						tmpList <- lapply(tmpLoci, function(x){ union(which(grepl(x,ppInt_links.detailed$protein1)),
																										which(grepl(x, ppInt_links.detailed$protein2))) })
						tmpIntersect <- recursiveIntersect(tmpList)																					
						# If we find some rows of the STRING object then we return the combined score value
						if(length(tmpIntersect) > 0){ 
							mean(ppInt_links.detailed$combined_score[tmpIntersect])
						} else {
							# This means there is no interaction found in the database, thus the value will be zero
							0
						}
					} else {
						# We couldn't find a STRING database entry for at least one of these mutations and thus we have no value
						NA
					}
					
		# I also want to be able to determine what internal nucleotide is represented by a mutation, this matters when mutations are in negative sense strands
		# since my nomenclature is strictly left to right of the sequence rather than based on actual gene position			
		tmp_refInfo <- strsplit(thisMutation,"_")[[1]]
		tmp_ntPosition <- if(pa14_genicDB$Strand[which(pa14_genicDB$geneID == tmp_refInfo[1])]  == -1){
								nchar(pa14_genicDB$Seq_nt[which(pa14_genicDB$geneID == tmp_refInfo[1])]) + 1 - as.numeric(tmp_refInfo[2])
							} else {
								as.numeric(tmp_refInfo[2])
							}
		# Here we define the basal nucleotide and character states of at the site related to this mutation
		tmp_ntCharacters <- names(siteMutations[[thisMutation]]$codonTable[[siteMutations[[thisMutation]]$ntPos]])
		# Here we adjust this information so it is a list of the PA14 reference character and all others 
		tmp_ntCharacters <- list("reference" = siteMutations[[thisMutation]]$codon["PA14",siteMutations[[thisMutation]]$ntPos],
									"alternates" = tmp_ntCharacters[which(tmp_ntCharacters != siteMutations[[thisMutation]]$codon["PA14",siteMutations[[thisMutation]]$ntPos])])
		# Now we adjust the characters so that if this is a reverse strand we want the compliments
		if(pa14_genicDB$Strand[which(pa14_genicDB$geneID == tmp_refInfo[1])]  == -1){
			# If an NA has been created we change it to gap character status
			tmp_ntCharacters$reference <- comp(tmp_ntCharacters$reference)
			tmp_ntCharacters$reference[which(is.na(tmp_ntCharacters$reference))] <- "-"
			tmp_ntCharacters$alternates <- comp(tmp_ntCharacters$alternates)
			tmp_ntCharacters$alternates[which(is.na(tmp_ntCharacters$alternates))] <- "-"
		}
		
		tmp_ntCharacters <- paste(tmp_ntCharacters$reference, paste(tmp_ntCharacters$alternates,collapse="_::_"), sep = "_-->_")
		
		# Now we ought to have all our pieces of information to put together our informational data.frame's row
		mutSummary_frame <- rbind(mutSummary_frame, data.frame("Mutation"=thisMutation,
																"Partner"= thisPair[which(thisPair != thisMutation)],
																"corrType"= tmp_corrType,
																"aaTypes" = paste(sapply(1:length(siteMutations[[thisMutation]]$aa), function(x){ paste(names(siteMutations[[thisMutation]]$aa[x]),siteMutations[[thisMutation]]$aa[x],sep="_") }),collapse="_::_"),
																"aa_physChng"= tmpPhys,
																"aa_shapeChng"= tmpShape,
																"aa_sizeChng"= round(tmpSize,4),
																"nt_intPosition" = tmp_ntPosition,
																"nt_polyMorphisms" = tmp_ntCharacters,
																"intSupp_STRING"= tmpSTRING,
																"pValue_FDR"= unique_sigPairs$pValue_FDR[thisRow], stringsAsFactors=FALSE, row.names=NULL) )
	} # This closes out the cycle through mutations in the row
} # This closes out going through each row
# We now save the summary frame we've created, and write it out as well
write.csv(mutSummary_frame, file=paste(postDir,"mutSummary_max_qValue", max_qValue,"_frame.csv",sep=""),row.names=FALSE)
save(mutSummary_frame, file=paste(postDir,"mutSummary_max_qValue", max_qValue,"_frame.RData",sep=""))







# This is a figure to see how the dsitribution of q-Values falls for my different sets and analyses
# It use some objects created along the way....
# We define the number of breaks we want in our histogram

# First we identify the pairs which have single, and double synoymous mutations
isPair_nonSyn <- apply(all_sigPairs[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
		tmpPair <- sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] })
		# We check that both of these are within site_isSyn (which drops certain instances....)
		if(!all(is.element(tmpPair,names(site_isSyn)))){
			return( NA )
		# We ask if any of these are nonSynonymous
		} else if(all(!site_isSyn[tmpPair])){
			return( TRUE )
		} else {
			return( FALSE )
		}
	})

isSingle_nonSyn <- apply(all_sigPairs[,siteRef_colSet], MARGIN = 1, function(thisPair){ 
		tmpPair <- sapply(thisPair,function(x){ strsplit(x,"::")[[1]][1] })
		# We check that both of these are within site_isSyn (which drops certain instances....)
		if(!all(is.element(tmpPair,names(site_isSyn)))){
			return( NA )
		# We ask if any of these are nonSynonymous
		} else if(any(!site_isSyn[tmpPair])){
			return( TRUE )
		} else {
			return( FALSE )
		}
	})	


pdf(paste(postDir,"qValue_distPlots.pdf",sep=""),height = 16, width= 12)
par(mfcol=c(2,1))
# We grab out all the unique_sigPairs from all_sigPairs in order to have a set we can modify without hassle of replication of script runs.
tmp_allPairs <- all_sigPairs
tmp_finalPairs <- unique_sigPairs
numBreaks <- 100
# Before plotting we need to convert all the counts to be ln counts (so that our small sets are visible
# and also convert any values of zero into some alternate value than can pass through log functions
tmp_allPairs$pValue_FDR[which(tmp_allPairs$pValue_FDR == 0)] <- min(tmp_allPairs$pValue_FDR[which(tmp_allPairs$pValue_FDR != 0)]) / 100
tmp_finalPairs$pValue_FDR[which(tmp_finalPairs$pValue_FDR == 0)] <- min(tmp_allPairs$pValue_FDR)

fullHist <- hist(log10(tmp_allPairs$pValue_FDR),breaks= numBreaks, plot = FALSE)

# Now we add some new elements to the histogram using the alpha colour value to make transparency
# as well we use the breaks value from the fullHist object so we can apply to our existing
nonsynHist_isPair <- hist(log10(tmp_allPairs$pValue_FDR[which(isPair_nonSyn)]), breaks = fullHist$breaks, plot = FALSE)
nonsynHist_isSingle <- hist(log10(tmp_allPairs$pValue_FDR[which(isSingle_nonSyn)]), breaks = fullHist$breaks, plot = FALSE)

# Now we adjust the counts to be log2
for(thisPlot in c("fullHist", "nonsynHist_isPair", "nonsynHist_isSingle")){
	tmpData <- eval(as.name(thisPlot))
	# We adjust single counts to be values of 1.1 so that we can differentiate them from zero values
	tmpData$counts[which(tmpData$counts == 1)] <- 1.1
	tmpData$counts[which(tmpData$counts != 0)] <- log2(tmpData$counts[which(tmpData$counts != 0)])
	assign(thisPlot, tmpData , pos=".GlobalEnv")
}

plot(fullHist, main = "All p <= 1e-4 vs. nonSyn Pairs", ylab = "Log2 Frequency", xlab = "Log10 pValue")
plot(nonsynHist_isSingle, col = alpha("blue",0.5), add=TRUE)
plot(nonsynHist_isPair, col = alpha("red",0.5), add=TRUE)

# Now I want to add the sites which failed the chi2 and LD tests, this means I can use the difference between
# tmp_allPairs$pValue_FDR[which(isPair_nonSyn)] and final_keepPairs
fullHist_nonSyn <- hist(log10(tmp_allPairs$pValue_FDR[which(isPair_nonSyn)]),breaks= numBreaks, plot = FALSE)
hist_finalKeep <- hist(log10(tmp_finalPairs$pValue_FDR), breaks = fullHist_nonSyn$breaks, plot = FALSE)
plot(fullHist_nonSyn, main = "All nonSyn Pairs vs. Final", ylab = "Log2 Frequency", xlab = "Log10 pValue", col=alpha("red",0.5))
plot(hist_finalKeep, col = alpha("purple"), add=TRUE)
dev.off()


# This is my network plots method but here so I can fiddle with it
# We start by creating the data.frame of info we want to plot
# We can set a qValue threshold if we want
max_qValue <- 1e-9
networkCols <- brewer.pal(4, name = "Set1")
keepPairs_network <- all_sigPairs[which(all_sigPairs$pValue_FDR <= max_qValue),c(siteRef_colSet,"pValue_FDR")]													
for(thisCol in siteRef_colSet){  keepPairs_network[,thisCol] <- sapply(keepPairs_network[,thisCol], function(x){ strsplit(x,"::")[[1]][1] }) }
limit_toCharacterised <- FALSE
if(limit_toCharacterised){ 
	keepPairs_network <- keepPairs_network[which(apply(keepPairs_network,MARGIN = 1,function(x){ 
													all(sapply(c("GI:","ENTREZ"), function(y){ 
																all(!grepl(y,unlist(x))) })) })),]
}
if (length(intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))) > 0 && is.element(0,keepPairs_network$pValue_FDR) ){
	min_qVal <- min(keepPairs_network$pValue_FDR[intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))])
	keepPairs_network$pValue_FDR[which(keepPairs_network$pValue_FDR == 0)] <- as.numeric(min_qVal) * 0.1
} else if (all(keepPairs_network$pValue_FDR == 0, na.rm = TRUE)) {
	keepPairs_network$pValue_FDR[which(!is.na(keepPairs_network$pValue_FDR))] <- 2e-16
}
adjVertices <- unique(as.vector(sapply(which(grepl("siteRef_",colnames(keepPairs_network))),function(j){
						sapply(which(!is.na(keepPairs_network$pValue_FDR)),function(i){ 
							return( keepPairs_network[i,j] ) 
						}) })))
adjMatrix <- matrix(0,ncol=length(adjVertices),nrow=length(adjVertices), dimnames=list(adjVertices, adjVertices))
for (thisRow in which(!is.na(keepPairs_network$pValue_FDR))){
	tmpVal <- as.numeric(keepPairs_network$pValue_FDR[thisRow])
	adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]]] <- tmpVal
	adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]]] <- tmpVal
}
aPlot <- graph.adjacency(adjMatrix, mode="undirected", weighted=TRUE)
E(aPlot)$weight <- -log10(E(aPlot)$weight)

#pdf(file=paste(postDir,"allPairs_max_qValue", max_qValue,"_limitCharacterised_", limit_toCharacterised,"_network.pdf",sep=""),width=60,height=60)
pdf(file=paste(postDir,"final_keepPairs_nonSyn_chi2_network.pdf",sep=""),width=60,height=60)
# Here is where we can do some colouring for visualisation means
vertCols <- sapply(adjVertices,function(thisVertex){
	thisReturn <- NULL
	if (is.element(substr(thisVertex,1,gregexpr("_",thisVertex)[[1]][length(gregexpr("_",thisVertex)[[1]])]-1), testedGenes$expect)){
		if(site_isSyn[thisVertex]){
			return( networkCols[1] )
		} else {
			return( networkCols[2] )
		} 	
	} else {
		if(site_isSyn[thisVertex]){
			return( networkCols[3] )
		} else {
			return( networkCols[4] )
		}
	} })
thickenLine <- function(inNums){ 3* seq(2,11,length.out = 5)[as.factor(inNums)] }
edgeBin <- cut(E(aPlot)$weight, breaks=5); edge_plotThickness <- thickenLine(edgeBin) 
plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= NULL, vertex.size = 16, edge.label.cex=1, edge.width= edge_plotThickness , 
	vertex.label.cex=5, vertex.color= vertCols, asp = 1.5)

legend("bottomleft",legend=levels(edgeBin),lwd= thickenLine(1:length(edgeBin)) ,title="Edge -log10(qValue)",bty="n",cex=5)
legend("topleft",legend=c("Cannonical_Gene_synMut", "Cannonical_Gene_nonsynMut","Novel_Gene_synMut","Novel_Gene_nonsynMut"),
		fill= networkCols,title="Mutation Type",bty="n",cex=5)
dev.off()

#pdf(paste(postDir,"allPairs_pairCounts_max_qValue",max_qValue,".pdf",sep=""),height=8, width = 6)
pdf(paste(postDir,"final_keepPairs_nonSyn_chi2_barplots.pdf",sep=""),height=8, width = 6)
par(mfcol=c(2,1))
# Let's review how many of our pairs in intragenic and in how many pairs they are found
tmpTable <- table(unlist(keepPairs_network[,siteRef_colSet]))
barplot(table(tmpTable), xlab="Number of Pairs with Mutation", main= paste("Pair Counts for pValue <= ",max_qValue,sep=""), ylab="Count")
# Here is the count of number of intragenic pairs
isIntragenic <- apply(keepPairs_network[,siteRef_colSet], MARGIN = 1,function(x){
					strsplit(strsplit(x[1],"::")[[1]][1],"_")[[1]][1] == strsplit(strsplit(x[2],"::")[[1]][1],"_")[[1]][1]
				})
barplot(table(isIntragenic), main = paste("Pair is Intragenic for pValue <= ",max_qValue,sep=""))
dev.off()


# This is for if I want to make an example of the test statistics for sites
# Here we define a list of vectors which are the leading and paired sites of interest
testSite <- list("gyrA_248"=c("parC_2006"),"rpoB_864"=c("parE_856","dnr_364"), "gyrB_1422"=c("gyrB_1443"), 
					 "parC_679"=c("parC_685"), "parC_685"=c("parC_733","parC_712"), "dnaN_495"=c("dnaN_504"),"morA_4041"=c("morA_4083") )
# Note there must be at least one colour for each pair to be plotted from the list above
histCols <- c("red","blue","purple","orange")  
numPlots <- length(unlist(testSite))
numCols <- numRows <- ceiling(sqrt(numPlots))
numBreaks <- seq(0,max(unlist(lapply(names(testSite),function(x){ unlist(lapply(testSite[[thisPlot1]], function(y){ pairedTests[[x]][[y]]$testStat })) }))) , length.out = 100)
numPoints <- length(pairedTests[[1]][[1]]$testStat)

# Now make a plot for the test stat of all values considered and tested.
pdf(paste(postDir,"pairedTests_specificTargets.pdf",sep=""),height=4 * numRows,width= 6 * numCols )
par(mfrow=c(numCols,numRows))
colCounter <- 1
# We plot some chisq distributions with the same number of breaks as our testStat data to compare
#tmp_mainHist <- hist(rchisq(10000,4),breaks= numBreaks,main="Test Statistic Distributions", xlab = "Chi2 Distribution (df = 4)")#, ylim = c(0,numPoints/5))

for(thisPlot1 in names(testSite)){
	for(thisPlot2 in testSite[[thisPlot1]]){
		hist(pairedTests[[thisPlot1]][[thisPlot2]]$testStat,breaks= tmp_mainHist$breaks, col = if(colCounter == 1){ histCols[3] } else if(is.element(colCounter,c(2,3))){ histCols[2] } else { histCols[4] }  ,
				main = paste(thisPlot1, thisPlot2, sep = " :: "), xlab = "Test Statistics", probability = TRUE )
		# We add a distribution of chi2 density as comparison
		lines(density(rchisq(100000,4)), lwd = 4, col = "red")
		colCounter <- colCounter + 1
	}
}
dev.off()


###################################################################################################
######################################## END OF BODY ##############################################
###################################################################################################

q(save="no")




































































































































































































































































































######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
# THESE ARE OLD PIECES THAT CAN BE USED TO ASSEMBLE MY NEW ROBUST SCRIPT

# This is just a handy set of calls that will allow me to get a feel for some information concerning a mutant site



#### This will allow me to visualise if there might be a correlation between mutations and the phenotype used to recode the data.
rm(list=ls())
# Here we load our data sets that will be fed to generate the confusion table
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/susceptStrains.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln_traitMat_5perc_wholeAlnTree.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")

# Now we will use columns of interest, from the postAnalysis of corr. evolution, to consider two pair members and
# see if they are correalted to the categorical information of strain susceptability.  Without this then the correlation identified
# may not have much of a biological basis, but requires interpretation.


# This is to get the binary vector of our phenotype information (the second categorical dimension)
susceptVec <- sapply(rownames(polyAln),function(x){ !(is.element(x,susceptStrains)) })

# We get out the siteID information for our sites of interest
siteNames <- c("gyrA_248","parC_2006")
tmpSites <- unlist(sapply(siteNames,function(x){ which(colnames(polyAln) == x) }))
# This is to visualise the alignemnts
cbind(susceptVec,traitMat[,tmpSites])
# This will give me an idea of the number of instance when they are together, appart, etc...
pairVec <- sapply(1:nrow(polyAln), function(x){ return( paste(traitMat[x,tmpSites],collapse="::") ) })
table(pairVec)/nrow(polyAln)



######## This will inform of the frequency of appearance of paired mutations, if they are together or separate
library(foreach)
library(doMC)
registerDoMC(cores=8) 
options(digits=3, stringsAsFactors=FALSE)
pairs_freqTable <- foreach(thisRow = 1:nrow(sigPairs), .combine="rbind") %dopar% {
#for( thisRow in 1:nrow(sigPairs)){
	siteNames <- unlist(sigPairs[thisRow,c("siteRef_1","siteRef_2")])
	tmpSites <- unlist(sapply(siteNames,function(x){ which(colnames(polyAln) == x) }))
	pairVec <- sapply(1:nrow(polyAln), function(x){ return( paste(traitMat[x,tmpSites],collapse="::") ) })
	tmpTable <- table(pairVec)/nrow(polyAln)
	tmpReturn <- sapply(c("0::0","0::1","1::0","1::1","2"),function(x){ 
					if (length(which(grepl(x,names(tmpTable)))) > 0){
						return( sum(tmpTable[which(grepl(x,names(tmpTable)))]) )
					} else {
						return( NA )
					} })

	return( data.frame("siteRef_1"=siteNames[1],"siteRef_2"=siteNames[2],
						"Neither"= tmpReturn[1],
						"onlySite2"= tmpReturn[2],
						"onlySite1"= tmpReturn[3],
						"Both"= tmpReturn[4],
						"Twos"= tmpReturn[5]) )
}
rownames(pairs_freqTable) <- NULL
save(pairs_freqTable, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/pairs_freqTable.RData")



##### This is to get the translated series of code information and see if the site is synonymous or not. #####
library(seqinr)
library(ape)
load("/Users/Jdench/Desktop/PARuns/Sequences/PA14DBase/PA14_DataBase.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/wholeAln_index.RData")
#################################################################################################
##################################### BEGIN PRE FUNCTIONS #######################################
#################################################################################################
# This is the translate function from seqinr, but we needed to fix part of the code
# So I've copied it here and will make a simple change
translate <- function (seq, frame = 0, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE) {
    #print("My version")
    if (any(seq %in% LETTERS)) {
        seq <- tolower(seq)
    }
    # I added a line to cahnge NA back into gap characters, this is because any NA passed to s2n causes a crashout
    if (sens == "R") {
        seq <- comp(rev(seq), ambiguous = ambiguous)
        seq[is.na(seq)] <- "-"
    } 
    seqn <- s2n(seq, levels = s2c("tcag"))
    l <- 3 * ((length(seq) - frame)%/%3)
    c1 <- seq(from = frame + 1, to = frame + l, by = 3)
    tra <- 16 * seqn[c1] + 4 * seqn[c1 + 1] + seqn[c1 + 2] + 
        1
    code <- s2c(SEQINR.UTIL$CODES.NCBI$CODES[numcode])
    result <- code[tra]
    result[is.na(result)] <- NAstring
    if (ambiguous) {
        toCheck <- which(result == NAstring)
        for (i in toCheck) {
            codon <- seq[c1[i]:(c1[i] + 2)]
            allcodons <- as.vector(outer(as.vector(outer(amb(codon[1]), 
                amb(codon[2]), paste, sep = "")), amb(codon[3]), 
                paste, sep = ""))
            allaminoacids <- sapply(allcodons, function(x) translate(s2c(x), 
                numcode = numcode, ambiguous = FALSE))
            if (all(allaminoacids == allaminoacids[1])) 
                result[i] <- allaminoacids[1]
        }
    }
    return(result)
}

# This will identify if all the objects in a vector fed in as v, are identical
all_identical <- function(v) all(sapply( as.list(v[-1]),FUN=function(z) {identical(z, v[1])}))

#################################################################################################
#####################################  END PRE FUNCTIONS ########################################
#################################################################################################

siteNames <- c("obgE_1084","parC_712","GI:116049432_723")
# Here we load an alignment so that it can be translated, the alignment loaded is based on the Locus_Tag of the gene
siteLocus <- unlist(sapply(siteNames,function(x){
	x <- strsplit(x,"_")[[1]][1]
	if (grepl("GI:",x)){
		return( as.character(dBase$Locus_Tag[which(grepl(x,dBase$Molecule_Xrefs))]) )
	} else {
		return( as.character(dBase$Locus_Tag[which(dBase$Gene == x)]) )
	} 
}))
# This is where we can load alignments from
alnDir <- "/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/rawAln/"
responseVec <- NULL

for(thisSite in 1:length(siteLocus)){
	tmpAln <- read.dna(paste(alnDir,siteLocus[thisSite],"_seqAln.fas",sep=""),format="fasta",as.character=TRUE)
	# These colnames will help us keep track of the Ref and Insert _ntPos of sites as we clean.
	colnames(tmpAln) <- sapply(which(wholeAln_index$Locus_Tag == siteLocus[thisSite]), function(x){ 
		# we create the colnames using the Ref or Insert _ntPos information from our index file
		paste(wholeAln_index$Gene[x], if (is.na(wholeAln_index$Ref_ntPos[x])){ wholeAln_index$Insert_ntPos[x] } else { wholeAln_index$Ref_ntPos[x] }, sep="_") }) 
	
	# We make certain there are only the same strains as contained within polyAln
	tmpAln <- tmpAln[unlist(sapply(rownames(polyAln),function(x){ which(rownames(tmpAln) == x) })),]
	
	# This is a quality/sanity check to be reviewed manually through the output file, it checks if the alignment was a multiple of 3
	if (length(tmpAln)[[1]] %% 3 != 0){ 
		print(rep("-",25))
		print(paste("The loaded alignment for ",siteLocus[thisSite]," - ",strsplit(siteLocus[thisSite],"_")[[1]][1]," was not a multiple of three, mannual review",sep=""))
		print(rep("-",25))
	}
	
	tmpTrans <- list()
	# We will identify any positions in this gene which our index defines as being insert positions
	tmp_indPos <- which(wholeAln_index$Locus_Tag == siteLocus[thisSite])
	tmpInserts <- unlist(sapply(tmp_indPos, function(x){ if(!is.na(wholeAln_index$Insert_ntPos[x])){ return( which(tmp_indPos == x) ) } }))
	for (tmpRow in 1:nrow(tmpAln)){	
		# We will remove any gap characters that belong to insert positions, these can be identified from the index
		tmpNucleotides <- tmpAln[tmpRow,]
		# If there are gaps at the insert positions in this sequence we will remove them before translating this information
		if (length(which(tmpNucleotides[tmpInserts] == "-")) > 0){ tmpNucleotides <- tmpNucleotides[-tmpInserts[which(tmpNucleotides[tmpInserts] == "-")]] }
		tmpTrans[[rownames(tmpAln)[tmpRow]]] <- translate(tmpNucleotides,sens= if(dBase$Strand[which(dBase$Locus_Tag == siteLocus[thisSite])] == 1){"F"} else {"R"}  )
	}
	# Now we look at the amino acid which relates to the ntPos, we need to know the strand to understand if we read forward or backward
	# for our divison of codons
	tmp_ntPos <- if (dBase$Strand[which(dBase$Locus_Tag == siteLocus[thisSite])] == 1){ 
		as.numeric(strsplit(siteNames[thisSite],"_")[[1]][2])
	} else {
		ncol(tmpAln) + 1 - as.numeric(strsplit(siteNames[thisSite],"_")[[1]][2])
	}
	#print(paste(" This site is actually the ",tmp_ntPos,"th nucleotide when read ",sep=""))
	tmp_aaPos <- ceiling( tmp_ntPos / 3 )
	print(paste(" This site in ",siteNames[thisSite]," is actually the ",tmp_aaPos,"th amino acid when read ",sep=""))
	
	# If there is more than one character than this was a non-synonymous change
	if (length(table(sapply(tmpTrans,function(x){x[tmp_aaPos] }))) == 1){
		print(paste(" Change ",siteNames[thisSite]," is synonymous ",sep=""))
		responseVec <- c(responseVec,"Synonymous")
	} else if (length(table(sapply(tmpTrans,function(x){x[tmp_aaPos] }))) > 1) {
		print(paste(" Change ",siteNames[thisSite]," is non-synonymous ",sep=""))
		responseVec <- c(responseVec,"Non_Synonymous")
	}
}







# This will allow me to identify, as per the reference database (dBase of PA14) which genes are nearby to genes 
# which hold mutations of interest.
rm(list=ls())
library(ape)
library(foreach)
library(doMC)
registerDoMC(cores=8)
library(pegas)  # this has the nuc.div pi calculating function (Nei's pi)
options(digits=3, stringsAsFactors=FALSE)
load("/Users/Jdench/Desktop/PARuns/Sequences/PA14DBase/PA14_DataBase.RData")
load('~/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/CompareSets_polyAln_0.8/q_1e-06/sigPairs_FALSE.RData')
load("/Users/Jdench/Desktop/PARuns/Sequences/PA14DBase/PA14_DataBase.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/wholeAln_index.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/susceptStrains.RData")

alnDir <- "/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/rawAln/"
workDir <- "/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/codeML/"

tTree <- read.tree("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/wholeAln_reduced.tree")

# All we need to do is extract the list of gene's from the sireRef_1, _2 information
allGenes <- unique(unlist(sapply(unlist(unique_sigPairs[,c("siteRef_1","siteRef_2")]),function(x){
	return( strsplit(x,"_")[[1]][1] ) })))
allLoci <- sapply(allGenes, function(x){
		as.character(if (grepl("GI:",x)){
			dBase$Locus_Tag[which(grepl(x,dBase$Molecule_Xrefs))]
		} else {
			dBase$Locus_Tag[which(dBase$Gene == x)]
		}) })

# Now we find the dBase$Locus_Tag for each of these genes, but confirming there is a _seqAln.fas assocaited file,
# accounting for GI information, and will look for the genes within +/- 50 genes of this Locus_Tag (regionSize), these will 
# be passed to PAML's CODEML in order to analyse dN/dS and see if there are regions with a signal that might show 
# that hitchiking has caused a association with a phenotype NOTE: This assumes there even is an association!
allFiles <- list.files(path=alnDir,pattern="_seqAln.fas")


if(!file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/nucDiv_lociList.RData")){
	regionSize <- 51
	lociList <- lapply(allLoci, function(x){
		# We want a set that has regionSize of loci held within it, so we look for which element of allFiles is tmpLoci
		# and grab a balanced number of flanking loci, adjusted for ends
		tmpInd <- which(grepl(x,allFiles))
		tmpRange <- c(tmpInd - (regionSize-1)/2, tmpInd + (regionSize-1)/2)
		if (tmpRange[1] < 1){
			tmpRange[2] <- (tmpRange[2] - tmpRange[1]) + 1
			tmpRange[1] <- 1
		} else if (tmpRange[2] > length(allFiles)){
			tmpRange[1] <- tmpRange[1] - (tmpRange[2] - length(allFiles))
			tmpRange[2] <- length(allFiles)
		}
		# We return a set of the loci that are 5 +/- from tmpLoci
		return( allFiles[tmpRange[1]:tmpRange[2]] ) })
		
	# This list will now allow us to relate, for each allGenes of interest, which PAML CODEML information needs to be plotted
	# in order to look for patterns of chaning pi (dN/dS), BUT now we need to call CODEML and get the dN/dS information
	nucDiv <- foreach(thisFile = unique(unlist(lociList)), .combine="rbind") %dopar%{ 
		# we load as a DNAbin type object to be handled by the nuc.div call
		tmpAln <- read.dna(file=paste(alnDir,allFiles[which(grepl(thisFile,allFiles))],sep=""), format="fasta")
		# These colnames will help us keep track of the Ref and Insert _ntPos of sites as we clean.
		colnames(tmpAln) <- sapply(which(wholeAln_index$Locus_Tag == sub("_seqAln.fas","",thisFile)), function(x){ 
			# we create the colnames using the Ref or Insert _ntPos information from our index file
			paste(wholeAln_index$Gene[x], if (is.na(wholeAln_index$Ref_ntPos[x])){ wholeAln_index$Insert_ntPos[x] } else { wholeAln_index$Ref_ntPos[x] }, sep="_") }) 
		# we keep only those sequence elements that are in tTree (these were the ones analysed and worth considering...
		tmpAln <- tmpAln[unlist(sapply(tTree$tip.label, function(x){ which(rownames(tmpAln) == x) })),]
		# To make this meaningful we also now remove gappy columns or sequences (less than 90% non-gap characters)
		# This was done since nuc.div was not able to calculate pi on files which had these very gappy sequences
		# Here we remove any sequences that are almost entirely gaps, then remove columns which are gappy in the moderate
		# quality sequences, and then make a final pass to remove any sequences that are still largely gaps
		tmpAln <- tmpAln[unlist(sapply(1:nrow(tmpAln), function(x){
			tmpTable <-  table(as.character(tmpAln[x,]))
			if (sum(tmpTable[which(names(tmpTable) != "-")])/ncol(tmpAln) >= 0.2) { 
				return( x ) 
			} })),]
		
		tmpAln <- tmpAln[,unlist(sapply(1:ncol(tmpAln), function(x){
			tmpTable <-  table(as.character(tmpAln[,x]))
			if (sum(tmpTable[which(names(tmpTable) != "-")])/nrow(tmpAln) >= 0.9) { 
				return( x ) 
			} }))]
			
		tmpAln <- tmpAln[unlist(sapply(1:nrow(tmpAln), function(x){
			tmpTable <-  table(as.character(tmpAln[x,]))
			if (sum(tmpTable[which(names(tmpTable) != "-")])/ncol(tmpAln) >= 0.9) { 
				return( x ) 
			} })),]
		# This will allow us to define a separate pi value for each phenotype
		tmp_susceptStrains <- unlist(sapply(susceptStrains,function(x){ which(rownames(tmpAln) == x) }))
			
		# We will write out this tidied tmpAln file for use downstream
		write.dna(tmpAln,file=paste(workDir,"cleanAln/",thisFile,sep=""),format="fasta", colsep="", nbcol=15)
		
		return( data.frame("Loci"=thisFile,"refGenome_start"=dBase$Start[which(dBase$Locus_Tag == sub("_seqAln.fas","",thisFile))],"refGenome_end"=dBase$End[which(dBase$Locus_Tag == sub("_seqAln.fas","",thisFile))],
							"pi_Suscept"=if(nrow(tmpAln) > 1 && length(tmp_susceptStrains) > 1){ nuc.div(tmpAln[tmp_susceptStrains,]) } else { 99999 },
							"pi_NotSuscept"=if(nrow(tmpAln) > 1 && length(setdiff(1:nrow(tmpAln), tmp_susceptStrains)) > 1){ nuc.div(tmpAln[-tmp_susceptStrains,]) } else { 99999 }) )
	}
	
	
	# I also want to look at intragenic pi, I will plot genes of interest with +/- 1 gene around them
	regionSize = 3
	lociList_intraGene <- lapply(allLoci, function(x){
		# We want a set that has regionSize of loci held within it, so we look for which element of allFiles is tmpLoci
		# and grab a balanced number of flanking loci, adjusted for ends
		tmpInd <- which(grepl(x,allFiles))
		tmpRange <- c(tmpInd - (regionSize-1)/2, tmpInd + (regionSize-1)/2)
		if (tmpRange[1] < 1){
			tmpRange[2] <- (tmpRange[2] - tmpRange[1]) + 1
			tmpRange[1] <- 1
		} else if (tmpRange[2] > length(allFiles)){
			tmpRange[1] <- tmpRange[1] - (tmpRange[2] - length(allFiles))
			tmpRange[2] <- length(allFiles)
		}
		# We return a set of the loci that are 5 +/- from tmpLoci
		return( allFiles[tmpRange[1]:tmpRange[2]] ) })
	
	nucDiv_intraGene <- NULL
	for (thisFile in unique(unlist(lociList_intraGene))) { 
		# we load as a DNAbin type object to be handled by the nuc.div call
		tmpAln <- read.dna(file=paste(alnDir,allFiles[which(grepl(thisFile,allFiles))],sep=""), format="fasta")
		# These colnames will help us keep track of the Ref and Insert _ntPos of sites as we clean.
		colnames(tmpAln) <- sapply(which(wholeAln_index$Locus_Tag == sub("_seqAln.fas","",thisFile)), function(x){ 
			# we create the colnames using the Ref or Insert _ntPos information from our index file
			paste(wholeAln_index$Gene[x], if (is.na(wholeAln_index$Ref_ntPos[x])){ wholeAln_index$Insert_ntPos[x] } else { wholeAln_index$Ref_ntPos[x] }, sep="_") }) 
		# we keep only those sequence elements that are in tTree (these were the ones analysed and worth considering...
		tmpAln <- tmpAln[unlist(sapply(tTree$tip.label, function(x){ which(rownames(tmpAln) == x) })),]
		# To make this meaningful we also now remove gappy columns or sequences (less than 90% non-gap characters)
		# This was done since nuc.div was not able to calculate pi on files which had these very gappy sequences
		# Here we remove any sequences that are almost entirely gaps, then remove columns which are gappy in the moderate
		# quality sequences, and then make a final pass to remove any sequences that are still largely gaps
		tmpAln <- tmpAln[unlist(sapply(1:nrow(tmpAln), function(x){
			tmpTable <-  table(as.character(tmpAln[x,]))
			if (sum(tmpTable[which(names(tmpTable) != "-")])/ncol(tmpAln) >= 0.2) { 
				return( x ) 
			} })),]
		
		keepCols <- unlist(sapply(1:ncol(tmpAln), function(x){
			tmpTable <-  table(as.character(tmpAln[,x]))
			if (sum(tmpTable[which(names(tmpTable) != "-")])/nrow(tmpAln) >= 0.9) { 
				return( x ) 
			} }))
		
		# Then we adjust tmpAln	
		tmpAln <- tmpAln[,keepCols]
			
		tmpAln <- tmpAln[unlist(sapply(1:nrow(tmpAln), function(x){
			tmpTable <-  table(as.character(tmpAln[x,]))
			if (sum(tmpTable[which(names(tmpTable) != "-")])/ncol(tmpAln) >= 0.9) { 
				return( x ) 
			} })),]
		
		# This will allow us to define a separate pi value for each phenotype
		tmp_susceptStrains <- unlist(sapply(susceptStrains,function(x){ which(rownames(tmpAln) == x) }))
		
		# We will write out this tidied tmpAln file for use downstream
		write.dna(tmpAln,file=paste(workDir,"cleanAln/",thisFile,sep=""),format="fasta", colsep="", nbcol=15)
		# We store the cleaned ntPos since we have cleaned our tmpAln so this might not truly relate to an original ntPosition but should give us a relatively close value, 
		# further we can look at if the length of the reference gene and what has actually be analysed to get a feel for how much might have been removed.  We can use our adj_ntPos for sites of interest
		nucDiv_intraGene <- rbind(nucDiv_intraGene, data.frame("Loci"=thisFile,"colnames_cleanedAln"= colnames(tmpAln), 
																"refGenome_start"=dBase$Start[which(dBase$Locus_Tag == sub("_seqAln.fas","",thisFile))],
																"refGenome_end"=dBase$End[which(dBase$Locus_Tag == sub("_seqAln.fas","",thisFile))],
																"pi_Suscept"=if(nrow(tmpAln) > 1 && length(tmp_susceptStrains) > 1){ 
																		sapply(1:ncol(tmpAln),function(x){
																			tmpCols <- intersect(tmp_susceptStrains,intersect(which(as.character(tmpAln[,x]) != "-"),which(tolower(as.character(tmpAln[,x])) != "n"))) 
																			nuc.div(tmpAln[tmpCols,x], pairwise.deletion = TRUE) }) 
																	} else { 
																		99999 
																	},
																"pi_NotSuscept"=if(nrow(tmpAln) > 1 && length(setdiff(1:nrow(tmpAln), tmp_susceptStrains)) > 1){
																		sapply(1:ncol(tmpAln),function(x){ 
																			tmpCols <- intersect(setdiff(1:nrow(tmpAln), tmp_susceptStrains),intersect(which(as.character(tmpAln[,x]) != "-"),which(tolower(as.character(tmpAln[,x])) != "n"))) 
																			return( nuc.div(tmpAln[tmpCols,x], pairwise.deletion = TRUE) )}) 
																	} else { 
																		99999 
																	}) )
	}
	save(nucDiv, nucDiv_intraGene, lociList, lociList_intraGene, allGenes, allLoci, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/nucDiv_lociList.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/nucDiv_lociList.RData")
}
# Now we will plot this nucDiv information to allow us to review if our sites of interest lay in regions of depressed nt diversity.
## SHOULD I BREAK THIS UP FOR EACH LOCI OF INTEREST TO SIMPLY LOOK AT INTS ENVIRONS???
pdf(file="~/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/piPlots.pdf",height=8,width=10)
nucDiv <- nucDiv[order(nucDiv$Loci),]
par(mfrow = c(2,1))
for (thisPi in c("pi_Suscept","pi_NotSuscept")){
	plot(x = NULL, y = NULL, xlim = c(min(nucDiv$refGenome_start),max(nucDiv$refGenome_end)), xlab = "Reference Genome Position",
			ylim = c(0,max(nucDiv[, c("pi_Suscept","pi_NotSuscept")])), ylab = "Nucleotide Diveristy (pi)",
			cex= 0.2, cex.lab = 1, cex.axis = 1, main = thisPi )
	# Now we add the Suscept and NotSuscept points individuall
	lines(x = unlist(sapply(1:nrow(nucDiv),function(x){ mean(nucDiv$refGenome_start[x],nucDiv$refGenome_end[x]) })), y = nucDiv[, thisPi], col= "black")
	
	# Now we create a bit of text and recolour the positions of our loci of interest
	for (thisLoci in allLoci){
		tmpIndex <- which(grepl(thisLoci,nucDiv$Loci))
		points(mean(nucDiv$refGenome_start[tmpIndex],nucDiv$refGenome_end[tmpIndex]),nucDiv[tmpIndex, thisPi], col = "red", pch= 20)
		# This will label the gene of interest which relates to this loci of interest, pos is to offset the text above
		text(mean(nucDiv$refGenome_start[tmpIndex],nucDiv$refGenome_end[tmpIndex]),nucDiv[tmpIndex, thisPi], labels= names(allLoci)[which(allLoci == thisLoci)], pos = 3, offset = 0.5, col = "red")
	}
}
dev.off()


nucDiv_intraGene <- nucDiv_intraGene[order(nucDiv_intraGene$Loci),]
# These are plotting colours we'll use downsteam
tmp_colSet <- list("pi_Suscept"="blue","pi_NotSuscept"="red")
for (thisLoci in allLoci){
	tmpRows <- sapply(lociList[[names(allLoci)[which(allLoci == thisLoci)] ]], function(x){ which(nucDiv_intraGene$Loci == x) }) 
	
	pdf(file=paste("~/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/",names(allLoci)[which(allLoci == thisLoci)],"pi_intraGene_Plots.pdf",sep=""),height=4,width=12)
	par(mfrow = c(1,2))
		
	for (thisPi in c("pi_Suscept","pi_NotSuscept")){
		# We define the part of this dataframe we are concerned with by finding those rows which relate to thisLoci and it's 
		# regionSize environs (which are stored in lociList_intraGene)
		plot(x = NULL, y = NULL, xlim = c(min(nucDiv_intraGene$refGenome_start[unlist(tmpRows)]),max(nucDiv_intraGene$refGenome_end[unlist(tmpRows)])), 
				xlab = "Reference Genome Position", ylim = c(0,max(nucDiv_intraGene[unlist(tmpRows), c("pi_Suscept","pi_NotSuscept")])), 
				ylab = "Nucleotide Diveristy (pi)",	cex= 0.2, cex.lab = 1, cex.axis = 1, main = paste(names(allLoci)[which(allLoci == thisLoci)], thisPi, sep=" :: ") )
	
		# I will now work out the x coordinates by taking the sum of the refGenome_start and the colnames stored nt_pos taken by the substr call
		tmp_xCord <- unlist(sapply(names(tmpRows), USE.NAMES=FALSE, function(x){ (nucDiv_intraGene$refGenome_start[ tmpRows[[x]][1] ] -1) + 
								as.numeric(substr(nucDiv_intraGene$colnames_cleanedAln[ tmpRows[[x]] ], regexpr("_",nucDiv_intraGene$colnames_cleanedAln[ tmpRows[[x]] ]) + 1, nchar(nucDiv_intraGene$colnames_cleanedAln[ tmpRows[[x]] ])))
							 }))
		# Now we add the Suscept and NotSuscept points individuall
		lines(tmp_xCord, nucDiv_intraGene[unlist(tmpRows), thisPi], col= tmp_colSet[[thisPi]], lwd = 2)
		
		# Now we create a bit of text and recolour the positions of our loci of interest
		for (thisMut in unlist(unique_sigPairs[,c("siteRef_1","siteRef_2")])[which(grepl(names(allLoci)[which(allLoci == thisLoci)], unlist(unique_sigPairs[,c("siteRef_1","siteRef_2")])))] ){
			tmpIndex <- which(nucDiv_intraGene$colnames_cleanedAln == thisMut)
			thisX <- (nucDiv_intraGene$refGenome_start[tmpIndex] - 1) + as.numeric(substr(thisMut, regexpr("_",thisMut) + 1, nchar(thisMut)))
			points(x= thisX , y= nucDiv_intraGene[tmpIndex, thisPi], col= "black", pch= 20, cex = 1.5)
			# This will label the gene of interest which relates to this loci of interest, pos is to offset the text above
			text(x= thisX, y= nucDiv_intraGene[tmpIndex, thisPi], labels= thisMut, pos = 3, offset = 1, col= "black")
		}	
	}
	dev.off()
}

load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/susceptStrains.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln_traitMat_5perc_wholeAlnTree.RData")
phenoVec <- sapply(rownames(polyAln),function(x){ 
	if(is.element(x, susceptStrains)){ 
		return( "Suscept" )  
	} else { 
		return( "Not_Suscept" ) } })
## This is a little review tool to be used when considering my nucDiv plots and finding that there is increased diversity at sites
## when the phenotype switches between Suscept and Not_Suscept.  Compare this information with the info in pairs_phenoMatch to help build a picture.
site_divTables <- list()
for (thisSite in unlist(unique_sigPairs[,c("siteRef_1","siteRef_2")]) ){
	#print( paste("Here is the table of states for ",thisSite,sep="") )
	tmpSite <- which(colnames(polyAln) == thisSite)
	tmpMat <- matrix(NA,nrow=2,ncol=2, dimnames=list(c("polyAln","traitMat"),c("Suscept", "Not_Suscept")))
	for (thisPheno in colnames(tmpMat)){
		#print( paste("For strains with this ",thisPheno," phenotype.",sep="") )
		for(thisData in rownames(tmpMat)){
			#print( paste("In ",thisData,sep="") )
			#print( table(eval(as.name(thisData))[unlist(sapply(1:length(phenoVec),function(x){ if(phenoVec[x] == thisPheno){ return(x) } })),tmpSite]) )
			tmpTable <- table(eval(as.name(thisData))[unlist(sapply(1:length(phenoVec),function(x){ if(phenoVec[x] == thisPheno){ return(x) } })),tmpSite])
			tmpMat[thisData,thisPheno] <-  paste(unlist(sapply(1:length(tmpTable), function(x){ paste(names(tmpTable)[x],tmpTable[x],sep=":") })), collapse=" -- ") 
		}
	}
	site_divTables[[thisSite]] <- tmpMat
}
save(site_divTables, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/site_divTables.RData")











# This will be used to plot the phenotype of strains which show at least one derived site in our sigPairs
rm(list=ls())
library(ape)
library(foreach)
library(doMC)
registerDoMC(cores=8)
options(digits = 6, stringsAsFactors=FALSE)

load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/susceptStrains.RData")
load('~/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/CompareSets_polyAln_0.8/q_1e-06/sigPairs_FALSE.RData')
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln_traitMat_5perc_wholeAlnTree.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")

phenoVec <- sapply(rownames(polyAln),function(x){ 
	if(is.element(x, susceptStrains)){ 
		return( "Suscept" )  
	} else { 
		return( "Not_Suscept" ) } })
		
if(!file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/pairs_phenoMatch.RData")){
	pairs_phenoMatch <- foreach(thisRow = 1:nrow(unique_sigPairs), .combine="rbind") %dopar% {
	#pairs_phenoMatch <- NULL
	#for(thisRow in 1:nrow(unique_sigPairs)) {
		
		siteNames <- unlist(unique_sigPairs[thisRow,c("siteRef_1","siteRef_2")])
		tmpSites <- unlist(sapply(siteNames,function(x){ which(colnames(polyAln) == x) }))
		# We will count the overlap between each independent mutation, and pairwise, with the "Not_Suscept" phenotype
		ind_NOTsusceptCorr <- sapply(tmpSites, function(x){
			# We count the proportion of times that have a non-zero state correlates with "Not_Suscept", then for "Suscept"
			return( length(intersect(which(traitMat[,x] != 0),which(phenoVec == "Not_Suscept")))/ length(which(traitMat[,x] != 0)) )
			})
		pair_NOTsusceptCorr <- sum(unlist(sapply(1:nrow(traitMat), function(x){
			if (is.element(1,traitMat[x,tmpSites])){
				if (phenoVec[x] == "Not_Suscept"){
					return( 1 ) 
				}
			} }) )) / sum(unlist(sapply(1:nrow(traitMat), function(x){ if (is.element(1,traitMat[x,tmpSites])){ return( 1 ) } })))
		# Now we want to perform some chi2 statistics on the correlation of phenotype and ind and paired states
		# this will involve us table-ing the fequency of states and then perform our chi2 analysis
		ind_stateArray <- array(0,dim=c(2,2,2),dimnames=list(c("Suscept","Not_Suscept"),c("0","1"),siteNames))
		pair_stateMat <- matrix(0,nrow=2,ncol=4,dimnames=list(c("Suscept","Not_Suscept"),c("0::0","1::0","0::1","1::1")))
		# Now we will go and create counts within our state tables
		for (thisRow in 1:nrow(traitMat)){
			# We deal with the independent states for each independent site, we use the dimnames to simply add counts in the array
			if (traitMat[thisRow,tmpSites[1]] != 2){
				ind_stateArray[phenoVec[thisRow],as.character(traitMat[thisRow,tmpSites[1]]), siteNames[1]] <- ind_stateArray[phenoVec[thisRow],as.character(traitMat[thisRow,tmpSites[1]]), siteNames[1]] + 1
			}
			if (traitMat[thisRow,tmpSites[2]] != 2){
				ind_stateArray[phenoVec[thisRow],as.character(traitMat[thisRow,tmpSites[2]]), siteNames[2]] <- ind_stateArray[phenoVec[thisRow],as.character(traitMat[thisRow,tmpSites[2]]), siteNames[2]] + 1
			}
			# Now we perform a similar act for the paired sites
			if (!(is.element(2,traitMat[thisRow,tmpSites]))){
				pair_stateMat[phenoVec[thisRow],paste(traitMat[thisRow,tmpSites],collapse="::")] <- pair_stateMat[phenoVec[thisRow],paste(traitMat[thisRow,tmpSites],collapse="::")] + 1
			}
		}
		
		# NOW WE CAN PERFORM A SUITE OF CHI2 TESTS FOR DIFFERENT QUESTIONS, I WANT TO PERFORM THE IND IN A STRAIGHTFORWARD MANNER,
		# BUT FOR THE PAIRED METHOD I SHOULD DISCOUNT THE 0::0 STATES SINCE I CANNOT ACCOUNT FOR HOW OFTEN A CHANGE WILL OCCUR BUT WHEN
		# A CHANGE DOES OCCUR IT MIGHT BE FAIR TO ASSUME THAT IF NOT CORREALTED TO PHENOTYPE THEN EQUAL NUMBER OF PHENO OCCURENCES
		# THESE CHI2 VALUES WILL BE USED TO SUGGEST IF THE paired (or sitewise)_NOTsusceptCorr VALUES ARE SIGNIFICANTLY DIFF FROM NULL
		ind_chi2 <- double(length=dim(ind_stateArray)[3])
		for (tmpSite in 1:dim(ind_stateArray)[3]){
			tmp_rowSum <- rowSums(ind_stateArray[,,tmpSite])
			tmp_colSum <- colSums(ind_stateArray[,,tmpSite])
			for (thisRow in 1:dim(ind_stateArray)[1]){
				for (thisCol in 1:dim(ind_stateArray)[2]){
					tmpExpect <- tmp_rowSum[thisRow] * tmp_colSum[thisCol] / sum(ind_stateArray[,,tmpSite])
					# We perform the chi2 calculation as ((O - E)^2 - 0.5 )/ E  - as per a multi dimensional chi2
					ind_chi2[tmpSite] <- ind_chi2[tmpSite] + if (tmpExpect == 0) { 0 } else {((abs(ind_stateArray[thisRow,thisCol,tmpSite] - tmpExpect) - 0.5)^2 /(tmpExpect))}
				}
			}
		}
		# This is the full pair chi2 value
		pair_chi2 <- double(length=1)
		tmp_rowSum <- rowSums(pair_stateMat)
		tmp_colSum <- colSums(pair_stateMat)
		for (thisRow in 1:dim(pair_stateMat)[1]){
			for (thisCol in 1:dim(pair_stateMat)[2]){
				tmpExpect <- tmp_rowSum[thisRow] * tmp_colSum[thisCol] / sum(pair_stateMat)
				# We perform the chi2 calculation as ((O - E)^2 - 0.5 )/ E  - as per a multi dimensional chi2
				pair_chi2 <- pair_chi2 + if (tmpExpect == 0) { 0 } else { ( (abs(pair_stateMat[thisRow,thisCol] - tmpExpect) - 0.5)^2 /(tmpExpect) ) }
			}
		}
		# Then I calculate a pair_chi2 value which excludes the 0::0 states, since I have no means to account for how often a change should happen
		# and my only intent is to measure that when a change does happen the changes are correlated to a particular phenotype, the null hypothesis being
		# that there is no correlation between state and phenotype
		pair_stateMat <- pair_stateMat[,-which(colnames(pair_stateMat) == "0::0")]
		pairAdj_chi2 <- double(length=1)
		tmp_rowSum <- rowSums(pair_stateMat)
		tmp_colSum <- colSums(pair_stateMat)
		for (thisRow in 1:dim(pair_stateMat)[1]){
			# Here we ignore the first column which is mean to represent the "0::0" states
			for (thisCol in 1:dim(pair_stateMat)[2]){
				tmpExpect <- tmp_rowSum[thisRow] * tmp_colSum[thisCol] / sum(pair_stateMat)
				# We perform the chi2 calculation as ((O - E)^2 - 0.5 )/ E  - as per a multi dimensional chi2
				pairAdj_chi2 <- pairAdj_chi2 + if (tmpExpect == 0) { 0 } else { ( (abs(pair_stateMat[thisRow,thisCol] - tmpExpect) - 0.5)^2 /(tmpExpect) ) }
			}
		}
		
		#pairs_phenoMatch <- rbind(pairs_phenoMatch , 
		return( data.frame("siteRef_1"=siteNames[1],"siteRef_2"=siteNames[2],
							"site1_NOTsusceptCorr"= ind_NOTsusceptCorr[which(names(tmpSites) == "siteRef_1")],
							"site1_chi2" = 1-pchisq(ind_chi2[1],1) ,
							"site2_NOTsusceptCorr" = ind_NOTsusceptCorr[which(names(tmpSites) == "siteRef_2")],
							"site2_chi2" = 1-pchisq(ind_chi2[2],1) ,
							"paired_NOTsusceptCorr"= pair_NOTsusceptCorr, 
							"pair_chi2" = 1-pchisq(pair_chi2,3),
							"pairAdj_chi2"= 1-pchisq(pairAdj_chi2,2)) )
	}
	save(pairs_phenoMatch,file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/pairs_phenoMatch.RData")
} else { 
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/q_e-06/pairs_phenoMatch.RData")
}


# Now we have some chi2 values for the sites identified being correalted to the phenotype, let's run the randomForest package and see if we can find
# where in the importance list these individual sites are found.  This will be a secondary data point to consider.
library(randomForest) 
library(tree)
# We reset the options for stringsAsFactors in order to satisfy using tree and randomForest (which use factors)
options(stringsAsFactors = TRUE)
load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")
# We establish the three phenotypes from the Kos (2015) data to help us calculate a better correalation between sites and phenotypes
kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
	if (is.element(x, c("PA14","PA01"))){ 
		return( "ND" )
	} else {
		as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
	} } ))


# We want to train our data before we cross validate it's decision tree building, so set the percentage of rows to use for training
# and then find the integer value of the number of rows to use.
trainPerc <- 0.4
trainSize <- round(nrow(polyAln[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)")),])* trainPerc,0)
chunkSize <- 1000
if (!(file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/dec_treeInfo.RData"))){
	# This will allow us to run through all of polyAln, find the sites most correlated with the phenotype assignment, then we can use this subset
	totalChunks <- ceiling(1:(ncol(polyAln)/chunkSize))
	dec_treeInfo <- list()
	polyAln_cols <- ncol(polyAln)
	for (thisChunk in totalChunks){
		tmpCols <- (1+((thisChunk - 1)* chunkSize)):(min(ncol(polyAln),thisChunk * chunkSize))
		tmpAln <- data.frame("Phenotype"=kos_phenoVec[which(kos_phenoVec != "ND")],polyAln[which(kos_phenoVec != "ND"), tmpCols], stringsAsFactors = TRUE)
		# We create a tree that could be used for training of a predictive tree but we do not inplement it directly here, but will store
		trainSet <- sample(1:nrow(tmpAln),trainSize,replace=FALSE)
		trainTree <- tree(Phenotype ~ .,data=tmpAln, split="deviance", subset=trainSet)
		# We create an intial decision tree from our current data set, this will allow us to perform cross.validation to prune for a better tree
		testTree <- tree(Phenotype ~ .,data=tmpAln, split="deviance")
		# Now we want to consider building a better tree that uses cross validation in order to establish the optimal set 
		# of nodse to be used.  Since we don't want the deviance to guide this criteria, instead we want the classification
		# error rate to guide our tree pruning, hence we use the FUN=prune.misclass
		crossvalTree <-cv.tree(testTree,FUN=prune.misclass)
		# While not necessarily the absolute best tree, we can prune our tree to be the one with the absolute lowest error ($dev)
		# NOTE: this may not be a significantly lower error rate than other possiblities, hence why error is plotted for review.
		#### However, since this is statistical machine learning in theory this smallest $dev value is the statistically best tree.
		pruneTree <- prune.misclass(testTree,best= crossvalTree$size[which(crossvalTree$dev == min(crossvalTree$dev))[1]])
		# We can look at the classification rate with our pruned tree by making predictions on our alignment, minus the training set
		predTree <- predict(pruneTree,tmpAln[-trainSet,],type="class")
	
		# What we want to keep and use at this stage is the sites that were identified as worth using in the decision tree, this is stored in the:
		# $frame$var list of the decision tree data.frames.  As at this stage we are being libral we will use the testTree informaiton, but keep the rest
		# as a measure of review and posterity.  We only need to keep any values that are not: <leaf> [these represent plotted leaves not nodes of decision]
		tmp_keepSites <- as.character(unique(testTree$frame$var))
		tmp_keepSites <- tmp_keepSites[which(tmp_keepSites != "<leaf>")]
		# NOTE: Tree will have changed any "GI:" into "GI.", but grepl will accept both of these special symbols and search them identically.
		
		dec_treeInfo[[thisChunk]] <- list("polyAln_col" = paste(tmpCols[1],tmpCols[length(tmpCols)],sep=":"), "trainTree"= trainTree, "testTree"= testTree, 
					"crossvalTree"= crossvalTree, "pruneTree" = pruneTree, "predTree" = predTree, "keepSites" = tmp_keepSites)
	}
	save(dec_treeInfo, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/dec_treeInfo.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/dec_treeInfo.RData")
}



# This reviews the types of mutations represented by any of the keepSites, this script is taken from above
if(!file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/all_keepSites_mutTypes.RData")){
	library(seqinr)
	library(ape)
	library(foreach)
	library(doMC)
	registerDoMC(cores = 4)
	
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/dec_treeInfo.RData")
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA14DBase/PA14_DataBase.RData")
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/wholeAln_index.RData")
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
	
	# Next to proceed we need to reload our polyAln and recover those portions that have been identified as keepSites
	all_keepSites <- unique(unlist(sapply(dec_treeInfo, function(x){ return( x$keepSites ) })))
	# We fix the GI. to be GI: values, this allows us to use this object for calling the assignment to tmpAln
	all_keepSites <- sub("GI.","GI:",all_keepSites)
	
	polyAln_sequences <- rownames(polyAln)
	rm(polyAln, dec_treeInfo, iBase, cBase)
	#################################################################################################
	##################################### BEGIN PRE FUNCTIONS #######################################
	#################################################################################################
	# This is the translate function from seqinr, but we needed to fix part of the code
	# So I've copied it here and will make a simple change
	translate <- function (seq, frame = 0, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE) {
	    #print("My version")
	    if (any(seq %in% LETTERS)) {
	        seq <- tolower(seq)
	    }
	    # I added a line to cahnge NA back into gap characters, this is because any NA passed to s2n causes a crashout
	    if (sens == "R") {
	        seq <- comp(rev(seq), ambiguous = ambiguous)
	        seq[is.na(seq)] <- "-"
	    } 
	    seqn <- s2n(seq, levels = s2c("tcag"))
	    l <- 3 * ((length(seq) - frame)%/%3)
	    c1 <- seq(from = frame + 1, to = frame + l, by = 3)
	    tra <- 16 * seqn[c1] + 4 * seqn[c1 + 1] + seqn[c1 + 2] + 
	        1
	    code <- s2c(SEQINR.UTIL$CODES.NCBI$CODES[numcode])
	    result <- code[tra]
	    result[is.na(result)] <- NAstring
	    if (ambiguous) {
	        toCheck <- which(result == NAstring)
	        for (i in toCheck) {
	            codon <- seq[c1[i]:(c1[i] + 2)]
	            allcodons <- as.vector(outer(as.vector(outer(amb(codon[1]), 
	                amb(codon[2]), paste, sep = "")), amb(codon[3]), 
	                paste, sep = ""))
	            allaminoacids <- sapply(allcodons, function(x) translate(s2c(x), 
	                numcode = numcode, ambiguous = FALSE))
	            if (all(allaminoacids == allaminoacids[1])) 
	                result[i] <- allaminoacids[1]
	        }
	    }
	    return(result)
	}
	#################################################################################################
	###################################### END PRE FUNCTIONS ########################################
	#################################################################################################
	
	# Here we load an alignment so that it can be translated, the alignment loaded is based on the Locus_Tag of the gene
	siteLocus <- unlist(sapply(all_keepSites,function(x){
		x <- strsplit(x,"_")[[1]][1]
		if (grepl("GI:",x)){
			return( as.character(dBase$Locus_Tag[which(grepl(x,dBase$Molecule_Xrefs))]) )
		} else {
			# There are instances of duplication among gene names, unfortunately.. so we keep the first one and handle possible errors
			# downstream, this is not ideal but moving forward without the original Locus_Tag I'd need to review all possible tmpAln's and split this more than I'd care to.
			return( as.character(dBase$Locus_Tag[which(dBase$Gene == x)[1]]) )
		} 
	}))
	# This is where we can load alignments from
	alnDir <- "/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/rawAln/"
	
	all_keepSites_mutTypes <- foreach(thisSite = 1:length(siteLocus), .combine="rbind") %dopar% {
		if ( is.na(siteLocus[thisSite]) ){ 
			return( data.frame("Site"=names(siteLocus[thisSite]),"mutType"="Error",stringsAsFactors=FALSE) )
		} else if ( grepl("Insert", names(siteLocus[thisSite])) ){
		 	return( data.frame("Site"=names(siteLocus[thisSite]),"mutType"="Insert",stringsAsFactors=FALSE) )
		} else {
			tmpAln <- read.dna(paste(alnDir,siteLocus[thisSite],"_seqAln.fas",sep=""),format="fasta",as.character=TRUE)
			# These colnames will help us keep track of the Ref and Insert _ntPos of sites as we clean.
			colnames(tmpAln) <- sapply(which(wholeAln_index$Locus_Tag == siteLocus[thisSite]), function(x){ 
				# we create the colnames using the Ref or Insert _ntPos information from our index file
				paste(wholeAln_index$Gene[x], if (is.na(wholeAln_index$Ref_ntPos[x])){ wholeAln_index$Insert_ntPos[x] } else { wholeAln_index$Ref_ntPos[x] }, sep="_") }) 
			
			# We make certain there are only the same strains as contained within polyAln
			tmpAln <- tmpAln[unlist(sapply(polyAln_sequences,function(x){ which(rownames(tmpAln) == x) })),]
	
			tmpTrans <- list()
			# We will identify any positions in this gene which our index defines as being insert positions
			tmp_indPos <- which(wholeAln_index$Locus_Tag == siteLocus[thisSite])
			tmpInserts <- unlist(sapply(tmp_indPos, function(x){ if(!is.na(wholeAln_index$Insert_ntPos[x])){ return( which(tmp_indPos == x) ) } }))
			for (tmpRow in 1:nrow(tmpAln)){	
				# We will remove any gap characters that belong to insert positions, these can be identified from the index
				tmpNucleotides <- tmpAln[tmpRow,]
				# If there are gaps at the insert positions in this sequence we will remove them before translating this information
				if (length(which(tmpNucleotides[tmpInserts] == "-")) > 0){ tmpNucleotides <- tmpNucleotides[-tmpInserts[which(tmpNucleotides[tmpInserts] == "-")]] }
				tmpTrans[[rownames(tmpAln)[tmpRow]]] <- translate(tmpNucleotides,sens= if(dBase$Strand[which(dBase$Locus_Tag == siteLocus[thisSite])] == 1){"F"} else {"R"}  )
			}
			# Now we look at the amino acid which relates to the ntPos, we need to know the strand to understand if we read forward or backward
			# for our divison of codons
			tmp_ntPos <- if (dBase$Strand[which(dBase$Locus_Tag == siteLocus[thisSite])] == 1){ 
				as.numeric(strsplit(all_keepSites[thisSite],"_")[[1]][2])
			} else {
				ncol(tmpAln) + 1 - as.numeric(strsplit(all_keepSites[thisSite],"_")[[1]][2])
			}
			
			if (tmp_ntPos > ncol(tmpAln)){
				# This is an error of our above random assignment of siteLocus information, this ought to affect very few instances
				return( "Error" )
				#responseVec <- c(responseVec, "Error")
			} else {
				tmp_aaPos <- ceiling( tmp_ntPos / 3 )
				# If there is more than one character, other than "X" (which is ambiguity) than this was a non-synonymous change
				tmpTable <- table(sapply(tmpTrans,function(x){ x[tmp_aaPos] }))
				if (length(tmpTable[which(names(tmpTable) != "X")]) == 1){
					return( data.frame("Site"=names(siteLocus[thisSite]),"mutType"="Synonymous",stringsAsFactors=FALSE) )
				} else if (length(tmpTable[which(names(tmpTable) != "X")]) > 1) {
					return( data.frame("Site"=names(siteLocus[thisSite]),"mutType"="Non_Synonymous",stringsAsFactors=FALSE) )
				}
			}
		}
	}
	# We save this since it is a costly computation
	save(all_keepSites_mutTypes, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/all_keepSites_mutTypes.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/all_keepSites_mutTypes.RData")
}	



### This is a small test to find out if, when using my large alignments, it's better that I generate larger forests in order to reduce 
### cross analysis variation in the importance of a value
if(!file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/forestSize_trials.RData")){
	library(foreach)
	library(doMC)
	library(randomForest)
	registerDoMC(cores = 4)
	
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
	load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")
	
	kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
	if (is.element(x, c("PA14","PA01"))){ 
		return( "ND" )
	} else {
		as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
	} } ))
	# Tidy some room
	rm( bacFrame )
	# We will build a trial alignment, and then perform several randomForest growths and see the variation in importance values defined.
	# We use a similar approach to what will be used in other sets, where we remove the classes that have imbalance (the ND and (I) types)
	trialAln <- data.frame("Phenotype"=kos_phenoVec[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))],
						polyAln[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)")), sample(1:ncol(polyAln),2000, replace=FALSE)], stringsAsFactors = TRUE)
	
	trainPerc <- 0.4
	trainSize <- round(nrow(polyAln[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)")),])* trainPerc,0)
	# again let's tidy some space
	rm(dec_treeInfo, polyAln)
	# We create a training set to be passed to tree and randomForest calls so that we can evaluate the performance between sets
	trainSet <- sample(1:nrow(trialAln),trainSize,replace=FALSE)
	phenoTest <- trialAln[-trainSet,"Phenotype"]
	
	startTime <- proc.time()
	# Now we build some random forest objects, these are replicate forests, with the default size of 500 trees
	defaultForests <- foreach(thisForest = 1:8) %dopar%{
		trainForest <- randomForest(Phenotype ~ ., data= trialAln, subset = trainSet, mtry = floor(sqrt(ncol(trialAln[,-1]))), importance = TRUE)
		predForest <- predict(trainForest, newdata= trialAln[-trainSet,])
		importForest <- randomForest(Phenotype ~ ., data= trialAln, mtry = floor(sqrt(ncol(trialAln[,-1]))), importance = TRUE)
		
		return( list("trainForest"= trainForest,"predForest"= predForest,"importForest"= importForest) )
	}
	defaultTime <-  proc.time() - startTime
	
	
	startTime <- proc.time()
	# Now we build some random forest objects, these are replicate forests, with a number of trees which varies based on the number of predictors.
	# To ensure better probability that a predictor has been tried among all other's, we will create a forest of twice the size of our predictors, this means
	# that each variable will have been tried, on average, the mtry number x2 of times meaning it's been sampled with up to 2*(mtry * (mtry -1)) other predictors
	# which should mean each other predictor more than one but less than twice.
	variableForests <- foreach(thisForest = 1:8) %dopar%{
		trainForest <- randomForest(Phenotype ~ ., data= trialAln, subset = trainSet, ntree=(2*ncol(trialAln[,-1])), mtry = floor(sqrt(ncol(trialAln[,-1]))), importance = TRUE)
		predForest <- predict(trainForest, newdata= trialAln[-trainSet,])
		importForest <- randomForest(Phenotype ~ ., data= trialAln, ntree=(2*ncol(trialAln[,-1])), mtry = floor(sqrt(ncol(trialAln[,-1]))), importance = TRUE)
		
		return( list("trainForest"= trainForest,"predForest"= predForest,"importForest"= importForest) )
	}
	variableTime <- proc.time() - startTime
	
	startTime <- proc.time()
	# Now we build some random forest objects, these are replicate forests, with a number of trees which varies based on the number of predictors.
	# To ensure better probability that a predictor has been tried among all other's, we will create a forest of twice the size of our predictors, this means
	# that each variable will have been tried, on average, the mtry number x10 of times meaning it's been sampled with up to 10*(mtry * (mtry -1)) other predictors
	# which should mean each other predictor more than one but less than twice.
	variableForests_large <- foreach(thisForest = 1:8) %dopar%{
		trainForest <- randomForest(Phenotype ~ ., data= trialAln, subset = trainSet, ntree=(10*ncol(trialAln[,-1])), mtry = floor(sqrt(ncol(trialAln[,-1]))), importance = TRUE)
		predForest <- predict(trainForest, newdata= trialAln[-trainSet,])
		importForest <- randomForest(Phenotype ~ ., data= trialAln, ntree=(10*ncol(trialAln[,-1])), mtry = floor(sqrt(ncol(trialAln[,-1]))), importance = TRUE)
		
		return( list("trainForest"= trainForest,"predForest"= predForest,"importForest"= importForest) )
	}
	variable_largeTime <- proc.time() - startTime
	
	
	startTime <- proc.time()
	# Now we build some bagging forest objects, these are replicate forests to see if bagging will be the best way to reduce our variance in importance values
	baggingForests <- foreach(thisForest = 1:8) %dopar%{
		trainForest <- randomForest(Phenotype ~ ., data= trialAln, subset = trainSet, mtry = ncol(trialAln[,-1]), importance = TRUE)
		predForest <- predict(trainForest, newdata= trialAln[-trainSet,])
		importForest <- randomForest(Phenotype ~ ., data= trialAln, mtry = ncol(trialAln[,-1]), importance = TRUE)
		
		return( list("trainForest"= trainForest,"predForest"= predForest,"importForest"= importForest) )
	}
	baggingTime <- proc.time() - startTime
	
	### HERE ###
	# Now we will extract the importance matrices for each of these and work out the variance in  MeanDecreaseGini for our sites
	meanMeasures <- sdMeasures <- matrix(NA,ncol=4,nrow=ncol(trialAln[,-1]),dimnames=list(colnames(trialAln[,-1]),c("defaultForests","variableForests","variableForests_large","baggingForests")))
	# for each of the types will go and calculate the sd for each 
	for (thisType in colnames(meanMeasures) ){
		# We will go through each of the sites and then perform sd on the MeanDecreaseGini for that site across all replicates (there are 8)
		tmpImportance <- foreach(thisSite = rownames(meanMeasures), .combine="rbind") %dopar% {
			tmpVal <- unlist(sapply(1:8,function(x){ 
				eval(as.name(thisType))[[x]]$importForest$importance[which(rownames(eval(as.name(thisType))[[x]]$importForest$importance) == thisSite),
																		"MeanDecreaseGini"] }))
			return( data.frame("Site"=thisSite, "Mean"=mean(tmpVal), "STDEV"=sd(tmpVal), stringsAsFactors=FALSE) )
		}
		# We re-order our tmpImportance based on how the Sites match our meanMeasures rownames
		corrOrder <- unlist(sapply(rownames(meanMeasures),function(x){ which(tmpImportance$Site == x) }))
		# Now we assign the information to our meanMeasures and sdMeasures matrices
		meanMeasures[,thisType] <- tmpImportance$Mean[corrOrder]
		sdMeasures[,thisType] <- tmpImportance$STDEV[corrOrder]
	}
	
	#This will compare the predictive quality of each forest type and return these values' mean and sd
	# NOTE: this has been coded to work for classification trees with 2 types to be classified.
	predictiveMeasures <- matrix(NA,ncol=4,nrow=4,dimnames=list(c("meanCorrect","sdCorrect","meanIncorrect","sdIncorrect"),c("defaultForests","variableForests","variableForests_large","baggingForests")))
	# for each of the types will go and calculate the sd for each 
	for (thisType in colnames(meanMeasures) ){
		# We will go through each of the predForest confusion table values and return the mean and sd values
		tmpVal <- unlist(sapply(1:8,function(x){ table(eval(as.name(thisType))[[x]]$predForest,trialAln$Phenotype[-trainSet]) }))
		predictiveMeasures["meanCorrect",thisType] <- mean(tmpVal[c(1,4),])/ length(trialAln$Phenotype[-trainSet])
		predictiveMeasures["sdCorrect",thisType] <- sd(tmpVal[c(1,4),])/ length(trialAln$Phenotype[-trainSet])
		predictiveMeasures["meanIncorrect",thisType] <- mean(tmpVal[c(2,3),])/ length(trialAln$Phenotype[-trainSet])
		predictiveMeasures["sdIncorrect",thisType] <- sd(tmpVal[c(2,3),])/ length(trialAln$Phenotype[-trainSet])
	}
	predictiveMeasures <- round(predictiveMeasures,4)
	# We save this since it is a costly computation
	save(trialAln, trainSet, defaultForests, defaultTime, variableForests, variableTime, variableForests_large, variable_largeTime, baggingForests, baggingTime, 
			meanMeasures, sdMeasures, predictiveMeasures, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/forestSize_trials.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/forestSize_trials.RData")
}
# Now lets visualise these values a bit and get a sense if the distribution of standard deviations changes across types
# But first we will normalise this standard deviation as a function of the mean it surrounds, this is because when the mean is small
# a small standard deviation "means" more.
tmp_plotCols <- list("defaultForests"="black","variableForests"="red","variableForests_large"="blue","baggingForests"="purple")
# This ought not to be strictly necessary but ;et's be explicit and calculate the correct order for comparing sd and mean Measures
tmp_corrOrder <- unlist(sapply(rownames(sdMeasures),function(x){ which(rownames(meanMeasures) == x) }))
# Here we calculate the normalised standard deviation
normalised_sdMeasures <- sdMeasures/meanMeasures[tmp_corrOrder,]
# Now we change any NaN values (caused by division by zero) to be zero values assuming (possibly incorrectl but unlikely) that sd was zero as well
normalised_sdMeasures[which(normalised_sdMeasures == "NaN")] <- 0
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/forestSize_trials_sdMeasures.pdf", height= 4 * ncol(sdMeasures), width=12)
par(mfrow=c(ncol(sdMeasures),2))
for(thisForest in colnames(sdMeasures)){
	tmpHist <- hist(normalised_sdMeasures[, thisForest], plot = FALSE)
	tmpHist$counts <- log(tmpHist$counts)
	# We make the counts into log format to make it easier to see when there are small values
	plot(tmpHist, ylim = log(c(1,1000)), xlim = c(0,max(normalised_sdMeasures)*1.1), col = tmp_plotCols[[thisForest]], 
			xlab = paste("Standard Deviation",sep=""), ylab="Frequency (ln)", main = NULL ) #main = paste("Distribution of stand. dev. for ",thisForest,sep="") )
	# Now we plot some text in the plot, this will be the predictive quality of the types and the run time
	legend("topright", cex=1,"", col= "black", bty ="n", horiz=FALSE, pch = 16,
			legend=c(paste("Correct Classification: ", predictiveMeasures[1, thisForest], ' +/- ', predictiveMeasures[2, thisForest],  sep=""),
					paste("In-Correct Classification: ", predictiveMeasures[3, thisForest], ' +/- ', predictiveMeasures[2, thisForest],  sep=""),
					paste("Run Time: ", eval(as.name(paste(sub("Forests","", thisForest,),"Time",sep="")))[3], sep="")) )
	
	# Now we plot the importForest object to show if the OOB error had calmed and where is was found to be, we simply plot and over-plot all 
	# 8 instances for thisForest
	plot( eval(as.name(thisForest))[[1]]$importForest, ylim = c(0.10,0.40), main=NULL )
	for (x in 2:8){
		plot( eval(as.name(thisForest))[[x]]$trainForest, add=TRUE, ylim = c(0.10,0.40), main=NULL )
	}
}
legend("top",cex=1 ,colnames(sdMeasures),bty="n",fill= unlist(tmp_plotCols), horiz=TRUE)
dev.off()





# This builds our decision trees and random forests using the sites "worth keeping" as per our initial round of CART [tree() calls]
if(!file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/tree_randomForest_objects.RData")){
	library(tree)
	library(randomForest)
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/dec_treeInfo.RData")
	load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")
	# We establish the three phenotypes from the Kos (2015) data to help us calculate a better correalation between sites and phenotypes
	kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
		if (is.element(x, c("PA14","PA01"))){ 
			return( "ND" )
		} else {
			as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
		} } ))
	# We want to train our data before we cross validate it's decision tree building, so set the percentage of rows to use for training
	# and then find the integer value of the number of rows to use.
	trainPerc <- 0.4
	trainSize <- round(nrow(polyAln[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)")),])* trainPerc,0)

	# Next to proceed we need to reload our polyAln and recover those portions that have been identified as keepSites
	all_keepSites <- unique(unlist(sapply(dec_treeInfo, function(x){ return( x$keepSites ) })))
	# We fix the GI. to be GI: values, this allows us to use this object for calling the assignment to tmpAln
	all_keepSites <- sub("GI.","GI:",all_keepSites)
	
	## From Stephane, and I agree, I will remove the strains which are of (I) - intermediate class, the reason being that there is a class imbalance
	## They represent only ~ 10% of the total alignment, and this can cause False Positives since the correaltion coefficient of the regressions can be 
	## affected by this.  NOTE:  We do not bother in the previous, as FP mean more sites to be refined, it does not mean we might lose sites of true interest.
	tmpAln <- data.frame("Phenotype"=kos_phenoVec[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))],
						polyAln[intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)")), unlist(sapply(all_keepSites,function(x){
							which(colnames(polyAln) == x)}))], stringsAsFactors = TRUE)
	# again let's tidy some space
	rm(dec_treeInfo, polyAln)
	# We create a training set to be passed to tree and randomForest calls so that we can evaluate the performance between sets
	trainSet <- sample(1:nrow(tmpAln),trainSize,replace=FALSE)
	phenoTest <- tmpAln[-trainSet,"Phenotype"]
	
	trainTree <- tree(Phenotype ~ .,data=tmpAln, split="deviance", subset=trainSet)
	crossvalTree <-cv.tree(trainTree,FUN=prune.misclass)
	# This prunes our tree and takes the first smallest $size argument, meaning we are accepting the largest tree which delivers the min crossval deviation
	pruneTree_large <- prune.misclass(trainTree, best= crossvalTree$size[which(crossvalTree$dev == min(crossvalTree$dev))[1]])
	predTree_large <- predict(pruneTree_large,tmpAln[-trainSet,],type="class")
	# Here we prune but take the smallest of the $size values, we don't verify that there was more than one answer, but in plotting this will become obvious
	pruneTree_small <- prune.misclass(trainTree, best= crossvalTree$size[which(crossvalTree$dev == min(crossvalTree$dev))[length(which(crossvalTree$dev == min(crossvalTree$dev)))]] )			
	predTree_small <- predict(pruneTree_small,tmpAln[-trainSet,],type="class")
	
	# Now we build some random forest objects, and explicitly set the number of predictors attempted in each tree growth based on the sqrt()
	# this is from the Statistical machine learning book which details these functions and the mathematics behind them, it is chosen since
	# when there are a large number of predictors (such as sequence information), using a smaller number of predictors at each growth reduces MSE
	trainForest <- randomForest(Phenotype ~ ., data=tmpAln, ntree=(10 * ncol(tmpAln[,-1])), subset = trainSet, mtry = floor(sqrt(ncol(tmpAln[,-1]))), importance = TRUE)
	predForest <- predict(trainForest, newdata= tmpAln[-trainSet,])
	# We use the information from above to generate better trees by increasing the size of our forests
	importForest <- randomForest(Phenotype ~ ., data=tmpAln, ntree=(10 * ncol(tmpAln[,-1])), mtry = floor(sqrt(ncol(tmpAln[,-1]))), importance = TRUE)
	
	# We save this since it is a costly computation
	save(tmpAln, trainSet, phenoTest, trainTree, crossvalTree, pruneTree_large, predTree_large, pruneTree_small, predTree_small, trainForest, predForest, importForest,file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/tree_randomForest_objects.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/tree_randomForest_objects.RData")
}

###### PLOTTING PLAY AREA TO REVIEW TREE AND RANDOMFOREST WORK ######
# These are some review calls which can be made on tree objects to review the decision trees made
plot(trainTree)
text(trainTree,pretty = 0)
# This will show us the cross validation information telling us an optimal tree size
plot(crossvalTree)
# This will allow us to review how our predictions on our non-training set compare to the phenotypes of that set, the table is effectively a confusion table
plot(predTree_large, phenoTest, xlab = "Predicted Phenotype", ylab = "Actual Phenotype")
table(predTree,phenoTest)

# These are some review calls which are practical for the randomForest objects, it shows the error rate of the forest 
plot(importForest)
varImpPlot(importForest)
# Review how good our predictions were from our trained randomForest 
plot(predForest, phenoTest)
table(predForest, phenoTest) # This has turned out worse than our predTree


# This will try and create a full list of the importance for all sites within polyAln, it will use the first four prioritised sites as items to normalise importance across runs
prioritySites <- rownames(importForest$importance)[order(importForest$importance[,"MeanDecreaseGini"], decreasing = TRUE)[1:4]]
# If there are any "GI:" type sites they will have been restrung as "GI.", we fix this
prioritySites <- sub("GI.","GI:", prioritySites)
# We first tidy a bit of space
rm(tmpAln, trainSet, phenoTest, trainTree, crossvalTree, pruneTree_large, predTree_large, pruneTree_small, predTree_small, trainForest, predForest, importForest, dec_treeInfo)
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
if (!(file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest.RData"))){
	# This will allow us to run through all of polyAln, find the sites most correlated with the phenotype assignment, then we can use this subset
	totalChunks <- ceiling(1:(ncol(polyAln)/chunkSize))
	fullForest <- list()
	polyAln_cols <- ncol(polyAln)
	# This will inform which elements of our data set are properly represented within the data and worth keeping, ND and I states are imbalanced classes
	keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
	for (thisChunk in totalChunks){
		tmpCols <- (1+((thisChunk - 1)* chunkSize)):(min(ncol(polyAln),thisChunk * chunkSize))
		tmpAln <- data.frame("Phenotype"=kos_phenoVec[keepRows],polyAln[keepRows, tmpCols], stringsAsFactors = TRUE)
		# Now we review iftmpAln includes all of our prioritySites, if it is missing any we add them
		addCols <- which(unlist(sapply(prioritySites, function(x){ !is.element(x,colnames(tmpAln)) })))
		if (length(addCols) > 0){ tmpAln <- data.frame(tmpAln, polyAln[keepRows, prioritySites[addCols]], stringsAsFactors = TRUE) }
		# We create a forest and ensure that our priority sites are included (this will normalise our importance values)
		fullForest[[thisChunk]] <- randomForest(Phenotype ~ ., data=tmpAln, ntree=(10 * ncol(tmpAln[,-1])), mtry = floor(sqrt(ncol(tmpAln[,-1]))), importance = TRUE)
	}
	save(fullForest, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest.RData")
}


# This builds a full forest without any normalising sites included, this is to be compared against that with the normalisation.	
if (!(file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNormal.RData"))){
	library(randomForest)
	library(foreach)
	library(doMC)
	registerDoMC(cores=4)
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
	load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")
	# We establish the three phenotypes from the Kos (2015) data to help us calculate a better correalation between sites and phenotypes
	kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
		if (is.element(x, c("PA14","PA01"))){ 
			return( "ND" )
		} else {
			as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
		} } ))
	rm( bacFrame) # tidy up some room
	# This will allow us to run through all of polyAln, find the sites most correlated with the phenotype assignment, then we can use this subset
	chunkSize <- 1000
	totalChunks <- 1:ceiling(ncol(polyAln)/chunkSize)
	fullForest_noNormal <- vector("list",length= max(totalChunks))
	polyAln_cols <- ncol(polyAln)
	# This will inform which elements of our data set are properly represented within the data and worth keeping, ND and I states are imbalanced classes
	keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
	## NOTE, from experience we can't use foreach, the file sizes are too large for transfer to happen properly.
	for (thisChunk in totalChunks){
	#fullForest <- foreach(thisChunk = totalChunks) %dopar%{	
		tmpCols <- (1+((thisChunk - 1)* chunkSize)):(min(ncol(polyAln),thisChunk * chunkSize))
		tmpAln <- data.frame("Phenotype"=kos_phenoVec[keepRows],polyAln[keepRows, tmpCols], stringsAsFactors = TRUE)
		# We create a forest and ensure that our priority sites are included (this will normalise our importance values)
		fullForest_noNormal[[thisChunk]] <- randomForest(Phenotype ~ ., data=tmpAln, ntree=(10 * ncol(tmpAln[,-1])), mtry = floor(sqrt(ncol(tmpAln[,-1]))), importance = TRUE)
		#return( randomForest(Phenotype ~ ., data=tmpAln, ntree=(10 * ncol(tmpAln[,-1])), mtry = floor(sqrt(ncol(tmpAln[,-1]))), importance = TRUE) )
	}
	save(fullForest_noNormal, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNormal.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNormal.RData")
}
# Now I will extract just the importance values from the fullForest_noNormal object
if (!(file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.RData"))){
	library(randomForest)
	library(foreach)
	library(doMC)
	registerDoMC(cores=4)
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNormal.RData")

	fullForest_noNrml_importance <- foreach(thisForest = fullForest_noNormal, .combine="rbind") %dopar% {
		return( thisForest$importance )
	}
	fullForest_noNrml_importance <- fullForest_noNrml_importance[order(fullForest_noNrml_importance[,"MeanDecreaseGini"], decreasing = TRUE),]

	save(fullForest_noNrml_importance, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.RData")
}
# We load this data to make the special coloured plot for phenotype and sequence alignments
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")	
kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
if (is.element(x, c("PA14","PA01"))){ 
	return( "ND" )
} else {
	as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
} } ))
keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
# these are plotting colours for our grid plot
ntCols <- list("A"="orange","C"="green","G"="purple","T"="yellow","-"="black")
phenoCols <- list("(S)"="blue","(R)"="red")
	
# Here we will plot the importance values to get a sense of the data set
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.pdf",height=12,width=6)
par(mfcol=c(2,1))
tmpHist <- hist(fullForest_noNrml_importance[,"MeanDecreaseGini"], plot = FALSE)
# There are too many sites at which the importance is effectively zero, we will rescale the counts
tmpHist$counts <- log10(tmpHist$counts)
# We need to re-assign the zero values as being -1, since counts can't be less than 1 before reaching zero this represents a zero value
tmpHist$counts[which(tmpHist$counts == -Inf)] <- -1
plot(tmpHist, xlab = "Mean Decrease Gini", ylab = "Frequency (log10)", main = NULL, ylim = c(-1,signif(max(tmpHist$counts),3)), xlim = c(0,max(tmpHist$breaks)*1.1) )
# I don't want to plot those values which are in less than the first 10,000 positions
plotThresh <- 10000
sigThresh <- 10
plot(y = fullForest_noNrml_importance[1: plotThresh,"MeanDecreaseGini"],x = log(1: plotThresh), xlab = "Ordered Importance (log)", ylab = "Mean Decrease Gini",
		type = "b", lwd = 2, pch = 15, bty = "n")
# Now we add some text to the plot for the first X most important sites
text(y = fullForest_noNrml_importance[1: sigThresh,"MeanDecreaseGini"],x = log(1: sigThresh), labels=rownames(fullForest_noNrml_importance)[1: sigThresh], pos = 3, offset = 1)
dev.off()

# Now we'll use the top X most interesting sites and look at a "plot" of their phenotype's by nucleotide states nicely coloured
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_Mostimportant.pdf",height=36,width=4)
plotMat <- cbind(kos_phenoVec[keepRows],toupper(polyAln[keepRows,sub("GI.","GI:",rownames(fullForest_noNrml_importance)[1: sigThresh])]))
# Now we create a blank plot space and then will simply plot coloured pch squares based on the phenotype and nt types
plot(x = NULL, y = NULL, xlim = c(1, sigThresh + 1), ylim = c(1, nrow(plotMat)), ylab = "Phenotype and Nucleotide Alignment", xlab = "", axes = FALSE)
for(thisType in c("phenoCols","ntCols")){
	for(thisState in names(eval(as.name(thisType))) ){
		tmpPlot <- which(plotMat == thisState,arr.ind = TRUE)
		points(x = tmpPlot[,"col"], y = tmpPlot[,"row"], pch = 15, cex = 1, col = eval(as.name(thisType))[[thisState]])
	}
}
# We will add text for the 
text(x = 0:sigThresh + 1, y = 1, labels = c("Phenotype",rownames(fullForest_noNrml_importance)[1: sigThresh]) ,pos = 1, offset = 2, srt = 90)
legend("bottom",cex=1 ,c(names(phenoCols),names(ntCols)),bty="n",fill= unlist(list(phenoCols, ntCols)), ncol = 4)
dev.off()



### HERE TO GET NEW NORMALISED IMPORTANCE VALUES ###
# Now we normalise all of our importance information, we will normalise to a value of 1 for all values in the importance matrix 
# however the normalisation has 4 options, each of our priorty sites, so we need to create a normalised value set for each of these
if (!(file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_rankChange.RData"))){
	# We store this separately as the fullForest is a large object for polyAln (120 MB)
	# Step 1: tidy some room
	rm(polyAln)
	registerDoMC(cores = 4)
	normal_fullForests <- list()
	for (thisSite in prioritySites){
		# we adjust thisSite for our internal calls since "GI:" are still stored as "GI." values
		tmpSite <- sub("GI:","GI.",thisSite)
		normal_fullForests[[thisSite]] <- foreach(thisForest = fullForest, .combine="rbind") %dopar% {
			# We will normalise based on the MeanDecreaseAccuracy since it's values are already normalised between 0,1 and the range might be lower
			tmp_normalisingValues <- thisForest$importance[tmpSite,]
			# Now we simply divide each row of the matrix by this vector of normalisation
			for(tmpRow in 1:nrow(thisForest$importance)){
				thisForest$importance[tmpRow,] <- thisForest$importance[tmpRow,] / tmp_normalisingValues
			}
			return( thisForest$importance )
		}
		
		# We remove our normalising elements by keeping only a single instance of our normalising value, and the mean values for our other sites
		for (cleanSite in prioritySites){
			tmp_cleanSite <- sub("GI:","GI.", cleanSite)
			# We find all the rows which are for this cleanSite 
			tmpRows <- which(rownames(normal_fullForests[[thisSite]]) == tmp_cleanSite)
			# We apply our mean function to the each column for all tmpRows
			normal_fullForests[[thisSite]][tmpRows[1], ] <- unlist(apply(normal_fullForests[[thisSite]][tmpRows, ], 2, mean))
			if (length(tmpRows) > 1){ normal_fullForests[[thisSite]] <- normal_fullForests[[thisSite]][-tmpRows[-1],] }
		}
		# Then we order each of these based on decreasing MeanDecreaseGini values
		normal_fullForests[[thisSite]] <- normal_fullForests[[thisSite]][order(normal_fullForests[[thisSite]][,"MeanDecreaseGini"],decreasing=TRUE),]
	}
	
	# Let's evaluate the change in rank value based on how we've normalised the data, this will give us a sense if there was any value in doing so
	# Now this requires that the normal_fullForests all be of the same length, which they ought to be
	fullForest_rankChange <- NULL
	identical_uniqueRows <- sapply(normal_fullForests,function(x){ unique(rownames(x)) })
	# If this has come out as a matrix then the size of the unique rowsnames in all sets is equal, this is a good start
	if (is.matrix(identical_uniqueRows)){
		# Now we review if there is any setdiff's between our identical unique rownames 
		if( all(sapply(1:(ncol(identical_uniqueRows) - 1), function(x){ length(setdiff(identical_uniqueRows[,x],identical_uniqueRows[,x+1]) )== 0 } )) ){
			# Then this means we've been able to extract the same indentical unique rows from each set and can work out the rank changes
			# We use one of the columns and grab the values stores, these are the names to be used
			allNames <- identical_uniqueRows[,1]
			# Tidy some space
			rm( identical_uniqueRows ) 
			# Now we go through each rowname and find the rank (or mean rank as SOME very few instance have multiple entries - these are poor quality and ignorable)
			fullForest_rankChange <- foreach(thisName = allNames, .combine="rbind") %dopar% {
				# This call exists to speed up how quickly we can get this matrix built, the challenge being that it is rather large ( ~492k rows)
				return( matrix( unlist(sapply(names(normal_fullForests),function(x){
									round(mean(which(rownames(normal_fullForests[[x]]) == thisName)),0) })), nrow = 1, 
								ncol = length(names(normal_fullForests)), dimnames=list(thisName,paste(names(normal_fullForests),"normalised",sep="_"))) )
			}
			
			# Now we calculate the mean change in rank for each of these sites, this is done by performing all pairwise rank comparisons and dividing by
			# the number of comparisons
			rankChange <- rep(0,nrow(fullForest_rankChange))
			divCounter <- 0
			for (i in 1:(ncol(fullForest_rankChange)-1) ){
				for(j in (i+1):ncol(fullForest_rankChange)){
					rankChange <- rankChange + abs(fullForest_rankChange[,i] - fullForest_rankChange[,j])
					divCounter <- divCounter + 1
				}
			}
			rankChange <- rankChange/ divCounter
			fullForest_rankChange <- cbind(fullForest_rankChange, round(rankChange,3) )
			colnames(fullForest_rankChange) <- c(paste(names(normal_fullForests),"normalised",sep="_"), "rankChange")
			
			# Now we will extract the stdev and mean values for each site, as generated by the normalised sets
			fullForest_Measures <- foreach(thisName = allNames, .combine="rbind") %dopar% {
				# This is the importance values as calculated by the MeanDecreaseGini for thisSite (reported as thisName) 
				# cross of 4 normalised forests.  A mean call is made here as there are some VERY RARE instances when a rownames is not unique
				# this only occurs from GI: type or other error named sites (They have V#### or Entrez##### types... safe to ignore)
				tmpMeasures <- unlist(sapply(names(normal_fullForests),function(x){
									mean(normal_fullForests[[x]][which(rownames(normal_fullForests[[x]]) == thisName),"MeanDecreaseGini"]) }))
				
				return( data.frame("Site"= thisName,"Mean"=round(mean(tmpMeasures),3),"STDEV"=round(sd(tmpMeasures),3), stringsAsFactors=FALSE) )
			}
		} else {
			fullForest_rankChange <- fullForest_Measures <- "Unique rownames between sets are not all identical, though the number of unique rownames is, review"
		}
	} else {
		fullForest_rankChange <- fullForest_Measures <- "There are a different number of unique rownames between sets, review"
	}
	
	# Now we save our information for review, though this step is relatively quick.	
	save(normal_fullForests, fullForest_rankChange, fullForest_Measures, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_rankChange.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_rankChange.RData")
}
library(locfit) # This is to perform the locfit density plots, best means of providing a kernl smoothed density plot.
# This will help visualise the amount of importance rank order changes which have occured across the 4 prioirty sites used to normalise our data
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_rankChange.pdf", height= 8, width=6)
par(mfrow=c(2,1))
# First we will plot the historgram of rankChange caused across all 4 sets 
tmpHist <- hist(fullForest_rankChange[, "rankChange"], plot = FALSE)
plot(tmpHist, ylim = c(1,signif(max(tmpHist$counts),3)), xlim = c(0,max(tmpHist$breaks)*1.1), 
			xlab = paste("Normalised Ranked Order Change ( [SUM_i in 1:n-1 | SUM_j in i+1:n | abs(rank_i - rank_j)] / n(n-1)/2 )",sep=""), 
			ylab="Frequency Counts", main = NULL ) #main = paste("Distribution of stand. dev. for ",thisForest,sep="") )
# Then we plot the amount of rank change as a function of the mean normalised importance value
# The expectation here being that the rank changes are largely a function of un-important sites being arbitrary and changing small amounts
	#plot(x = fullForest_Measures$STDEV/fullForest_Measures$Mean, y = fullForest_rankChange[, "rankChange"], 
	#		xlab = "Normalised Std. Deviation of Importance (Mean Decrease Gini)", 
	#		ylab = "Normalised Ranked Order Change ( [SUM_i in 1:n-1 | SUM_j in i+1:n | abs(rank_i - rank_j)] / n(n-1)/2 )" )
# We plot something similar to the above but only the density of this scatter plot, this is largely because the number of points in the 
# above makes for a cumbersome and unweildly plot.
plotFrame <- data.frame("x"=fullForest_Measures$STDEV/fullForest_Measures$Mean, "y" = fullForest_rankChange[corrOrder, "rankChange"],stringsAsFactors=FALSE)
plot(locfit(y~lp(x,nn=0.7), data= plotFrame,maxk=100000, kern="tcub"),
		xlab = "Normalised Std. Deviation of Importance (Mean Decrease Gini)", 
		ylab = "Normalised Ranked Order Change ( [SUM_i in 1:n-1 | SUM_j in i+1:n | abs(rank_i - rank_j)] / n(n-1)/2 )" )

dev.off()

# Now let's look at the importance plots for our sets:
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_importComapre.pdf",height=6,width=20)
par(mfcol=c(1,5))
mainNames <- c("noNormal", names(normal_fullForests))
plotThresh <- 20; aCounter <- 1
for (thisForest in c(list("noNormal"= fullForest_noNrml_importance),normal_fullForests)){
	plot(y = thisForest[1: plotThresh,"MeanDecreaseGini"], x = 1: plotThresh, xlab = "Ordered Importance", xlog = FALSE,
	ylab = "Mean Decrease Gini", type = "b", lwd = 2, pch = 15, bty = "n", main = mainNames[aCounter])
	aCounter <- aCounter + 1
	# Now we add text for the top five sites 
	text(y = thisForest[1:5,"MeanDecreaseGini"],x = 1:5, labels=rownames(thisForest)[1:5], pos = 3, offset = 1)
}
dev.off()

########################################################################################################
########################################################################################################
########################################################################################################
############################################# FINAL REVIEW #############################################
########################################################################################################
########################################################################################################
########################################################################################################

### From a personal review of the error between all my methods, as inplemented below, but found also in: compare_machineLearning.dirty.v.1.r,
### I will use the importance found from adaBoost as the cross validated method has the best error rate other than a simple tree(), which
### only suggests that 2 sites gyrA_248 and parC_2006 are important (not a big surprise).
		# This is a script to compare the statistical prediction models generated by CART/tree, randomForest and adaBoost methods,
		# These are based on runs which use the libraries "tree", "randomForest" and "adabag" and their information is produced by scripts in:
		###### DEPENDENCIES: #######
		# 	adaBoost.dirty.v.2.r	siteReview_mannualInspect_dirty.v.1.r
		
		library(tree)
		library(randomForest)
		library(adabag)
		
		# We make certain we've got a clean session
		system("rm .RData")
		rm(list=ls())
		
		
		# This hold the tree objects
		load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/tree_randomForest_objects.RData")
		# This is the randomForest object
		load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.RData")
		# This holds the adaboost object
		load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/adaBoost/adaBoost_fullResults_trimmed.RData")
		
		# This returns the total error rate as presented by a confusion table
		returnError <- function(in_confTable){
			return( round(sum(in_confTable[c(2,3)])/ sum(in_confTable),4) )
		}
		
		# Here is a quick review of the predictive qualities of each:
		# TREE objects
		table(predTree_large, phenoTest)
		print(paste("Error rate is: ", returnError(table(predTree_large, phenoTest)),sep=""))
		table(predTree_small, phenoTest)
		print(paste("Error rate is: ", returnError(table(predTree_small, phenoTest)),sep=""))
		# randomForest objects
		## NOTE: this forest was built from all the "keepSites" as defined by tree objects NOT just the top importance ones....
		table(predForest, phenoTest)
		print(paste("Error rate is: ", returnError(table(predForest, phenoTest)),sep=""))
		# adaBoost objects
		adaPred_topImp$confusion
		adaPred_topImp$error
		adaBoost_topCV$confusion
		adaBoost_topCV$error
		
		
		
		
		# Now we go and compare the importance values assigned by the randomForest vs. adaboost methods, THIS IS THE THING TO CHANGE
		topPerc <- 1
		
		giniBins <- stats::quantile(fullForest_noNrml_importance[,"MeanDecreaseGini"],seq(0,1,length.out=(100/topPerc)+1))
		topGini_sites <- fullForest_noNrml_importance[which(fullForest_noNrml_importance[,"MeanDecreaseGini"] >= giniBins[length(giniBins) - 1]),]
		topGini_sites <- topGini_sites[order(topGini_sites[,"MeanDecreaseGini"], decreasing = TRUE),]
		
		importanceBin <- stats::quantile(adaForest_importance,seq(0,1,length.out=(100/topPerc)+1))
		#importanceBin <- stats::quantile(adaForest_importance,seq(0,1,length.out=(100/topPerc)+1))
		topImportance_sites <- adaForest_importance[which(adaForest_importance >= importanceBin[length(importanceBin) -1])]
		topImportance_sites <- topImportance_sites[order(topImportance_sites, decreasing = TRUE)]
		
		# Now we cbind the siteRef information for the two approaches
		importCompare <- cbind(rownames(topGini_sites),names(topImportance_sites))
		colnames(importCompare) <- c("randomForest","adaBoost")
		round(length(intersect(importCompare[,1], importCompare[,2]))/nrow(importCompare),4)
		table(importCompare[,1]== importCompare[,2])
		which(importCompare[,1]== importCompare[,2])
		importCompare[which(importCompare[,1]== importCompare[,2]),]
		importCompare[1:25,]
		
		
		## NOTE: From my visualisation the exact rank order of top importance sites are quite different and they don't agree very rigorously on
		##		the topPerc number of sites


### Moving forward then I will use the importance determination of adaBoost and extract those sites worth preserving for further consideration
rm(list=ls())
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/adaBoost/adaBoost_fullResults_trimmed.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitesInfo_0.01.RData")
useImportance <- adaBoost_topImp$importance[order(adaBoost_topImp$importance,decreasing=TRUE)]
# We remove any importance values that are zero, this is believed to affect the analysis of data cluserting since it minimises dissimilarity
useImportance <- useImportance[which(useImportance != 0)]
rm(list=setdiff(ls(), c("useImportance","unique_sigPairs","sitesInfo")))
# Start by getting a visual "feel" for the importance values extracted
par(mfrow=c(3,1))
plot(useImportance[1:100])
plot(useImportance[1:40])
plot(useImportance[1:10])
# We will now use the PAM - partitioning around medoids function in "cluster" to determine distinct groups of importance 
# For this we'll need to know the number of clusters that exist, this is done through a silhouette calculation in HOPAC
library(cluster)
library(hopach)

findClusters <- function(inData, imprThresh = 0){
	# imprThresh is a value between 0 and 1 which informs our clustering to know how much imporvement much be generated by clustering in order
	# for a cluster to be collpased, higher values means that more clusters will be preserved, DEFAULT IS ZERO AND WORK WELL
	
	# We calculate the number of clusters that should be generated by considering euclidean distance between our data and requiring 
	# that there be at least a imprThresh amount of improvement 
	impClusters <- hopach(inData, impr = imprThresh, d = "euclid")
	retLabels <- pam(inData,k= impClusters$clustering$k ,diss = FALSE, cluster.only = TRUE)
	return( retLabels )
}

impLabels <- findClusters(useImportance)

# Now we decide which parts of our impSets to keep based on which sets do NOT contain importance values less than our keepThresh
keepThresh <- 1
keepLabels <- unique(impLabels[which(useImportance >= keepThresh)])
# Now we will create a list of the siteRef information contained in each clustered set:
impSets <- lapply(keepLabels, function(x){ return( names(impLabels)[which(impLabels == x)] ) })

# Now we can look for correlated pairs found by BayesTraits analysis which contain at least one instance of a siteRef from impSets
# This is done by looking at the siteRef_# information as found in the unique_sigPairs found in the sitesInfo_XXXXXX.RData object
# The sub('GI.',....) option is due to the way genomic island (GI) information is modified through this process by function calls.
keepSites <- unique(intersect(sub('GI.','GI:',unlist(impSets)),unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))]) ))













# This is a method to review the presence of genes of interest within our data sets, I'm going to consider the 
# three different methods for identifying site correlation with phenotype (tree, randomForest, adaBoost), as well as 
# my polymorphic alignment (to see what kind of diversity exists to have been reviewed) and BayesTraits analysis.
library(cluster)
library(hopach)
library(Vennerable)
library(RColorBrewer) # This is for colour pallete selection
rm(list=ls())
options(stringsAsFactors = FALSE)
# Tree info on all the sites used in decision trees from non-overlapping sets of sites passed to function, 
# then refined tree and randomForest objects on those "keepSites"
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/all_keepSites_mutTypes.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/tree_randomForest_objects.RData")
# randonForest on all sites without any normalisation of importance between sets  
### NOTE I HAVE TWO VERSIONS ONE INCLUDES MANY REPLCIATES AND IS THE MEAN and STDEV OF THE VALUES THIS WAS DONE TO REDUCE THE SIZE
### OF MY GINI > 1 SET AND BUILD MORE CONFIDENCE IN THESE RESULTS, THIS NEW VERSION COMES FROM: randomForest_byChunks_dirty.v.1.r
#load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.RData")
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNormal_v2_mean.RData")
	fullForest_noNrml_importance <- fullForest_noNormal_mean
	rm(fullForest_noNormal_mean)
	fullForest_noNrml_importance <- fullForest_noNrml_importance[order(fullForest_noNrml_importance[,"MeanDecreaseGini"],decreasing=TRUE),]
			
# adaBoost object
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/adaBoost/adaBoost_fullResults_trimmed.RData")
# polymorphism alignment in character and binary state
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln_traitMat_5perc_wholeAlnTree.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
# BayesTraits post-analysis information from polyAln_0.8
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitesInfo_0.01.RData")

# See my related MyAccessories script or above for info on this function.
findClusters <- function(inData, imprThresh = 0){
	impClusters <- hopach(inData, impr = imprThresh, d = "euclid")
	retLabels <- pam(inData,k= impClusters$clustering$k ,diss = FALSE, cluster.only = TRUE)
	return( retLabels )
}

# This is a function to return the overlap between sets, it requires a named vector, list, matrix or data frame type object.
# The inSets should be defined as a vector of the set identifiers we want compared.
# The byCol argument tells us the dimension to search for the inSets definitions if the data type is matrix or data.frame
# the in_ mainSet and otherSet variables should be strings that match the set references which you want to compare, where the main set
# will be compared against otherSets, these two variables can have overlap.
findOverlap <- function(inData, inSets = NULL, byCol = TRUE, in_mainSet = NULL, in_otherSet = NULL){
	# First we define what is the data type, The order here matters since matrix and data frame also apply to lower types querried here.
	tmpType <- if(is.matrix(inData)){ 
		"matrix" 
	} else if (is.data.frame(inData)) { 
		"data.frame" 
	} else if (is.list(inData)) { 
		"list" 
	} else if (length(names(inData)) > 0){
		"named"
	} else {
		print("Error: findOverlap - type not suitable or identifiable")
		return ( NULL )
	}
	# Now if the inSets were not define explicitly we inherit them as all the names of the data.
	if (is.null(inSets)){
		# This handles the data type and defines the sets  
		if (is.element(tmpType,c("matrix","data.frame"))){
			if (byCol){ inSets <- colnames(inData) } else { inSets <- rownames(inData) }
		} else if(is.element(tmpType,c("list","named"))){
			inSets <- names(inData)
		}
	}
	# This checks that we've identified at least two sets otherwise we return an error
	if (length(inSets) < 2){ 
		print("Error: findOverlap - less than 1 set was defined or identifiable from inData")
		return( NULL )
	}
	
	# Now we create vectors, using the inSets identifiers, of all the elements related to that type
	setLists <- lapply(inSets,function(thisSet){ 
		# This handles the data type and returns the objects that will make it's vector 
		if (is.element(tmpType,c("matrix","data.frame"))){
			if (byCol){ return(inData[,thisSet]) } else { return(inData[thisSet,]) }
		} else if(is.element(tmpType,c("list","named"))){
			return( unlist(inData[which(names(inData) == thisSet)]) )
		}
	})
	
	# Now let's review if the user has defined any parituclar mainSet to use, if not we return full pairwise indexes
	mainIndex <- if (is.null(in_mainSet)){
						1:(length(setLists) - 1)
					} else {
						# Otherwise the mainIndex will be those values of the set elements 
						unlist(sapply(in_mainSet,function(x){ which(inSets == x) }))
					}
	
	# Now we will create all setwise comparissons and return a list of the intersecting objects
	tmpReturn <- list()
	for(mainSet in mainIndex){
		# Now we review and define the otherSet based on user Input
		otherIndex <- if (is.null(in_otherSet)){
						(mainSet + 1):length(setLists)
					} else {
						# Otherwise the mainIndex will be those values of the set elements 
						setdiff(unlist(sapply(in_otherSet,function(x){ which(inSets == x) })), mainIndex)
					}
		# Now that we know the indexes of the otherSet we use the expand.grid function to create all possible combinations of these values
		indexSets <- expand.grid(lapply(otherIndex, function(x){ return(otherIndex) }))
		# However this will create a set that has redundancy so we'll tidy it to have unique instances
		indexSets <- unique(unlist(sapply(1:nrow(indexSets),function(x){ 
			# We turn this into a string after ordering this indexSet in order to allow us to easily call the unique function
			tmpString <- unique(unlist(indexSets[x,]))
			tmpString <- tmpString[order(tmpString)]
			return( paste(tmpString, collapse="_:_") )
			})))
		# Now we unstring the indexSets and return them as lists of indexes to be evaluated
		indexSets <- lapply(indexSets, function(x){ as.numeric(unlist(strsplit(x,"_:_"))) })
		# Now we start looking at the overlap between the mainSet and all otherSets in combination defined now by the indexSets
		for(otherSet in 1:length(indexSets)){
			tmpOverlap <- setLists[[mainSet]]
			# This is our recursion of intersections
			for(tmpSet in indexSets[[otherSet]]){
				tmpOverlap <- intersect(tmpOverlap, setLists[[tmpSet]])
			}
			# Now that we've gone through all the indexes of otherSet we will return this 
			tmpReturn[[paste(inSets[c(mainSet, indexSets[[otherSet]])],collapse="_vs_")]] <- tmpOverlap
		}
	}
	# Now we can return the intersects found
	return( tmpReturn )
}

# This is a function to return a vector of paste(collapsed) elements within columns of a matrix or data.frame()
# inData is the matrix/data.frame, inCols are the columns which you want concatenated, inSep is how they should be collapsed
extractPair <- function(inData, inCols, inSep = "__:__"){
	unique(vapply(1:nrow(inData), FUN.VALUE = vector(mode="character",length=1), FUN = function(tmpRow){
											paste(inData[tmpRow, inCols],collapse=inSep) }))
}


###############################################################################################################
########## These are the objects that ought to be changed in order to subset my analyses differently ##########
# This list will inform what my thresholds are for considering information.
keepThresh <- list("Importance"=1,"MeanGini"=1,"pValue"=1e-2,"setPerc"=0.005)
# These are the genes of interest that I want to consider
expGenes <- c("gyrA","gyrB","parC","parE","nfxB","morA")
# First will be to see what kind of representation of mutations of interest exists across my analyses
analSets <- c("tree","randomForest","adaBoost","polyAln","BayesTraits")
###############################################################################################################
###############################################################################################################
# First thing is to reduce our unique_sigPairs into only those that pass our keepThresh$pValue, but recall this is for pValue_FDR
unique_sigPairs <- unique_sigPairs[which(unique_sigPairs$pValue_FDR <= keepThresh$pValue),]
# This object uses my knowledge of each analsis set and will return the mutations which exist that are related to expGenes
if (!file.exists("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/dataSets_clust_mutPresence.RData")){
#mutPresence <- lapply(analSets, function(x){	
	mutPresence <- as.list(analSets, all.names=TRUE)
	names(mutPresence) <- analSets
	for (x in analSets){
		if(x == "tree"){
			#tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,all_keepSites_mutTypes$Site, value=TRUE) })
			#names(tmpList) <- expGenes	
			#return( tmpList )
			tmpList <- all_keepSites_mutTypes$Site
			mutPresence[[x]] <- tmpList
		} else if (x == "randomForest"){
			tmpRows <- 1:round(nrow(fullForest_noNrml_importance) * keepThresh$setPerc,0)
			tmpClust <- NULL
			# We only want to look for clusters in the highest 10% of these 400+k sites
			if(!is.element("randomForest_clust",ls())){
				tmpClust <- findClusters(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"])
				assign("randomForest_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
			} else {
				tmpClust <- randomForest_clust
			}
			tmpLabels <- unique(tmpClust[which(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"] >= keepThresh$MeanGini)])
			#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) }))) )
			mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) })))
		} else if(x =="adaBoost"){
			tmpClust <- NULL
			if(!is.element("adaBoost_clust",ls())){
				tmpClust <- findClusters(adaBoost_topImp$importance)
				assign("adaBoost_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
			} else {
				tmpClust <- adaBoost_clust
			}
			tmpLabels <- unique(tmpClust[which(adaBoost_topImp$importance >= keepThresh$Importance)])
			#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) }))) )
			mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) })))
		} else if(x=="polyAln"){
			tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,colnames(polyAln), value=TRUE) })
			names(tmpList) <- expGenes	
			#return( tmpList )
			mutPresence[[x]] <- tmpList
		} else if(x=="BayesTraits"){
			tmpList <- lapply(expGenes,function(thisGene){ 
						unlist(strsplit(grep(thisGene,unique(unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))])), value=TRUE),"::")) })
			names(tmpList) <- expGenes	
			#return( tmpList )
			mutPresence[[x]] <- tmpList
		} 
	}
	#names(mutPresence) <- analSets
	save(mutPresence, randomForest_clust, adaBoost_clust, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/dataSets_clust_mutPresence.RData")
} else {
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/dataSets_clust_mutPresence.RData")
	# We re-run building mutPresence information so that we can change the keepThresh values and extract different analyses.
	mutPresence <- as.list(analSets, all.names=TRUE)
	names(mutPresence) <- analSets
	for (x in analSets){
		if(x == "tree"){
			#tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,all_keepSites_mutTypes$Site, value=TRUE) })
			#names(tmpList) <- expGenes	
			#return( tmpList )
			tmpList <- all_keepSites_mutTypes$Site
			mutPresence[[x]] <- tmpList
		} else if (x == "randomForest"){
			tmpRows <- 1:round(nrow(fullForest_noNrml_importance) * keepThresh$setPerc,0)
			tmpClust <- NULL
			# We only want to look for clusters in the highest 10% of these 400+k sites
			if(!is.element("randomForest_clust",ls())){
				tmpClust <- findClusters(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"])
				assign("randomForest_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
			} else {
				tmpClust <- randomForest_clust
			}
			tmpLabels <- unique(tmpClust[which(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"] >= keepThresh$MeanGini)])
			#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) }))) )
			mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) })))
		} else if(x =="adaBoost"){
			tmpClust <- NULL
			if(!is.element("adaBoost_clust",ls())){
				tmpClust <- findClusters(adaBoost_topImp$importance)
				assign("adaBoost_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
			} else {
				tmpClust <- adaBoost_clust
			}
			tmpLabels <- unique(tmpClust[which(adaBoost_topImp$importance >= keepThresh$Importance)])
			#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) }))) )
			mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) })))
		} else if(x=="polyAln"){
			tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,colnames(polyAln), value=TRUE) })
			names(tmpList) <- expGenes	
			#return( tmpList )
			mutPresence[[x]] <- tmpList
		} else if(x=="BayesTraits"){
			tmpList <- lapply(expGenes,function(thisGene){ 
						unlist(strsplit(grep(thisGene,unique(unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))])), value=TRUE),"::")) })
			names(tmpList) <- expGenes	
			#return( tmpList )
			mutPresence[[x]] <- tmpList
		} 
	}	
}
# Now let's look at the overlap between my sets of analysis:  NOTE: I can use the doWeights =TRUE but it's not as visually informative.
plot(compute.Venn(Venn(list("tree"= unlist(mutPresence$tree),"randomForest"= unlist(mutPresence$randomForest),
							"adaBoost"= unlist(mutPresence$adaBoost), "polyAln"=unlist(mutPresence$polyAln), 
							"BayesTraits"=unlist(mutPresence$BayesTraits))), doWeights = FALSE,doEuler=TRUE))
# NOTE: There are adaBoost and randomForest objects not in polyAln or BayesTraits since I subset these results for expGenes.

# Here we define the overlap of siteRef between all my sets
baseOverlaps <- findOverlap(mutPresence)
#### NOW let's get a list of overlaps which show me all the BayesTraits pairs that have at least one site overlaps with correlation analysis sets.
# First we'll find which rows of unique_sigPairs have some sites that exist with overlaps on BayesTraits analysis
sigPair_overlaps <- findOverlap(list("BayesTraits"=unique(unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))])),
										#"tree"=all_keepSites_mutTypes$Site,
										"randomForest"=mutPresence$randomForest,
										"adaBoost"=mutPresence$adaBoost),in_mainSet = "BayesTraits")
sigPairs_rows <- unique(unlist(sapply(unique(unlist(sigPair_overlaps)),function(x){ union(which(unique_sigPairs[,"siteRef_1"] == x), which(unique_sigPairs[,"siteRef_2"] == x)) })))

# From experience adaBoost is worthy as it is the most refined method of my having subsetted the value of these siteRef's so we will
# consider all unique siteRef's which are correlated with some form of adaBoost analysis, 
keep_single_siteRefs <- unique(unlist(sigPair_overlaps[which(grepl("adaBoost",names(sigPair_overlaps)))]))
# Now we go and get the unique_sigPairs rows which have at least one of their siteRef values as a keep_single_siteRef
overlap_sngSite <- unique_sigPairs[unlist(sapply(keep_single_siteRefs,function(x){ union(which(unique_sigPairs$siteRef_1 == x),which(unique_sigPairs$siteRef_2 == x)) })),]
						

# Now we extract the sigPairs which have dlbSupp but have at least one site supported by multiple phenotype correlation methods
ovarlap_dblSite <- sapply(sigPair_overlaps[which(grepl("adaBoost",names(sigPair_overlaps)))],function(thisSet){ 
						# From within each set we will then return the unique_sigPairs row which corresponds to a dbl supported pair
						retFrame <- NULL
						for(thisRef in thisSet){
							# Now for each single siteRef value found in this thisSet set we will extract the rows
							# and then search if both elements in that row are elements of this thisSet
							for(thisCol in siteRef_colSet){
								tmpRows <- which(unique_sigPairs[,thisCol] == thisRef)
								if(length(tmpRows) > 0){
									for(thisRow in tmpRows){
										if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],thisSet)){
											retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
										}
									}
								}
							}
						}
						return( unique(retFrame) )  }, simplify = FALSE)

# To find sigPairs when both elements are within my analysis sets we will start with the rows which contain single siteRef's and then see
# if the other siteRef can be found
indiv_phenoSets <- list(#"tree"=mutPresence$tree,
					"randomForest"=mutPresence$randomForest,"adaBoost"=mutPresence$adaBoost)

# These are sigPairs which have at least one site supported by at least one phenotypic correlate method
sigPairs_sngSupp <- sapply(indiv_phenoSets,function(thisSet){ 
						# From within each set we will then return the unique_sigPairs row which corresponds to a dbl supported pair
						retFrame <- NULL
						for(thisRef in thisSet){
							# Now for each single siteRef value found in this thisSet set we will extract the rows
							# and then search if both elements in that row are elements of this thisSet
							for(thisCol in siteRef_colSet){
								tmpRows <- which(unique_sigPairs[,thisCol] == thisRef)
								if(length(tmpRows) > 0){
									retFrame <- rbind(retFrame,unique_sigPairs[tmpRows,])
								}
							}
						}
						return( unique(retFrame) )  }, simplify = FALSE)


# These are sigPairs which are both within the same individual phenotype correlation analysis
sigPairs_dblSupp <- sapply(indiv_phenoSets,function(thisSet){ 
						# From within each set we will then return the unique_sigPairs row which corresponds to a dbl supported pair
						retFrame <- NULL
						for(thisRef in thisSet){
							# Now for each single siteRef value found in this thisSet set we will extract the rows
							# and then search if both elements in that row are elements of this thisSet
							for(thisCol in siteRef_colSet){
								tmpRows <- which(unique_sigPairs[,thisCol] == thisRef)
								if(length(tmpRows) > 0){
									for(thisRow in tmpRows){
										if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],thisSet)){
											retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
										}
									}
								}
							}
						}
						return( unique(retFrame) )  }, simplify = FALSE)


# These are sigPairs where both sites are within at least one method but no consensus is required.  Here we should be passing a vector of all pheno correlates
# the overlapping single sites, hence the difference to some of the above
tmp_dblIndi <- sapply(unique(unlist(indiv_phenoSets)),function(thisSite){ 
						retFrame <- NULL
						# Now for each thisSite we extract the unique_sigPair rows that contain it by considering each siteRef_colSet 
						for(thisCol in siteRef_colSet){
							tmpRows <- which(unique_sigPairs[,thisCol] == thisSite)
							if(length(tmpRows) > 0){
								for(thisRow in tmpRows){
									# We look if the alternative site in this row is.element of all our individual phenotype correlations
									if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],unique(unlist(indiv_phenoSets)))){
										retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
									}
								}
							}
						}
						return( unique(retFrame) )  }, simplify=FALSE)
sigPairs_dblIndi <- NULL
for(thisSet in tmp_dblIndi){ sigPairs_dblIndi <- unique(rbind(sigPairs_dblIndi, thisSet)) }
# These are sigPairs where at least one member agreed upon by both phenotype correlation methods, but both sites are within at least
# one method.  Here we should be passing a vector of the overlapping single sites, hence the difference to some of the above
tmp_sngOver_dblIndi <- sapply(sigPair_overlaps$BayesTraits_vs_randomForest_vs_adaBoost,function(thisSite){ 
						retFrame <- NULL
						# Now for each thisSite we extract the unique_sigPair rows that contain it by considering each siteRef_colSet 
						for(thisCol in siteRef_colSet){
							tmpRows <- which(unique_sigPairs[,thisCol] == thisSite)
							if(length(tmpRows) > 0){
								for(thisRow in tmpRows){
									# We look if the alternative site in this row is.element of all our individual phenotype correlations
									if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],unique(unlist(indiv_phenoSets)))){
										retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
									}
								}
							}
						}
						return( unique(retFrame) )  },simplify=FALSE)
sigPairs_sngOver_dblIndi <- NULL
for(thisSet in tmp_sngOver_dblIndi){ sigPairs_sngOver_dblIndi <- unique(rbind(sigPairs_sngOver_dblIndi, thisSet)) }
# This is just the list of all the pairs found in unique_sigPairs
all_sigPairs = extractPair(unique_sigPairs, siteRef_colSet)

# Now let's look for sitePairs that are defined when we consider having at least one site with multiple phenot correaltion support
# and both sites being found within at least one form of phenotypic support, AND or both sites of a pair having overlapping support
### NOTE: The value of the sngOverlap being compared is it ONLY considers those single sites which have at least adaBoost support
suppPairs_list <- list("sngSupport"= unique(unlist(lapply(sigPairs_sngSupp, function(x){ extractPair(x, siteRef_colSet) }))),
						"sngOverlap"= extractPair(overlap_sngSite, siteRef_colSet),
						"dblOverlap"= unique(unlist(lapply(ovarlap_dblSite, function(x){ extractPair(x, siteRef_colSet) }))),
						"dblSupport" = unique(unlist(lapply(sigPairs_dblSupp, function(x){ extractPair(x, siteRef_colSet) }))),
						"dblIndi" = extractPair(sigPairs_dblIndi, siteRef_colSet),
						"sngOver_dblIndi" = extractPair(sigPairs_sngOver_dblIndi, siteRef_colSet))
suppPairs_overlaps <- findOverlap(suppPairs_list)



# Here we will use the protein protein interactions information to help visualise which of our protein interactions are supported.
# These are the new pa01_ & pa14_genicDB and ppInt_actions & _links.detailed objects that I can use.... NOTE: for GI information I need
# to work around a bit better as the GI of PA01 != GI: of PA14
load("/Users/Jdench/Desktop/PARuns/ppInteractions/ppInt_actions_links.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA01DBase/pa01_genicDB.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA14DBase/pa14_genicDB.RData")
# CHALLENGE IS THAT THE ppInt INFROMATION IS FOR THE REFERENCE STRAIN PA01 BUT ALL MY ANALYSIS HAS USED THE PA14 REFERENCE
# THIS MEANS WE WILL HAVE TO USE BLAST TO TRY AND FIND THE SEQUENCE ANALOGS BETWEEN GI: SEQUENCES.
######## BLAST definitions #############
# This is the base directory for blast
baseBlast <- paste("/Users/Jdench/Desktop/Blast/ncbi-blast-2.2.30+/",sep="")
# This is the explicit location in which we will find the blast programs we might want to run
blastLoc <- paste(baseBlast,"bin/",sep="")
# This is the location of blast databases and where we will write our temporary ones
blastDB <- paste(baseBlast,"db/Temp/",sep="")
library(seqinr)

# This is a function to perform a recursive intersection call on a list of information
recursiveIntersect <- function(inList){
	if(length(inList) > 1){
		retIntersect <- inList[[1]]
		for(x in 2:length(inList)){
			retIntersect <- intersect(retIntersect,inList[[x]])
		}
		return( retIntersect )
	} else {
		# If there is only one elements than we return that list element as it intersects itself by default.
		return( inList )
	}
}

# We will extract those of our proposed pairs we want to explore and look if there are ppInt that exist
# We will use all unique sngOverlap, bldOverlap and dblSupport sites
ppInt_pairs <- unlist(lapply(unique(unlist(suppPairs_list)),function(thisPair){
					# We only really need the gene information not the actual site of the mutation.
					tmpPair <- unlist(lapply(strsplit(thisPair,"__:__")[[1]],function(x){ substr(x,1,regexpr("_",x)-1) }))
					# We now look if there are PA###### Locus_Tag's identifiable for all elements of thisPair
					tmpLoci <- unlist(sapply(tmpPair,function(x){ 
						# We want to handle if nothing is found by returning an NA value
						tmpIndex <- which(pa01_genicDB$geneID == x)
						if(length(tmpIndex) == 0){
							return( NA )
						} else {
							return( pa01_genicDB$Locus_Tag[which(pa01_genicDB$geneID == x)] )
						} })) 
					# Now we handle circumstance where nothing was found bu geneID reference by using a blast search of the sequence
					if(is.element(NA,tmpLoci)){
						# First we check if there is a pa01_genicDB database else we make one
						if (!file.exists(paste(blastDB,"pa01_genicDB.nsq",sep=""))){
							print(paste("Missing the blastdb for pa01_genicDB trying to create it",sep=" "))
							# We need to make certain there is a Seq_nt sequence otherwise we won't include that to our blastDB
							tmpIndex_seqNT <- which(nchar(pa01_genicDB$Seq_nt) > 0)
							write.fasta(as.list(pa01_genicDB$Seq_nt[tmpIndex_seqNT]),names=pa01_genicDB$Locus_Tag[tmpIndex_seqNT],file=paste(blastDB,"pa01_genicDB.fas",sep=""))
							system(paste(blastLoc,"makeblastdb -in ",blastDB,"pa01_genicDB.fas -dbtype nucl -input_type fasta -out ",blastDB,"pa01_genicDB",sep=""))
							unlink( paste(blastDB,"pa01_genicDB.fas",sep="") )
						}
						# Now we will blast the Seq_nt for this gene from pa14 to the pa01 db
						for(thisGene in tmpPair[which(is.na(tmpLoci))]){
							# We hand the ENTREZ values as some instances of my siteRef's may not have the same caps values, this shouldn't happen for anything other than Entrez: values.
							tmpSeq <- pa14_genicDB$Seq_nt[which(toupper(pa14_genicDB$geneID)== toupper(thisGene))]
							tmpFile <- paste(blastDB,"tmp_",gsub("[[:punct:]]","_",thisGene),".fas",sep="")
							# If we found a sequence then we will perform a blast search for this
							write.fasta(as.list(tmpSeq),names=thisGene,file=tmpFile)
							system(paste(blastLoc,"blastn -query ", tmpFile," -db ",blastDB,"pa01_genicDB -out ", blastDB,"tmpQuery_blast.table -outfmt '6 sseqid sstart send qstart qend qlen evalue pident sseq'",sep=""))
							# Now we look if we found anything by searching is a _blast.table file exists and if it has content
							if(file.exists(paste(blastDB,"tmpQuery_blast.table",sep="")) && length(readLines(paste(blastDB,"tmpQuery_blast.table",sep=""))) > 0){
								# We read in the blast.table information and run a bit of quality control by verifying the percIdent value and eValue
								tmpFrame <- read.table(paste(blastDB,"tmpQuery_blast.table",sep=""),header=FALSE, stringsAsFactors = FALSE)
								# These colnames are based on my knowledge of how I am making the blast call
								colnames(tmpFrame) <- c("refID","refStart","refEnd","qStart","qEnd","qLen","eValue","PercIdent","refSeq")
								if(tmpFrame$eValue <= 1e-40 || tmpFrame$PercIdent >= 90){
									tmpLoci[thisGene] <- tmpFrame$refID[union(which(tmpFrame$eValue == min(tmpFrame$eValue)),which(tmpFrame$PercIdent == max(tmpFrame$PercIdent)))[1]]
								}
								# We no longer need this tmpQuery_blast.table information
								unlink( paste(blastDB,"tmpQuery_blast.table",sep="") )
							}
							# we also no longer need the tmp_thisGene.fas file
							unlink( tmpFile )	
						}
					}
					# Now in theory we have a PA01 locus_tag value for each element of tmpPair so we look if there is a ppInt using the action
					# information since it seems that this is a subset of all possible ppInt, and only those with some activity support.
					# NOTE: In using the links.detailed information I could assume that a combined score being greater than 400 is meaningfull, from experience NOT all
					# possible protein-protein interactions exist in this object so we could even have performed a raw exists check (I did and got limited return values)
					# I chose not to use a cut-off for this combined score since my objective is to find possibly poorly characterised pairs of mutations, thus we should not prejudice our
					# ppInt found pairs via a requirement of strong supporting information
					
					#  We check that there are no more NA values in the tmpLoci
					if(!is.element(NA,tmpLoci)){
						#tmpList <- lapply(tmpLoci, function(x){ union(which(grepl(x,ppInt_actions$item_id_a)),
						#																				which(grepl(x, ppInt_actions$item_id_a))) })																										
						tmpList <- lapply(tmpLoci, function(x){ union(which(grepl(x,ppInt_links.detailed$protein1)),
																										which(grepl(x, ppInt_links.detailed$protein2))) })
						tmpIntersect <- recursiveIntersect(tmpList)																					
						# If we have found any tmp_intersectAction than we will return thisPair, else we do not
						#if(length(tmp_intersectAction) > 0){ 
						if(length(tmpIntersect) > 0){ 
							# This is a visualising sanity check 
							#print(thisPair)
							#print(ppInt_links.detailed[tmpIntersect,])
							# Return thisPair since we've found an intersect of protein protein interactions.	
							return( thisPair ) 
						}
					}
					# If we can't find PA##### values for both of thisPair's genes OR we can;t find a tmpIntersection then we return nothing, which defaults to NULL
					# Which itself will be removed when we unlist this lapply call
	}))


######## END OF MY BASIC OVERLAP LISTS ALL BELOW MAY BE MORE SPECIALISED #########






################################## FINAL PLOTS ##################################
###### These are the final Venn plots that I will want to use for convincing Rees and maybe as SI in publication that my choices were sound.
# These first two are meant to show that our methods of analysis provided numerous choices:
# 1) The first is all possible BayesTraits sites and all Phenotype correlation sites
all_siteRef <- list("BayesTraits" = unique(unlist(unique_sigPairs[, siteRef_colSet])),
					"Phenotype" = unique(unlist(indiv_phenoSets)) )
tmpVenn <- compute.Venn(Venn(all_siteRef), doWeights = TRUE,doEuler=TRUE, type = "circles")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/final_1_sitePairs_all_siteRefs_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()
# 2) The second is to show that for simply considering sites that are correlated to resistance phenotype there were many options.
tmpVenn <- compute.Venn(Venn(indiv_phenoSets), doWeights = TRUE,doEuler=TRUE, type = "circles")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/final_2_siteRefs_indiv_phenoSets_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()
# The last two are to show how my choices were made, this will show that of the choices I made:
# 3) There were some sites in expected genes but we were also able to define novel sites that could be correlated for either the phenotype or through evolution (compensatory then?)
	# This is for single siteRef #final3_list <- list("sng_phenoOverlap"=unique(unlist(strsplit(suppPairs_list$sngOverlap,"__:__"))), "expGene_polySite"=unlist(mutPresence$polyAln))
final3_list <- list("sng_phenoOverlap"= unique(suppPairs_list$sngOverlap), "dbl_phenoSupport"= suppPairs_list$dblSupport, "expGene_polySite"= unique(unlist(lapply(expGenes,function(x){ suppPairs_list$allPairs[which(grepl(x,suppPairs_list$allPairs))] }))))
tmpVenn <- compute.Venn(Venn(final3_list), doWeights = TRUE,doEuler=TRUE, type = "circles")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/final_3_phenoOverlap_expPoly_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()
# 4) This shows that we have a small refined set which are built from having single overlapping pheno support, and double sites pheno detection
final4_list <- list("sng_phenoOverlap"= unique(suppPairs_list$sngOverlap), "dbl_phenoSupport"= suppPairs_list$dblIndi, "ppInt_support"= ppInt_pairs)
tmpVenn <- compute.Venn(Venn(final4_list), doWeights = TRUE ,doEuler=TRUE, type = "circles")#, type = "squares")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/final_4_sngOverlap_dblSupport_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()


# Now we will finally be using can be defined by finding the overlap between ppInt and other support mechanism:
keepPairs <- findOverlap(suppPairs_list)
save(keepPairs, suppPairs_list, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/keepPairs.RData")
## FROM REVIEW OF THE keepPairs list and final4_list Venn plot, I will consider valuable any of those sitePairs which have ppInt support,
## and will consider building a complete empirical fitness landscape from the site pairs with overlap between all three refinement methods
## and will include some "usual suspects" to simply explore the story a bit more thoroughly.  See my lab manual (ringed 3 subject Hillary notebook, on 28-03-2016)

# This last plot is of the pairs that I will use and is just a means of visualising the alignment for the sites and site pairs being considered.
load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")	
kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
if (is.element(x, c("PA14","PA01"))){ 
	return( "ND" )
} else {
	as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
} } ))
# We define the rows worth considering and setup some plotting information to be used downstream, we will keep the pairs that have ppInt support and double
# phnoetypic support PLUS one site pair which we are proposing to be part of a complete empirical fitness landscape we will construct - "rpoB_864__:__gyrA_259"
keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
plotSites <- list("sngSites"=unique(unlist(strsplit(c(keepPairs$dbl_phenoSupport_vs_ppInt_support,"rpoB_864__:__gyrA_259"), "__:__"))),
					"sitePairs"= unlist(strsplit(c(keepPairs$dbl_phenoSupport_vs_ppInt_support,"rpoB_864__:__gyrA_259"), "__:__")) )
plotLength <- max(sapply(plotSites,length))
numPlots <- 2
ntNames <- c("A","C","G","T","-")
phenoNames <- c("(R)","(S)")
# I'm pallate of colours from a qualitative set that by knowledge accepts the length of my set (max = 9)
# We set the phenoCols first as from visualising the set we know that blue and red are the first two choices
allCols <- brewer.pal(n= length(c(phenoNames,ntNames)),name= "Set1")
names(allCols) <- c(phenoNames,ntNames)

# This is the phenotype alignments for the sites of interest which I'm considering now as the "sngOverlap_vs_dblSupport" group
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/final_5_keepPairs_alignment.pdf",height=32,width= numPlots * max(4,ceiling(plotLength/5)))
par(mfcol = c(1, numPlots))
for(thisPlot in plotSites){
	# This first plot is of the 
	plotMat <- cbind(kos_phenoVec[keepRows],toupper(polyAln[keepRows, thisPlot]))
	# Now we create a blank plot space and then will simply plot coloured pch squares based on the phenotype and nt types
	plot(x = NULL, y = NULL, xlim = c(1, plotLength + 1), ylim = c(1, nrow(plotMat)), ylab = "Phenotype and Nucleotide Alignment", xlab = "", axes = FALSE)
	# We use the colour sets and their names to be how we will search through this, these are not robust since the lists are manually generated.
	for(thisState in names(allCols)){
		tmpPlot <- which(plotMat == thisState,arr.ind = TRUE)
		points(x = tmpPlot[,"col"], y = tmpPlot[,"row"], pch = 15, cex = 1, col = allCols[thisState])
	}
	# We will add text to the bottom of the plot to indicate what information is being shown
	text(x = 0: length(thisPlot) + 1, y = 1, labels = c("Phenotype", thisPlot) ,pos = 1, offset = 2, srt = 90)
}
legend("bottom",cex=1 ,names(allCols),bty="n",fill= allCols, ncol = 4)
dev.off()





### Here we look at the keepPairs that have been identified across a suite of sets and we will look at the qValues that have been accepted
### This is a raw means to considering which site pairs to move forward with since it's been asked that I increase the stringency of my 
### randomForest cutoff so that the sets appear to more comparable.
rm(list=ls())
# This is a function to return a vector of paste(collapsed) elements within columns of a matrix or data.frame()
# inData is the matrix/data.frame, inCols are the columns which you want concatenated, inSep is how they should be collapsed
extractPair <- function(inData, inCols, inSep = "__:__"){
	unique(vapply(1:nrow(inData), FUN.VALUE = vector(mode="character",length=1), FUN = function(tmpRow){
											paste(inData[tmpRow, inCols],collapse=inSep) }))
}

# This is the sites in our unique_sigPairs data.frame that we want to search for relation.
siteRef_colSet <- c("siteRef_1","siteRef_2")

# These will be our unique_sigPairs
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitesInfo_0.01.RData")
# We extract all the pairs in our unique_sigPairs
allPairs <- extractPair(unique_sigPairs, siteRef_colSet)

# Here chose the sets that I want to include
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/keepPairs_MeanGini_1_pVal_minus2/keepPairs.RData")
assign("keepPairs_gini1p2",keepPairs,pos=".GlobalEnv")
gini1p2_sigPairs <- unique_sigPairs[unlist(lapply(keepPairs_gini1p2$sng_phenoOverlap_vs_ppInt_support, function(x){ which(allPairs == x)})),]

load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/keepPairs_MeanGini_2_pVal_minus2/keepPairs.RData")
assign("keepPairs_gini2p2",keepPairs,pos=".GlobalEnv")
gini2p2_sigPairs <- unique_sigPairs[unlist(lapply(keepPairs_gini2p2$sng_phenoOverlap_vs_ppInt_support, function(x){ which(allPairs == x)})),]

load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/keepPairs_MeanGini_2_pVal_minus3/keepPairs.RData")
assign("keepPairs_gini2p3",keepPairs,pos=".GlobalEnv")
gini2p3_sigPairs <- unique_sigPairs[unlist(lapply(keepPairs_gini2p3$sng_phenoOverlap_vs_ppInt_support, function(x){ which(allPairs == x)})),]

# Now we just plot these three as histograms to see where they lay.
pdf("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/compare_gini1_2_qVals.pdf",height=6,width=12)
par(mfcol = c(1,3))
hist(gini1p2_sigPairs$pValue_FDR)
hist(gini2p2_sigPairs$pValue_FDR)
hist(gini2p3_sigPairs$pValue_FDR)
dev.off()










####################################################################################################################################
###################################################### OLD BUT NOT UNFAIR ##########################################################
####################################################################################################################################

# FROM HERE I SHOULD EXTRACT THE ACTUAL unique_sigPairs ROWS AND KEEP PERFORMING ANALYSIS WITH THESE #
# I THINK IT MAKES MORE SENSE TO SIMPLY EXTRACT ALL THE siteRef VALUES FOUND IN EACH ANALYSIS SET AND THEN GET OVERLAPS #




# So for all rows which had a site that was identified in our phenotype correlation analyses we will want to find overlaps with the site pair, 
# the strplit is to handle how I've merged redundant siteRef information in upstream application
all_sitePairs <- lapply(sigPairs_rows,function(x){ return( unlist(strsplit(unlist(unique_sigPairs[x,c("siteRef_1","siteRef_2")]),"::")) ) })
names(all_sitePairs) <- paste("row", sigPairs_rows,sep="")
all_sitePairs_overlaps <- lapply(all_sitePairs, function(x){ findOverlap(list("sitePairs"=unlist(x),"tree"=all_keepSites_mutTypes$Site,"randomForest"=mutPresence$randomForest,"adaBoost"=mutPresence$adaBoost), in_mainSet = "sitePairs", in_otherSet = c("tree","randomForest","adaBoost")) })
	
# Now we can compare which of these sitePairs has the most support, in a raw sense this is the sets with the longest "length"
overlapLength <- sapply(all_sitePairs_overlaps, function(x){ length(unlist(x)) })
table(overlapLength)
# Same but we are looking for the number of pair overlaps
# I CAN BE SMARTER ABOUT HOW I CHOOSE WHAT I KEEP....
# So there are two trains of thought to be had here, 

 #I) We will look at those sites which have been supported by at least two different approaches to find correlation
considered_analSets <- c("sitePairs_vs_tree_vs_randomForest","sitePairs_vs_tree_vs_adaBoost",
				"sitePairs_vs_tree_vs_randomForest_vs_adaBoost","sitePairs_vs_randomForest_vs_adaBoost")
# We will split our summary here by whether there was dual or tripple support:
suppStrengths <- list("dual"= unlist(lapply(considered_analSets, function(x){ if(length(gregexpr("_vs_",x)[[1]]) == 2){ return(x) } })), 
						"triple"= unlist(lapply(considered_analSets, function(x){ if(length(gregexpr("_vs_",x)[[1]]) == 3){ return(x) } })))
# We now extract the unique_sigPairs rows which relate to having dual and triple support
keep_multiSupp <- sapply(names(suppStrengths), function(thisStr){ 
	retData <- NULL
	for(thisAnal in suppStrengths[[thisStr]]){
		tmpRows <- as.numeric(unlist(sapply(1:length(all_sitePairs_overlaps), function(thisRow){ 
			if(length(all_sitePairs_overlaps[[thisRow]][[thisAnal]]) > 0){ 
				return( sub("row","",names(all_sitePairs_overlaps)[thisRow]) ) 
			} })))
		retData <- rbind(retData, unique_sigPairs[tmpRows,])
	}
	return( retData )
	}, simplify = FALSE)

#II) we can either look for when both sites are correlated to the resistance phenotype by a considered analysis method
#### NOTE: One thing to keep in mind is that tree() alone is not a very meaningfull set to consider
# Therefore we should only consider when there is support for a solo outside of tree(), or as well as tree()
# This is hard coded for this analysis
considered_analSets <- c("sitePairs_vs_tree_vs_randomForest","sitePairs_vs_tree_vs_adaBoost",
				"sitePairs_vs_tree_vs_randomForest_vs_adaBoost","sitePairs_vs_randomForest","sitePairs_vs_randomForest_vs_adaBoost",
				"sitePairs_vs_adaBoost")
paired_overlapLength <- sapply(all_sitePairs_overlaps, function(x){ length(which(sapply(considered_analSets, function(y){ length(x[[y]]) }) == 2)) })
table(paired_overlapLength)
keep_dblPairs <- unique_sigPairs[as.numeric(sub("row","",names(all_sitePairs_overlaps)[which(paired_overlapLength >= 1)])),]

##### NOW LETS LOOK AT THE SITES THAT HAVE DBL AND TRIPLE SUUPORT, AND FOR INTERSECTS WITH THE keep_dblPairs set
keepSets_pairs <- list("dual"= sapply(1:nrow(keep_multiSupp$dual),function(x){ paste(keep_multiSupp$dual[x,c("siteRef_1","siteRef_2")], collapse="__:__") }),
						"triple"= sapply(1:nrow(keep_multiSupp$triple),function(x){ paste(keep_multiSupp$dual[x,c("siteRef_1","siteRef_2")], collapse="__:__") }),
						"dblPairs"= sapply(1:nrow(keep_dblPairs),function(x){ paste(keep_multiSupp$dual[x,c("siteRef_1","siteRef_2")], collapse="__:__") }) )

plot(compute.Venn(Venn(keepSets_pairs), doWeights = FALSE,doEuler=TRUE))
sitePair_overlaps <- findOverlap(keepSets_pairs, in_mainSet = "dblPairs", in_otherSet = c("dual","triple"))
######################################################################################################################################################## 
######## THIS IS CHOICES BASED ON USING SITES WITH PRIORITY BASED ON HAVING BOTH BAYESTRAITS siteRef VALUES FOUND BY PHENO CORRELATION ANALYSIS #######
######################################################################################################################################################## 
# I've found that all my double pairs are themselves found in at least two or three analysis sets
# This is a sanity check to ensure that I've not made a mistake in the intermediate funneling of information
unique_sigPair_siteRef_strings <- unlist(lapply(1:nrow(unique_sigPairs),function(x){ paste(unique_sigPairs[x,c("siteRef_1","siteRef_2")], collapse="__:__") }))

keepSets_pairRows <- sapply(keepSets_pairs,simplify = FALSE, FUN = function(x){
	unlist(lapply(x,function(thisPair){ which(unique_sigPair_siteRef_strings == thisPair) })) })
# This checks that the unique_sigPairs elements are what we've defined them as being
sanityCheck <- NULL
sanityCheck <- c(sanityCheck, all(sapply(1:length(keepSets_pairRows),function(x){ 
	all(sapply(1:length(keepSets_pairs[[x]]), function(y){ 
		keepSets_pairs[[x]][y] == paste(unique_sigPairs[keepSets_pairRows[[x]][y],c("siteRef_1","siteRef_2")], collapse="__:__") })) })) )
# This checks that the pairs we're reporting are really those we've been able to find in the analyses, '
# we've defined them as being site pairs both found in at least two sets of analysis, therefore:
sanityCheck <- c(sanityCheck, all(sapply(keepSets_pairs$dblPairs, function(x){ is.element(x, unique_sigPair_siteRef_strings) })) )

if (all(sanityCheck)){
	# I am keeping this as I've found that the dblPairs are in dual and triple correlation analysis sets
	keep_sigPairs <- keep_dblPairs
	save(keep_sigPairs, file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/keep_sigPairs.RData")
} else {
	print("My sanity check has shown as problem with confirming my analyses having proper overlaps")
}

#### NEXT will be to look at the alignment of my sites of interest and make certain I'm satisfied with these results, get the actual character alignments visualised.
# We load some data to allow us to have the phenotype vector
load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")	
kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
if (is.element(x, c("PA14","PA01"))){ 
	return( "ND" )
} else {
	as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
} } ))
# We define the rows worth considering
keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
# Then we will create a list of all the siteRef's we are considering, this is from a keep_sigPairs object made above
keep_siteRef <- lapply(1:nrow(keep_sigPairs),function(x){ unlist(strsplit(unlist(keep_sigPairs[x,c("siteRef_1","siteRef_2")]),"::")) })
# We know the number of sites based on the length of the unlisted keep_siteRef list
plotLength <- length(unique(unlist(keep_siteRef)))
# these are plotting colours for our grid plot
ntCols <- list("A"="orange","C"="green","G"="purple","T"="yellow","-"="black")
phenoCols <- list("(S)"="blue","(R)"="red")
# Now we'll use the top X most interesting sites and look at a "plot" of their phenotype's by nucleotide states nicely coloured
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/keep_sigPairs_phenoAlign_v2.pdf",height=32,width=max(4,ceiling(plotLength/5)))
plotMat <- cbind(kos_phenoVec[keepRows],toupper(polyAln[keepRows,unique(unlist(keep_siteRef))]))
# Now we create a blank plot space and then will simply plot coloured pch squares based on the phenotype and nt types
plot(x = NULL, y = NULL, xlim = c(1, plotLength + 1), ylim = c(1, nrow(plotMat)), ylab = "Phenotype and Nucleotide Alignment", xlab = "", axes = FALSE)
# We use the colour sets and their names to be how we will search through this, these are not robust since the lists are manually generated.
for(thisType in c("phenoCols","ntCols")){
	for(thisState in names(eval(as.name(thisType))) ){
		tmpPlot <- which(plotMat == thisState,arr.ind = TRUE)
		points(x = tmpPlot[,"col"], y = tmpPlot[,"row"], pch = 15, cex = 1, col = eval(as.name(thisType))[[thisState]])
	}
}
# We will add text to the bottom of the plot to indicate what information is being shown
text(x = 0: plotLength + 1, y = 1, labels = c("Phenotype",unique(unlist(keep_siteRef))) ,pos = 1, offset = 2, srt = 90)
legend("bottom",cex=1 ,c(names(phenoCols),names(ntCols)),bty="n",fill= unlist(list(phenoCols, ntCols)), ncol = 4)
dev.off()
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 





####### FROM SOME REVIEW THE suppPairs_overlaps$sngOverlap_vs_dblSupport group might be most worthy of future work. #######
# NOTE: this takes quite a bit of time due to the number of allPairs that exist.
vennPlot <- compute.Venn(Venn(suppPairs_list), doWeights = FALSE,doEuler=TRUE, type = "ellipses")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sngOverlap_vs_dblSupport_Venn.pdf",height=8,width=8)
plot(vennPlot)
dev.off()

# Now let us also plot the overlap with expGene sites that exist among BayesTraits analysis, phenoType significance and my final set.
# This should plot as a series of nested positions where the polyAln and indiv_phenoSites are whole sets for the related element.
expGene_forVenn <- list("all_polyAln" = unlist(mutPresence$polyAln), "all_BayesTraits" = unique(unlist(mutPresence$BayesTraits)),
						"indiv_phenoSites" = unique(unlist(indiv_phenoSets)),
						"suppPairs_overlaps$sngOverlap_vs_dblSupport" = unique(unlist(strsplit(suppPairs_overlaps$sngOverlap_vs_dblSupport,"__:__"))) )
vennPlot <- compute.Venn(Venn(expGene_forVenn), doWeights = FALSE,doEuler=TRUE, type = "ellipses")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/all_siteRefs_Venn.pdf",height=8,width=8)
plot(vennPlot)
dev.off()

# To allow me to have some weighted Venn diagrams I'll plot sets of three or less sets to be compared
# The first will be a look at the sitePairs we will consider compared to all possible sitePairs, and those which have some pheno correlation
tmpPlot_list<- list("all_strongPairs"= extractPair(unique_sigPairs[which(unique_sigPairs$pValue_FDR <= 1e-4),], siteRef_colSet),
					"sngOverlap"= suppPairs_list$sngOverlap,
					"suppPairs_overlaps$sngOverlap_vs_dblSupport" = suppPairs_overlaps$sngOverlap_vs_dblSupport)
tmpVenn <- compute.Venn(Venn(tmpPlot_list), doWeights = FALSE,doEuler=TRUE, type = "circles")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitePairs_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()

# Next we will show how many of the pairs have at least one expGene site, both siteRefs in expGenes, and all sitePairs with at least one predicted site
exp_indiv_siteRef <- unique(unlist(lapply(expGenes,function(x){ unlist(indiv_phenoSets)[which(grepl(x,unlist(indiv_phenoSets)))] })))
exp_indiv_frame <- unique_sigPairs[unlist(sapply(exp_indiv_siteRef,function(x){ union(which(unique_sigPairs$siteRef_1 == x),which(unique_sigPairs$siteRef_2 == x)) })),]

expGene2_forVenn <- list("expSng_keepPairs"= unique(unlist(lapply(expGenes,function(x){ grep(x,suppPairs_overlaps$sngOverlap_vs_dblSupport,value=TRUE) }))),
							"expDbl_keepPairs" = unique(unlist(lapply(1:length(suppPairs_overlaps$sngOverlap_vs_dblSupport),function(x){ 
								tmpSet <- strsplit(suppPairs_overlaps$sngOverlap_vs_dblSupport[x],"__:__")[[1]]
								if(length(setdiff(tmpSet,unique(unlist(sapply(expGenes,function(y){ grep(y, tmpSet,value=TRUE) }))))) == 0){ 
									return(suppPairs_overlaps$sngOverlap_vs_dblSupport[x]) 
								} }))),
							"expSng_phenoPairs"= extractPair(exp_indiv_frame, siteRef_colSet) )
							 
tmpVenn <- compute.Venn(Venn(expGene2_forVenn), doWeights = FALSE,doEuler=TRUE, type = "circles")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitePairs_expGenes2_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()							

# Now let's plot all possible BayesTraits pairs vs. all possible phenotype Pairs
all_sitePairs <- list("BayesTraits"= suppPairs_list$allPairs,
						"Phenotype"= extractPair(unique_sigPairs[unlist(sapply(unlist(indiv_phenoSets),function(x){ union(which(unique_sigPairs$siteRef_1 == x),which(unique_sigPairs$siteRef_2 == x)) })),], siteRef_colSet) )
tmpVenn <- compute.Venn(Venn(all_sitePairs), doWeights = FALSE,doEuler=TRUE, type = "circles")
pdf(file="/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitePairs_allPairs_Venn.pdf",height=8,width=8)
plot(tmpVenn)
dev.off()



















######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
######################################################################################################################################################## 
# This is just a bunch of my "siteReview_mannualInspect_dirty.v.3.r" code which I've taken out and am cleaning up a bit so that I can run
# more batch like jobs.


# This is a method to review the presence of genes of interest within our data sets, I'm going to consider the 
# three different methods for identifying site correlation with phenotype (tree, randomForest, adaBoost), as well as 
# my polymorphic alignment (to see what kind of diversity exists to have been reviewed) and BayesTraits analysis.
rm(list=ls())

library(cluster)
library(hopach)
library(Vennerable)
library(igraph) # This is to visualise network of connections
library(RColorBrewer) # This is for colour pallete selection
options(stringsAsFactors = FALSE)
# Tree info on all the sites used in decision trees from non-overlapping sets of sites passed to function, 
# then refined tree and randomForest objects on those "keepSites"
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/all_keepSites_mutTypes.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/tree_randomForest_objects.RData")
# randonForest on all sites without any normalisation of importance between sets 
	#### NOTE: I now have two versions of the randomForest assessment
	#load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/fullForest_noNrml_importance.RData")
	load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/randomForest/fullForest_noNormal_v2_mean.RData")
	fullForest_noNrml_importance <- fullForest_noNormal_mean
	rm(fullForest_noNormal_mean)
	fullForest_noNrml_importance <- fullForest_noNrml_importance[order(fullForest_noNrml_importance[,"MeanDecreaseGini"],decreasing=TRUE),]			
# adaBoost object
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_wholeGenomes/newTrials/adaBoost/adaBoost_fullResults_trimmed.RData")  
# polymorphism alignment in character and binary state
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln_traitMat_5perc_wholeAlnTree.RData")
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/polyAln.RData")
# BayesTraits post-analysis information from polyAln_0.8
load("/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/sitesInfo_0.01.RData")

# These are sets for the threshholds over which we will analyse our data
giniSet <- c(1,1.5)
impSet <- c(1)
pVal_set <- c(1e-2, 1e-3)


# See my related MyAccessories script or above for info on this function.
findClusters <- function(inData, imprThresh = 0){
	impClusters <- hopach(inData, impr = imprThresh, d = "euclid")
	retLabels <- pam(inData,k= impClusters$clustering$k ,diss = FALSE, cluster.only = TRUE)
	return( retLabels )
}


##### OK HERE IS WHERE IS SPLIT TO MAKE SOME MORE BATCH TYPE WORK POSSIBLE ####
#### DEFINE A SUITE OF KEEPTHRESH VALUES TO LOOP OVER
baseDir <- "/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8/postAnal/" 

for(thisGini in giniSet){
	for(this_pVal in pVal_set){
		for(thisImportance in impSet){
			# We now set the working directory based on these values
			workDir <- paste(baseDir,"keepPairs_Gini_", thisGini,"_Import_",thisImportance,"_pVal_", this_pVal,"/",sep="")
			if(!file.exists(workDir)){ system(paste("mkdir ",workDir,sep="")) }
			###############################################################################################################
			########## These are the objects that ought to be changed in order to subset my analyses differently ##########
			# This list will inform what my thresholds are for considering information.
			keepThresh <- list("Importance"= thisImportance,"MeanGini"= thisGini,"pValue"= this_pVal,"setPerc"=0.005)
			# These are the genes of interest that I want to consider
			expGenes <- c("gyrA","gyrB","parC","parE","nfxB","morA", "grlA")
			# First will be to see what kind of representation of mutations of interest exists across my analyses
			analSets <- c("tree","randomForest","adaBoost","polyAln","BayesTraits")
			###############################################################################################################
			###############################################################################################################
			# First thing is to reduce our unique_sigPairs into only those that pass our keepThresh$pValue, but recall this is for pValue_FDR
			unique_sigPairs <- unique_sigPairs[which(unique_sigPairs$pValue_FDR <= keepThresh$pValue),]
			# This object uses my knowledge of each analsis set and will return the mutations which exist that are related to expGenes
			if (!file.exists(paste(baseDir,"dataSets_clust_mutPresence.RData",sep=""))){
			#mutPresence <- lapply(analSets, function(x){	
				mutPresence <- as.list(analSets, all.names=TRUE)
				names(mutPresence) <- analSets
				for (x in analSets){
					if(x == "tree"){
						#tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,all_keepSites_mutTypes$Site, value=TRUE) })
						#names(tmpList) <- expGenes	
						#return( tmpList )
						tmpList <- all_keepSites_mutTypes$Site
						mutPresence[[x]] <- tmpList
					} else if (x == "randomForest"){
						tmpRows <- 1:round(nrow(fullForest_noNrml_importance) * keepThresh$setPerc,0)
						tmpClust <- NULL
						# We only want to look for clusters in the highest 10% of these 400+k sites
						if(!is.element("randomForest_clust",ls())){
							tmpClust <- findClusters(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"])
							assign("randomForest_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
						} else {
							tmpClust <- randomForest_clust
						}
						tmpLabels <- unique(tmpClust[which(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"] >= keepThresh$MeanGini)])
						#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) }))) )
						mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) })))
					} else if(x =="adaBoost"){
						tmpClust <- NULL
						if(!is.element("adaBoost_clust",ls())){
							tmpClust <- findClusters(adaBoost_topImp$importance)
							assign("adaBoost_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
						} else {
							tmpClust <- adaBoost_clust
						}
						tmpLabels <- unique(tmpClust[which(adaBoost_topImp$importance >= keepThresh$Importance)])
						#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) }))) )
						mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) })))
					} else if(x=="polyAln"){
						tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,colnames(polyAln), value=TRUE) })
						names(tmpList) <- expGenes	
						#return( tmpList )
						mutPresence[[x]] <- tmpList
					} else if(x=="BayesTraits"){
						tmpList <- lapply(expGenes,function(thisGene){ 
									unlist(strsplit(grep(thisGene,unique(unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))])), value=TRUE),"::")) })
						names(tmpList) <- expGenes	
						#return( tmpList )
						mutPresence[[x]] <- tmpList
					} 
				}
				#names(mutPresence) <- analSets
				save(mutPresence, randomForest_clust, adaBoost_clust, file=paste(baseDir,"dataSets_clust_mutPresence.RData",sep=""))
			} else {
				load(paste(baseDir,"dataSets_clust_mutPresence.RData",sep=""))
				# We re-run building mutPresence information so that we can change the keepThresh values and extract different analyses.
				mutPresence <- as.list(analSets, all.names=TRUE)
				names(mutPresence) <- analSets
				for (x in analSets){
					if(x == "tree"){
						#tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,all_keepSites_mutTypes$Site, value=TRUE) })
						#names(tmpList) <- expGenes	
						#return( tmpList )
						tmpList <- all_keepSites_mutTypes$Site
						mutPresence[[x]] <- tmpList
					} else if (x == "randomForest"){
						tmpRows <- 1:round(nrow(fullForest_noNrml_importance) * keepThresh$setPerc,0)
						tmpClust <- NULL
						# We only want to look for clusters in the highest 10% of these 400+k sites
						if(!is.element("randomForest_clust",ls())){
							tmpClust <- findClusters(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"])
							assign("randomForest_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
						} else {
							tmpClust <- randomForest_clust
						}
						tmpLabels <- unique(tmpClust[which(fullForest_noNrml_importance[tmpRows,"MeanDecreaseGini"] >= keepThresh$MeanGini)])
						#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) }))) )
						mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( rownames(fullForest_noNrml_importance[tmpRows,])[which(tmpClust == x)] ) })))
					} else if(x =="adaBoost"){
						tmpClust <- NULL
						if(!is.element("adaBoost_clust",ls())){
							tmpClust <- findClusters(adaBoost_topImp$importance)
							assign("adaBoost_clust", tmpClust, pos=".GlobalEnv") # Save this as it's timely to compute....
						} else {
							tmpClust <- adaBoost_clust
						}
						tmpLabels <- unique(tmpClust[which(adaBoost_topImp$importance >= keepThresh$Importance)])
						#return( gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) }))) )
						mutPresence[[x]] <- gsub('GI.','GI:',unlist(lapply(tmpLabels, function(x){ return( names(adaBoost_topImp$importance)[which(tmpClust == x)] ) })))
					} else if(x=="polyAln"){
						tmpList <- lapply(expGenes,function(thisGene){ grep(thisGene,colnames(polyAln), value=TRUE) })
						names(tmpList) <- expGenes	
						#return( tmpList )
						mutPresence[[x]] <- tmpList
					} else if(x=="BayesTraits"){
						tmpList <- lapply(expGenes,function(thisGene){ 
									unlist(strsplit(grep(thisGene,unique(unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))])), value=TRUE),"::")) })
						names(tmpList) <- expGenes	
						#return( tmpList )
						mutPresence[[x]] <- tmpList
					} 
				}	
			}
			# Now we could visualise the overlap between my sets of analysis:  NOTE: I can use the doWeights =TRUE but it's not as visually informative.
				#plot(compute.Venn(Venn(list("tree"= unlist(mutPresence$tree),"randomForest"= unlist(mutPresence$randomForest),
				#							"adaBoost"= unlist(mutPresence$adaBoost), "polyAln"=unlist(mutPresence$polyAln), 
				#							"BayesTraits"=unlist(mutPresence$BayesTraits))), doWeights = FALSE,doEuler=TRUE))
			# NOTE: There are adaBoost and randomForest objects not in polyAln or BayesTraits since I subset these results for expGenes.
			
			# Here we define the overlap of siteRef between all my sets
			baseOverlaps <- findOverlap(mutPresence)
			#### NOW let's get a list of overlaps which show me all the BayesTraits pairs that have at least one site overlaps with correlation analysis sets.
			# First we'll find which rows of unique_sigPairs have some sites that exist with overlaps on BayesTraits analysis
			sigPair_overlaps <- findOverlap(list("BayesTraits"=unique(unlist(unique_sigPairs[,which(grepl("siteRef_",colnames(unique_sigPairs)))])),
													#"tree"=all_keepSites_mutTypes$Site,
													"randomForest"=mutPresence$randomForest,
													"adaBoost"=mutPresence$adaBoost),in_mainSet = "BayesTraits")
			sigPairs_rows <- unique(unlist(sapply(unique(unlist(sigPair_overlaps)),function(x){ union(which(unique_sigPairs[,"siteRef_1"] == x), which(unique_sigPairs[,"siteRef_2"] == x)) })))
			
			# From experience adaBoost is worthy as it is the most refined method of my having subsetted the value of these siteRef's so we will
			# consider all unique siteRef's which are correlated with some form of adaBoost analysis, 
			keep_single_siteRefs <- unique(unlist(sigPair_overlaps[which(grepl("adaBoost",names(sigPair_overlaps)))]))
			# Now we go and get the unique_sigPairs rows which have at least one of their siteRef values as a keep_single_siteRef
			overlap_sngSite <- unique_sigPairs[unlist(sapply(keep_single_siteRefs,function(x){ union(which(unique_sigPairs$siteRef_1 == x),which(unique_sigPairs$siteRef_2 == x)) })),]
									
			# This is the sites in our unique_sigPairs data.frame that we want to search for relation.
			siteRef_colSet <- c("siteRef_1","siteRef_2")
			# Now we extract the sigPairs which have dlbSupp but have at least one site supported by multiple phenotype correlation methods
			ovarlap_dblSite <- sapply(sigPair_overlaps[which(grepl("adaBoost",names(sigPair_overlaps)))],function(thisSet){ 
									# From within each set we will then return the unique_sigPairs row which corresponds to a dbl supported pair
									retFrame <- NULL
									for(thisRef in thisSet){
										# Now for each single siteRef value found in this thisSet set we will extract the rows
										# and then search if both elements in that row are elements of this thisSet
										for(thisCol in siteRef_colSet){
											tmpRows <- which(unique_sigPairs[,thisCol] == thisRef)
											if(length(tmpRows) > 0){
												for(thisRow in tmpRows){
													if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],thisSet)){
														retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
													}
												}
											}
										}
									}
									return( unique(retFrame) )  }, simplify = FALSE)
			
			# To find sigPairs when both elements are within my analysis sets we will start with the rows which contain single siteRef's and then see
			# if the other siteRef can be found
			indiv_phenoSets <- list(#"tree"=mutPresence$tree,
								"randomForest"=mutPresence$randomForest,"adaBoost"=mutPresence$adaBoost)
			
			# These are sigPairs which have at least one site supported by at least one phenotypic correlate method
			sigPairs_sngSupp <- sapply(indiv_phenoSets,function(thisSet){
									# From within each set we will then return the unique_sigPairs row which corresponds to a dbl supported pair
									retFrame <- NULL
									for(thisRef in thisSet){
										# Now for each single siteRef value found in this thisSet set we will extract the rows
										# and then search if both elements in that row are elements of this thisSet
										for(thisCol in siteRef_colSet){
											tmpRows <- which(unique_sigPairs[,thisCol] == thisRef)
											if(length(tmpRows) > 0){
												retFrame <- rbind(retFrame,unique_sigPairs[tmpRows,])
											}
										}
									}
									return( unique(retFrame) )  }, simplify = FALSE)
			
			
			# These are sigPairs which are both within the same individual phenotype correlation analysis
			sigPairs_dblSupp <- sapply(indiv_phenoSets,function(thisSet){ 
									# From within each set we will then return the unique_sigPairs row which corresponds to a dbl supported pair
									retFrame <- NULL
									for(thisRef in thisSet){
										# Now for each single siteRef value found in this thisSet set we will extract the rows
										# and then search if both elements in that row are elements of this thisSet
										for(thisCol in siteRef_colSet){
											tmpRows <- which(unique_sigPairs[,thisCol] == thisRef)
											if(length(tmpRows) > 0){
												for(thisRow in tmpRows){
													if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],thisSet)){
														retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
													}
												}
											}
										}
									}
									return( unique(retFrame) )  }, simplify = FALSE)
			
			
			# These are sigPairs where both sites are within at least one method but no consensus is required.  Here we should be passing a vector of all pheno correlates
			# the overlapping single sites, hence the difference to some of the above
			tmp_dblIndi <- sapply(unique(unlist(indiv_phenoSets)),function(thisSite){ 
									retFrame <- NULL
									# Now for each thisSite we extract the unique_sigPair rows that contain it by considering each siteRef_colSet 
									for(thisCol in siteRef_colSet){
										tmpRows <- which(unique_sigPairs[,thisCol] == thisSite)
										if(length(tmpRows) > 0){
											for(thisRow in tmpRows){
												# We look if the alternative site in this row is.element of all our individual phenotype correlations
												if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],unique(unlist(indiv_phenoSets)))){
													retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
												}
											}
										}
									}
									return( unique(retFrame) )  }, simplify=FALSE)
			sigPairs_dblIndi <- NULL
			for(thisSet in tmp_dblIndi){ sigPairs_dblIndi <- unique(rbind(sigPairs_dblIndi, thisSet)) }
			# These are sigPairs where at least one member agreed upon by both phenotype correlation methods, but both sites are within at least
			# one method.  Here we should be passing a vector of the overlapping single sites, hence the difference to some of the above
			tmp_sngOver_dblIndi <- sapply(sigPair_overlaps$BayesTraits_vs_randomForest_vs_adaBoost,function(thisSite){ 
									retFrame <- NULL
									# Now for each thisSite we extract the unique_sigPair rows that contain it by considering each siteRef_colSet 
									for(thisCol in siteRef_colSet){
										tmpRows <- which(unique_sigPairs[,thisCol] == thisSite)
										if(length(tmpRows) > 0){
											for(thisRow in tmpRows){
												# We look if the alternative site in this row is.element of all our individual phenotype correlations
												if(is.element(unique_sigPairs[thisRow,setdiff(siteRef_colSet,thisCol)],unique(unlist(indiv_phenoSets)))){
													retFrame <- rbind(retFrame,unique_sigPairs[thisRow,])
												}
											}
										}
									}
									return( unique(retFrame) )  },simplify=FALSE)
			sigPairs_sngOver_dblIndi <- NULL
			for(thisSet in tmp_sngOver_dblIndi){ sigPairs_sngOver_dblIndi <- unique(rbind(sigPairs_sngOver_dblIndi, thisSet)) }
			# This is just the list of all the pairs found in unique_sigPairs
			all_sigPairs = extractPair(unique_sigPairs, siteRef_colSet)
			
			# Now let's look for sitePairs that are defined when we consider having at least one site with multiple phenot correaltion support
			# and both sites being found within at least one form of phenotypic support, AND or both sites of a pair having overlapping support
			### NOTE: The value of the sngOverlap being compared is it ONLY considers those single sites which have at least adaBoost support
			suppPairs_list <- list("sngSupport"= unique(unlist(lapply(sigPairs_sngSupp, function(x){ extractPair(x, siteRef_colSet) }))),
									"sngOverlap"= extractPair(overlap_sngSite, siteRef_colSet),
									"dblOverlap"= unique(unlist(lapply(ovarlap_dblSite, function(x){ extractPair(x, siteRef_colSet) }))),
									"dblSupport" = unique(unlist(lapply(sigPairs_dblSupp, function(x){ extractPair(x, siteRef_colSet) }))),
									"dblIndi" = extractPair(sigPairs_dblIndi, siteRef_colSet),
									"sngOver_dblIndi" = extractPair(sigPairs_sngOver_dblIndi, siteRef_colSet))
			suppPairs_overlaps <- findOverlap(suppPairs_list)
			
			
			
			# Here we will use the protein protein interactions information to help visualise which of our protein interactions are supported.
			# These are the new pa01_ & pa14_genicDB and ppInt_actions & _links.detailed objects that I can use.... NOTE: for GI information I need
			# to work around a bit better as the GI of PA01 != GI: of PA14
			load("/Users/Jdench/Desktop/PARuns/ppInteractions/ppInt_actions_links.RData")
			load("/Users/Jdench/Desktop/PARuns/Sequences/PA01DBase/pa01_genicDB.RData")
			load("/Users/Jdench/Desktop/PARuns/Sequences/PA14DBase/pa14_genicDB.RData")
			# CHALLENGE IS THAT THE ppInt INFROMATION IS FOR THE REFERENCE STRAIN PA01 BUT ALL MY ANALYSIS HAS USED THE PA14 REFERENCE
			# THIS MEANS WE WILL HAVE TO USE BLAST TO TRY AND FIND THE SEQUENCE ANALOGS BETWEEN GI: SEQUENCES.
			######## BLAST definitions #############
			# This is the base directory for blast
			baseBlast <- paste("/Users/Jdench/Desktop/Blast/ncbi-blast-2.2.30+/",sep="")
			# This is the explicit location in which we will find the blast programs we might want to run
			blastLoc <- paste(baseBlast,"bin/",sep="")
			# This is the location of blast databases and where we will write our temporary ones
			blastDB <- paste(baseBlast,"db/Temp/",sep="")
			library(seqinr)
			
			# This is a function to perform a recursive intersection call on a list of information
			recursiveIntersect <- function(inList){
				if(length(inList) > 1){
					retIntersect <- inList[[1]]
					for(x in 2:length(inList)){
						retIntersect <- intersect(retIntersect,inList[[x]])
					}
					return( retIntersect )
				} else {
					# If there is only one elements than we return that list element as it intersects itself by default.
					return( inList )
				}
			}
			
			# We will extract those of our proposed pairs we want to explore and look if there are ppInt that exist
			# We will use all unique sngOverlap, bldOverlap and dblSupport sites
			ppInt_pairs <- unlist(lapply(unique(unlist(suppPairs_list)),function(thisPair){
								# We only really need the gene information not the actual site of the mutation.
								tmpPair <- unlist(lapply(strsplit(thisPair,"__:__")[[1]],function(x){ substr(x,1,regexpr("_",x)-1) }))
								# We now look if there are PA###### Locus_Tag's identifiable for all elements of thisPair
								tmpLoci <- unlist(sapply(tmpPair,function(x){ 
									# We want to handle if nothing is found by returning an NA value
									tmpIndex <- which(pa01_genicDB$geneID == x)
									if(length(tmpIndex) == 0){
										return( NA )
									} else {
										return( pa01_genicDB$Locus_Tag[which(pa01_genicDB$geneID == x)] )
									} })) 
								# Now we handle circumstance where nothing was found bu geneID reference by using a blast search of the sequence
								if(is.element(NA,tmpLoci)){
									# First we check if there is a pa01_genicDB database else we make one
									if (!file.exists(paste(blastDB,"pa01_genicDB.nsq",sep=""))){
										print(paste("Missing the blastdb for pa01_genicDB trying to create it",sep=" "))
										# We need to make certain there is a Seq_nt sequence otherwise we won't include that to our blastDB
										tmpIndex_seqNT <- which(nchar(pa01_genicDB$Seq_nt) > 0)
										write.fasta(as.list(pa01_genicDB$Seq_nt[tmpIndex_seqNT]),names=pa01_genicDB$Locus_Tag[tmpIndex_seqNT],file=paste(blastDB,"pa01_genicDB.fas",sep=""))
										system(paste(blastLoc,"makeblastdb -in ",blastDB,"pa01_genicDB.fas -dbtype nucl -input_type fasta -out ",blastDB,"pa01_genicDB",sep=""))
										unlink( paste(blastDB,"pa01_genicDB.fas",sep="") )
									}
									# Now we will blast the Seq_nt for this gene from pa14 to the pa01 db
									for(thisGene in tmpPair[which(is.na(tmpLoci))]){
										# We hand the ENTREZ values as some instances of my siteRef's may not have the same caps values, this shouldn't happen for anything other than Entrez: values.
										tmpSeq <- pa14_genicDB$Seq_nt[which(toupper(pa14_genicDB$geneID)== toupper(thisGene))]
										tmpFile <- paste(blastDB,"tmp_",gsub("[[:punct:]]","_",thisGene),".fas",sep="")
										# If we found a sequence then we will perform a blast search for this
										write.fasta(as.list(tmpSeq),names=thisGene,file=tmpFile)
										system(paste(blastLoc,"blastn -query ", tmpFile," -db ",blastDB,"pa01_genicDB -out ", blastDB,"tmpQuery_blast.table -outfmt '6 sseqid sstart send qstart qend qlen evalue pident sseq'",sep=""))
										# Now we look if we found anything by searching is a _blast.table file exists and if it has content
										if(file.exists(paste(blastDB,"tmpQuery_blast.table",sep="")) && length(readLines(paste(blastDB,"tmpQuery_blast.table",sep=""))) > 0){
											# We read in the blast.table information and run a bit of quality control by verifying the percIdent value and eValue
											tmpFrame <- read.table(paste(blastDB,"tmpQuery_blast.table",sep=""),header=FALSE, stringsAsFactors = FALSE)
											# These colnames are based on my knowledge of how I am making the blast call
											colnames(tmpFrame) <- c("refID","refStart","refEnd","qStart","qEnd","qLen","eValue","PercIdent","refSeq")
											if(tmpFrame$eValue <= 1e-40 || tmpFrame$PercIdent >= 90){
												tmpLoci[thisGene] <- tmpFrame$refID[union(which(tmpFrame$eValue == min(tmpFrame$eValue)),which(tmpFrame$PercIdent == max(tmpFrame$PercIdent)))[1]]
											}
											# We no longer need this tmpQuery_blast.table information
											unlink( paste(blastDB,"tmpQuery_blast.table",sep="") )
										}
										# we also no longer need the tmp_thisGene.fas file
										unlink( tmpFile )	
									}
								}
								# Now in theory we have a PA01 locus_tag value for each element of tmpPair so we look if there is a ppInt using the action
								# information since it seems that this is a subset of all possible ppInt, and only those with some activity support.
								# NOTE: In using the links.detailed information I could assume that a combined score being greater than 400 is meaningfull, from experience NOT all
								# possible protein-protein interactions exist in this object so we could even have performed a raw exists check (I did and got limited return values)
								# I chose not to use a cut-off for this combined score since my objective is to find possibly poorly characterised pairs of mutations, thus we should not prejudice our
								# ppInt found pairs via a requirement of strong supporting information
								
								#  We check that there are no more NA values in the tmpLoci
								if(!is.element(NA,tmpLoci)){
									#tmpList <- lapply(tmpLoci, function(x){ union(which(grepl(x,ppInt_actions$item_id_a)),
									#																				which(grepl(x, ppInt_actions$item_id_a))) })																										
									tmpList <- lapply(tmpLoci, function(x){ union(which(grepl(x,ppInt_links.detailed$protein1)),
																													which(grepl(x, ppInt_links.detailed$protein2))) })
									tmpIntersect <- recursiveIntersect(tmpList)																					
									# If we have found any tmp_intersectAction than we will return thisPair, else we do not
									#if(length(tmp_intersectAction) > 0){ 
									if(length(tmpIntersect) > 0){ 
										# This is a visualising sanity check 
										#print(thisPair)
										#print(ppInt_links.detailed[tmpIntersect,])
										# Return thisPair since we've found an intersect of protein protein interactions.	
										return( thisPair ) 
									}
								}
								# If we can't find PA##### values for both of thisPair's genes OR we can;t find a tmpIntersection then we return nothing, which defaults to NULL
								# Which itself will be removed when we unlist this lapply call
				}))
			
			
			######## END OF MY BASIC OVERLAP LISTS ALL BELOW MAY BE MORE SPECIALISED #########
			
			
			
			
			
			
			################################## FINAL PLOTS ##################################
			###### These are the final Venn plots that I will want to use for convincing Rees and maybe as SI in publication that my choices were sound.
			# These first two are meant to show that our methods of analysis provided numerous choices:
			# 1) The first is all possible BayesTraits sites and all Phenotype correlation sites as well as mutations in cannonical genes
			all_siteRef <- list("BayesTraits" = unique(unlist(unique_sigPairs[, siteRef_colSet])),
								"Phenotype" = unique(unlist(indiv_phenoSets)),
								"Connonical_Genes" = unique(unlist(mutPresence$polyAln)) )
			tmpVenn <- compute.Venn(Venn(all_siteRef), doWeights = FALSE,doEuler=TRUE, type = "circles")
			pdf(file=paste(workDir,"final_1_sitePairs_all_siteRefs_Venn.pdf",sep=""),height=8,width=8)
			plot(tmpVenn)
			dev.off()
			# 2) The second is to show that for simply considering sites that are correlated to resistance phenotype there were many options.
			tmpVenn <- compute.Venn(Venn(indiv_phenoSets), doWeights = TRUE,doEuler=TRUE, type = "circles")
			pdf(file=paste(workDir,"final_2_siteRefs_indiv_phenoSets_Venn.pdf",sep=""),height=8,width=8)
			plot(tmpVenn)
			dev.off()
			# The last two are to show how my choices were made, this will show that of the choices I made:
			# 3) There were some sites in expected genes but we were also able to define novel sites that could be correlated for either the phenotype or through evolution (compensatory then?)
				# This is for single siteRef #final3_list <- list("sng_phenoOverlap"=unique(unlist(strsplit(suppPairs_list$sngOverlap,"__:__"))), "expGene_polySite"=unlist(mutPresence$polyAln))
			final3_list <- list("sngSupp"= suppPairs_list$sngSupport,"sng_phenoOverlap"= unique(suppPairs_list$sngOverlap), "dbl_phenoSupport"= suppPairs_list$dblIndi, "expGene_polySite"= unique(unlist(lapply(expGenes,function(x){ all_sigPairs[which(grepl(x, all_sigPairs))] }))))
			tmpVenn <- compute.Venn(Venn(final3_list), doWeights = FALSE ,doEuler=TRUE, type = "squares")#, doWeights = TRUE,doEuler=TRUE, type = "circles")
			pdf(file=paste(workDir,"final_3_phenoOverlap_expPoly_Venn.pdf",sep=""),height=8,width=8)
			plot(tmpVenn)
			dev.off()
			# 4) This shows that we have a small refined set which are built from having single overlapping pheno support, and double sites pheno detection
			final4_list <- list("sngSupp"= suppPairs_list$sngSupport,"sng_phenoOverlap"= unique(suppPairs_list$sngOverlap), "dbl_phenoSupport"= suppPairs_list$dblIndi, "ppInt_support"= ppInt_pairs)
			tmpVenn <- compute.Venn(Venn(final4_list), doWeights = FALSE ,doEuler=TRUE, type = "squares")
			pdf(file=paste(workDir,"final_4_sngOverlap_dblSupport_Venn.pdf",sep=""),height=8,width=8)
			plot(tmpVenn)
			dev.off()
			# 5) This is our classic previously agreed upon visualiastion for our overlaps worth pursuing.
			final5_list <- list("sng_phenoOverlap"= unique(suppPairs_list$sngOverlap), "dbl_phenoSupport"= suppPairs_list$dblIndi, "ppInt_support"= ppInt_pairs)
			# If all the elements across the lists are redundant within one single largest list we will not plot weigthts, because we can't....	
			tmpVenn <- compute.Venn(Venn(final5_list), doEuler=TRUE, type = "circles", 
						doWeights = if (all(sapply(unique(unlist(final5_list)),function(x){ 
										is.element(x, final5_list[[which(length(final5_list) == max(length(final5_list)))]]) }))){ FALSE } else { TRUE } )
			pdf(file=paste(workDir,"final_5_sngOverlap_dblSupport_Venn.pdf",sep=""),height=8,width=8)
			plot(tmpVenn)
			dev.off()
			
			# Now we will finally be using can be defined by finding the overlap between ppInt and other support mechanism:
			keepPairs <- findOverlap(as.list(c(suppPairs_list,list("ppInt"= ppInt_pairs))), in_mainSet = "ppInt", in_otherSet = c("sngOverlap","dblIndi","sngSupport"))
			save(keepPairs, suppPairs_list, mutPresence, keepThresh, ppInt_pairs, thisGini, this_pVal, thisImportance,
				 file=paste(workDir,"keepPairs.RData",sep=""))
			## FROM REVIEW OF THE keepPairs list and final4_list Venn plot, I will consider valuable any of those sitePairs which have ppInt support,
			## and will consider building a complete empirical fitness landscape from the site pairs with overlap between all three refinement methods
			## and will include some "usual suspects" to simply explore the story a bit more thoroughly.  See my lab manual (ringed 3 subject Hillary notebook, on 28-03-2016)
			
			# This last plot is of the pairs that I will use and is just a means of visualising the alignment for the sites and site pairs being considered.
			load("/Users/Jdench/Desktop/PARuns/Kos2015/Kos2015_table.RData")	
			kos_phenoVec <- unlist(sapply(rownames(polyAln),function(x){ 
			if (is.element(x, c("PA14","PA01"))){ 
				return( "ND" )
			} else {
				as.character(bacFrame$Levo_res[which(as.character(bacFrame$Isolate) == x)])
			} } ))
			# We define the rows worth considering and setup some plotting information to be used downstream, we will keep the pairs that have ppInt support and double
			# phnoetypic support PLUS one site pair which we are proposing to be part of a complete empirical fitness landscape we will construct - "rpoB_864__:__gyrA_259"
			keepRows <- intersect(which(kos_phenoVec != "ND"),which(kos_phenoVec != "(I)"))
			plotSites <- list("sngSites"=unique(unlist(strsplit(keepPairs$ppInt_vs_sngOverlap, "__:__"))),
								"sitePairs"= unlist(strsplit(keepPairs$ppInt_vs_sngOverlap, "__:__")) )
			plotLength <- max(sapply(plotSites,length))
			numPlots <- 2
			ntNames <- c("A","C","G","T","-")
			phenoNames <- c("(R)","(S)")
			# I'm pallate of colours from a qualitative set that by knowledge accepts the length of my set (max = 9)
			# We set the phenoCols first as from visualising the set we know that blue and red are the first two choices
			allCols <- brewer.pal(n= length(c(phenoNames,ntNames)),name= "Set1")
			names(allCols) <- c(phenoNames,ntNames)
			
			# This is the phenotype alignments for the sites of interest which I'm considering now as the "sngOverlap_vs_dblSupport" group
			pdf(file=paste(workDir,"final_6_keepPairs_alignment.pdf",sep=""),height=32,width= numPlots * max(4,ceiling(plotLength/5)))
			par(mfcol = c(1, numPlots))
			for(thisPlot in plotSites){
				# This first plot is of the 
				plotMat <- cbind(kos_phenoVec[keepRows],toupper(polyAln[keepRows, thisPlot]))
				# Now we create a blank plot space and then will simply plot coloured pch squares based on the phenotype and nt types
				plot(x = NULL, y = NULL, xlim = c(1, plotLength + 1), ylim = c(1, nrow(plotMat)), ylab = "Phenotype and Nucleotide Alignment", xlab = "", axes = FALSE)
				# We use the colour sets and their names to be how we will search through this, these are not robust since the lists are manually generated.
				for(thisState in names(allCols)){
					tmpPlot <- which(plotMat == thisState,arr.ind = TRUE)
					points(x = tmpPlot[,"col"], y = tmpPlot[,"row"], pch = 15, cex = 1, col = allCols[thisState])
				}
				# We will add text to the bottom of the plot to indicate what information is being shown
				text(x = 0: length(thisPlot) + 1, y = 1, labels = c("Phenotype", thisPlot) ,pos = 1, offset = 2, srt = 90)
			}
			legend("bottom",cex=1 ,names(allCols),bty="n",fill= allCols, ncol = 4)
			dev.off()
			
			
			# We will also now include a plot of the network of correlated evolution that is proposed by our keepPairs$ppInt_vs_sngOverlap
			networkCols <- brewer.pal(3, name = "Set1")
			# We start by creating the data.frame of info we want to plot
			keepPairs_network <- data.frame("sigPair" = keepPairs$ppInt_vs_sngOverlap ,matrix(unlist(strsplit(keepPairs$ppInt_vs_sngOverlap, "__:__")), byrow = TRUE, ncol = 2, dimnames=list(NULL,c("siteRef_1","siteRef_2"))))
			# We add in the pValue_FDR for this sigPair, the location can be found by using the "all_sigPairs" object which is the same size as nrow(unique_sigPairs)
			keepPairs_network$pValue_FDR <- unique_sigPairs[vapply(keepPairs_network$sigPair, FUN.VALUE = vector(mode="numeric",length=1), USE.NAMES = FALSE,
																FUN = function(thisPair){ which(all_sigPairs == thisPair) }),"pValue_FDR"]
																
			# Now we adjust our pValue_FDR and "fix" any values which are 0, since adjacency matrix will not handle these normally
			# We will set all values of 0 to be 10% smaller than the smalles value, we also handle if the value of NA exists and we ignore it.
			if (length(intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))) > 0 && is.element(0,keepPairs_network$pValue_FDR) ){
				min_qVal <- min(keepPairs_network$pValue_FDR[intersect(which(keepPairs_network$pValue_FDR != 0),which(!is.na(keepPairs_network$pValue_FDR)))])
				keepPairs_network$pValue_FDR[which(keepPairs_network$pValue_FDR == 0)] <- as.numeric(min_qVal) * 0.1
			} else if (all(keepPairs_network$pValue_FDR == 0, na.rm = TRUE)) {
				# This means that our qValues are all equal to zero so we assign some bening small value
				keepPairs_network$pValue_FDR[which(!is.na(keepPairs_network$pValue_FDR))] <- 2e-16
			}
			# In order to use the iGraph adjacency matrix function with weighted edges we need to create this matrix
			# We will extract the siteRef information from rows which have qValues (so are not NA)
			adjVertices <- unique(as.vector(sapply(which(grepl("siteRef_",colnames(keepPairs_network))),function(j){
									sapply(which(!is.na(keepPairs_network$pValue_FDR)),function(i){ 
										return( keepPairs_network[i,j] ) 
									}) })))						
			adjMatrix <- matrix(0,ncol=length(adjVertices),nrow=length(adjVertices), dimnames=list(adjVertices, adjVertices))
			# Now we go through and extract qValue information from rows which are not NA (meaning this data had no instance)
			for (thisRow in which(!is.na(keepPairs_network$pValue_FDR))){
				tmpVal <- as.numeric(keepPairs_network$pValue_FDR[thisRow])
				adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]]] <- tmpVal
				adjMatrix[keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[2]], keepPairs_network[thisRow,which(grepl("siteRef_",colnames(keepPairs_network)))[1]]] <- tmpVal
			}
			aPlot <- graph.adjacency(adjMatrix, mode="undirected", weighted=TRUE)
			# In order to view these edges properly we will exxagerate the edge weights
			# We use the iGraph function of E()
			E(aPlot)$weight <- -log10(E(aPlot)$weight)
			pdf(file=paste(workDir,"final_7_keepPairs_network.pdf",sep=""),width=30,height=10)
			# This layout method of plotting is from Stephane's example of the iGraph graph.adjacency script see: "iGraph_network_Stephane.R"
			# Many of these plotting options are graphical parameters of igraph's plot.igraph
			
			# Here we create a vector to identify how we will colour our vertices, this uses information about if the gene holding a mutation
			# is one our our expected genes
			vertCols <- sapply(adjVertices,function(thisVertex){
				thisReturn <- NULL
				# Is the siteRef gene information one of the cannonical resistance genes?
				if (is.element(substr(thisVertex,1,gregexpr("_",thisVertex)[[1]][length(gregexpr("_",thisVertex)[[1]])]-1), expGenes)){
					return( networkCols[1] )	
				} else 
					return( networkCols[2] )
				})
					
			# Now we will handle the edge.widths by creating a binned set of sizes, we will make this set with 6 bins
			# Since this will create a vector of factors we use the factor level to represent the widths
			thickenLine <- function(inNums){ 2* seq(2,11,length.out = 10)[as.factor(inNums)] }
			edgeBin <- cut(E(aPlot)$weight, breaks=10); edge_plotThickness <- thickenLine(edgeBin) 
			# NOTE: this is a number of graphical parameters for the plot.igraph option as found through: "http://igraph.org/r/doc/plot.common.html"
			#plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= round(E(aPlot)$weight,3), vertex.size = 20, edge.label.cex=3, edge.width=edgeBin, 
			#	vertex.label.cex=3, vertex.color= vertCols, asp = 0.66, rescale = .5)
			plot(aPlot, layout=layout.fruchterman.reingold(aPlot),edge.label= NULL, vertex.size = 20, edge.label.cex=3, edge.width= edge_plotThickness , 
				vertex.label.cex=2.75, vertex.color= vertCols, asp = 0.66)
			
			legend("bottomleft",legend=levels(edgeBin),lwd= thickenLine(1:length(edgeBin)) ,title="Edge -log10(qValue)",bty="n",cex=3)
			legend("topleft",legend=c("Cannonical_Gene", "Novel_Gene"),fill= networkCols[1:2],title="Mutation Type",bty="n",cex=3)
			dev.off()
			
		# These should be closing out the original three looping elements for importance, gini and pVal	
		}
	}
}


q(save="no")