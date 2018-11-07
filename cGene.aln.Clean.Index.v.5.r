# This script used to review a suite of aligned genes which are part of a single genome, or that are desired to be concatenated
# It will identify any columns that represent insertions, create an index of a site's gene position, and concatenate the whole.


########## DEPENDENCIES ##########
# This is an updated version, and is dependent on cGene_aln.Blast.v.12.r
# This script requires use of the perl script catfasta2phyml.pl which must be stored manually in the folder of "muscleDir"
# This also now uses the new p14_genicDB information

# load library to read alignments
library(ape)
library(seqinr)
library(doMC)
library(foreach)

# First clean the session to make certain there are no other objects or items loaded
# This removes any hidden previous session info
system("rm .RData")
# This also removes any objects
rm(list=ls())

###################################################################################################
################# BEGIN OF PREDEFINITIONS #########################################################
###################################################################################################
farmHPCVL <- FALSE
######## Folder and working definitions ##########
baseDir <- if(farmHPCVL){
				"/home/hpc3058/jonathan/"
			} else {
				"/Users/Jdench/Desktop/"
			}
workDir <- paste(baseDir,"PARuns/Sequences/PA_wholeGenomes/",sep="")
# This is where I will put the the raw alignments of genes before any "cleaning"
muscleDir <- paste(workDir,"muscleAln/",sep="")
tmpDir <- paste(workDir,"tmpFiles/",sep="")
# This is a directory to which output should be written, it can be changed as required for intent BUT can be workDir...
outDir <- "/Users/Jdench/Desktop/PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8_retry/postAnal/compareLD/"
######## Genic Sequence defintions ###########
# This is the location of the file which contains the genomic database which holds our genes of interest
gDBdir <- paste(baseDir,"PARuns/Sequences/PA14DBase/",sep="")
# This is the location of the file which informs of the genes we want to search, and their sequences
# This is an RData file with the object pa14_genicDB holding all the genic sequences
load(paste(gDBdir,"pa14_genicDB.RData",sep=""))

# This will identify which genes this script will handle, the indexes should be changed as needed
workGenes <- as.character(pa14_genicDB$Locus_Tag)

######### Constants ############
registerDoMC(cores=4)
# These are a set of sequences we do not want to include in our final analysis or set
removeSeq <- c("PmendocinaNK01","PfluorescensSBW25","Pae86.genome","PresinovoransNBRC106553","Pae101.genome","PA_7","AZPAE14941")
# This is the sequence name for the sequence we define as reference.
refSeq <- "PA14"
# This is the suffix naming format used for our filenames that will relate to the Locus_Tag
fileSuffix <- "_seqAln.fas"
# This is the filename for the object we want to save
saveName <- "wholeAln_index.RData"
polyName <- "polyAln"
# This is the name to be used for our whole alignment concatenations
wholeName <- "wholeAln"
# This is the locaiton of the program called fasttree
fastTree_loc <- if(farmHPCVL){
				paste("/home/hpc3058/bin/",sep="")
			} else {
				paste(baseDir,"Jonathan_Dench/3rd_Party_Software/FastTree2.1.7/",sep="")
			}

# These are the characters we will accept as evidence that a position is well sequenced, not a gap, etc...
ntChars <- c("a","c","g","t")
# This will identify now many sequences or sites must be of sufficient quality for us to accept this gene/sequence
# as represented, core, worth holding onto....
min_numSeqs <- 0.8
min_numSites <- 0.75

###################################################################################################
################# END OF PREDEFINITIONS ###########################################################
###################################################################################################


###################################################################################################
########################### BEGIN OF FUNCTIONS ####################################################
###################################################################################################

# This is a function which assumes we have a data.frame reference database loaded (called pa14_genicDB)
# and will return the shape of the information to be collected by this script, defaults assume there was
# an error and will return information which will appear in that form, otherwise these arguments
# should be supplied to this call
returnFrame <- function(tmpName, tmpRow,tmpQuality="FALSE", tmpSize = NA, tmpSeq_pos= NA, tmpRef_pos = NA, tmpInsert_pos = NA){
	# Now we return the information we've compiled, default values are set so that we can track
	# where in quality analysis a file may have failed to pass
	return ( data.frame("Locus_Tag"=tmpName,"Gene"=pa14_genicDB$geneID[tmpRow],"Quality"=tmpQuality,"Size"= tmpSize, "Seq_ntPos"=tmpSeq_pos,
						"Ref_ntPos"=tmpRef_pos, "Insert_ntPos"=tmpInsert_pos, stringsAsFactors=FALSE) )
}

###################################################################################################
############################# END OF FUNCTIONS ####################################################
###################################################################################################



###################################################################################################
########################### BEGIN OF BODY #########################################################
###################################################################################################

# First we get a list of all the alignment files that we are analysing
allFiles <- list.files(muscleDir,pattern="_seqAln.fas")
allFiles <- allFiles[order(allFiles)]


if (!file.exists(paste(outDir,wholeName,".RData",sep="")) && !file.exists(paste(outDir,saveName,sep=""))){
	# Our goal is to create a data.frame that will include the Locus_Tag, Gene, reference internal nucleotide position for each element 
	# we will concatenate, further we will run a small bit of quality verification for reasons of sanity checking.
	wholeAln_index <- foreach(thisFile = allFiles, .combine="rbind") %dopar% {
		
	#wholeAln_index <- NULL
	#for (thisFile in allFiles){
		# We read thisFile 
		tmpAln <- read.dna(paste(muscleDir,thisFile,sep=""),format="fasta",as.character=TRUE)
		# This is the row in pa14_genicDB which relates to thisFile
		seqName <- sub(fileSuffix,"",thisFile)
		tmpRow <- which(pa14_genicDB$Locus_Tag == sub(fileSuffix,"",thisFile))
		
		#Sanity checking 
		#print(seqName)
		#print(tmpRow)
		
		# We remove any rows which have the rownames of removeSeq
		if (length(intersect(removeSeq,rownames(tmpAln))) > 0) { tmpAln <- tmpAln[-unlist(sapply(removeSeq,function(x){ which(rownames(tmpAln) == x) })),] }
		# We convert all geneChar characters into gap characters, and we reset it as a matrix
			#tmpAln <- matrix(mapply(function(x){ if(is.element(toupper(x),geneChar)){ return("-") } else { return(x) } }, tmpAln), nrow=nrow(tmpAln), dimnames=list(rownames(tmpAln),NULL))
		# Now we perform our first quality analysis and check if the reference strain sequence exists
		# within tmpAln, and then if the non-gap sequence matches the reference sequence in our database (pa14_genicDB)
		if (!is.element(refSeq,rownames(tmpAln))){ 
			return( returnFrame(tmpName = seqName, tmpRow = tmpRow) )
			#wholeAln_index <- rbind(wholeAln_index, returnFrame(tmpName = seqName, tmpRow = tmpRow))
		} else {
			# Our brute force measure of quality will be that the non-gap character remaining
			# portion of our refSeq column == the nucleotide sequence for the nucleotide sequence of pa14_genicDB
			tmpRef <- tmpAln[which(rownames(tmpAln) == refSeq),]
			# If there are gap characters we remove them to evaluate quality, and save those positions
			# to be considered as insertion points, NOTE: This assumes that there are not gap cahracters
			# defined in our reference sequence stored in the database
			if ( is.element("-",s2c(pa14_genicDB$Seq_nt[tmpRow])) ){
				# We return that there is an issue with thisFile because the reference has gaps
				# This exists to highlight issues and has no formal workaround as it has not yet been
				# found to be a problem
				return ( returnFrame(tmpName = seqName, tmpRow = tmpRow, tmpSize=0) )
				#wholeAln_index <- rbind(wholeAln_index, returnFrame(tmpName = seqName, tmpRow = tmpRow, tmpSize=0) )
			}
			insertPos <- 0
			if (is.element("-",tmpRef)) { 
				insertPos <-  which(tmpRef == "-")
				tmpRef <- tmpRef[-which(tmpRef == "-")]
			} 
			# We check if tmpRef is the same as our reference nucleotide strand (as per pa14_genicDB$Seq_nt[tmpRow])
			# If not we return a mostly null data.frame but with sufficient information to indicate where 
			# thisFile fell out of quality.
			if (paste(tmpRef,collapse="") != pa14_genicDB$Seq_nt[tmpRow]){
				return ( returnFrame(tmpName = seqName, tmpRow = tmpRow, tmpSize=0,tmpInsert_pos = insertPos) )
				#wholeAln_index <- rbind(wholeAln_index, returnFrame(tmpName = seqName, tmpRow = tmpRow, tmpSize=0,tmpInsert_pos = insertPos) )
			}
			# At this point the sequence has passed all our quality check points and we will then identify
			# all the information about tmpAln and relate if to our reference sequence
			
			# So we record the size of the file that will be concatenated 		
			thisSize <- ncol(tmpAln)
			# Now we need to identify, for each column of tmpAln, if it is a site in the reference
			# or if it is an insertion position
			tmpAln_ntPos <- matrix(NA,ncol=2,nrow=ncol(tmpAln))
			aCounter <- bCounter <- 1
			for (j in 1:ncol(tmpAln)){
				# We simply ask if the current position is within insertPos, if so we return
				# a tmpAln_ntPos in the second column (for insertion positions)
				if (is.element(j,insertPos)){
					tmpAln_ntPos[j,2] <- bCounter
					bCounter <- bCounter + 1
				} else {
					# we return a tmpAln_ntPos in the first column (for reference positions)
					tmpAln_ntPos[j,1] <- aCounter
					aCounter <- aCounter + 1
				}
			}
			# Now we simply get our data frame to return, submitting all parts we have calculated
			return( returnFrame(tmpName = seqName, tmpRow=tmpRow,tmpQuality="TRUE", tmpSize = thisSize, tmpSeq_pos = seq(1,ncol(tmpAln)), 
							tmpRef_pos = tmpAln_ntPos[,1], tmpInsert_pos = tmpAln_ntPos[,2]) )
			#wholeAln_index <- rbind(wholeAln_index, returnFrame(tmpName = seqName, tmpRow=tmpRow,tmpQuality="TRUE", tmpSize = thisSize, tmpSeq_pos = seq(1,ncol(tmpAln)), 
			#				tmpRef_pos = tmpAln_ntPos[,1], tmpInsert_pos = tmpAln_ntPos[,2]) )
		}
	}
	
	# Now we identify any bad genes that have been found, these can be reviewed for later and must be remembered to be excluded
	badRows <- which(wholeAln_index$Quality == FALSE)
	badGenes <- wholeAln_index[badRows,]
	# If there are any bad rows we remove them from the wholeAln_index object
	if (length(badRows) > 0){
		wholeAln_index <- wholeAln_index[-badRows,]
	}
	
	# Now we add the absolute nucleotide position so this can be related to our actual whole genome alignment
	wholeAln_index$Abs_ntPos <- seq(1,nrow(wholeAln_index))
	
	# Now we build the wholeAln file, using only the good genes, thus those identified/remaining in our wholeAln_index object 
	
	# Now we run the catfasta2phyml.pl script to create our concatenated alignment, however, if there badRows (i.e. genes) we need to remove 
	# those fasta files from our consideration.
	if (length(badRows) > 0){
		# We must go through each of the unique genes in wholeAln_index, we will temprorarily remove them from muscledir, then replace them after the catfasta2phyml.pl run
		tmpFiles <- unlist(sapply(unique(badGenes$Locus_Tag),function(x){ paste(x,fileSuffix,sep="")}))
		sapply(tmpFiles, function(x){ system(paste("mv ",muscleDir,x," ",tmpDir,sep="")) })
		# We run the concatenation, exactly as below
		system(paste("cd ",muscleDir," ; ./catfasta2phyml.pl -f *",fileSuffix," > ",outDir, wholeName,".fas",sep=""))
		# We return the tmpFiles to their original position
		sapply(tmpFiles, function(x){ system(paste("mv ",tmpDir,x," ", muscleDir,sep="")) })
	} else {
		# We can simply use all the gene and thus all the files in the muscleDir
		system(paste("cd ",muscleDir," ; ./catfasta2phyml.pl -f *",fileSuffix," > ",outDir, wholeName,".fas",sep=""))
	}
	
	# Now we load the whole alignment and then write out and RData object
	#wholeAln <- read.dna(paste(workDir,wholeName,".fas",sep=""),format="fasta",as.character=TRUE)
	# This work around had to be made to accomodate very large alignments as read.dna function did not handle long vectors
	tmpAln <- read.fasta(paste(outDir,wholeName,".fas",sep=""))
	# We identify all those positions in the alignment which are insert positions, to do this we look for NA in the Ref_ntPos column
	removeCols <- which(is.na(wholeAln_index$Ref_ntPos))
	wholeAln <- matrix(NA,ncol=length(tmpAln[[1]]),nrow=length(tmpAln))
	wholeAln_names <- names(tmpAln)
	for (thisRow in 1:nrow(wholeAln)){ wholeAln[thisRow,] <- tmpAln[[thisRow]] }
	rm( tmpAln )
	
	# We asses the quality of this alignment and our index by ensuring they are of the same size!
	alnToindex_Quality <- nrow(wholeAln_index) == ncol(wholeAln)
	# Sanity check placed in our Rout file
	print( paste("Our alignment and index data have good quality: ",alnToindex_Quality,sep="") )	

	# We now save and RData version of the whole alignment as well our index file and objects
	save(wholeAln,wholeAln_names, file=paste(outDir,wholeName,".RData",sep=""))
	save(wholeAln_index, badGenes, wholeName, alnToindex_Quality,removeSeq,file=paste(outDir,saveName,sep=""))
}
#########################################################################################
####################### Building Tree from the wholeName alignment ######################
#########################################################################################
# We build a tree for this alignment, but is from the wholeAlignment so I am 
#if (!file.exists(paste(workDir,wholeName,".tree",sep=""))){
	# While not necessary in all cases I do want trees made from these alignment so this will be done here and now:
#	system(paste(fastTree_loc,"fasttree -gtr -gamma -nt < ",workDir,wholeName,".fas > ",workDir,wholeName,".tree",sep=""))
#}
# We tidy some memory
rm( wholeAln,wholeAln_names, wholeName, alnToindex_Quality )

# We need the index file and badGenes file, so if it is not in the ls(), we will load it
if (length(which(grepl("wholeAln_index",ls()))) == 0){
	load(paste(outDir,saveName,sep=""))
}
# Now we review if there are any badGenes, this will have inherited the data frame class
if (nrow(badGenes) > 0){
	tmpFiles <- unlist(sapply(unique(badGenes$Locus_Tag),function(x){ paste(x,fileSuffix,sep="")}))
	# Then we must update the allFiles object to not consider these files/genes
	allFiles <- setdiff(allFiles,tmpFiles)
}

# We will load one file in order to get the nrow for our polyAln
tmpAln <- read.dna(paste(muscleDir,allFiles[1],sep=""),format="fasta",as.character=TRUE) 
# This will be the standardising order by which we will return the sequences, this is used for corrOrder to ensure that 
# our pieces are all put together for the proper sequences.

rowName_vec <- rownames(tmpAln)
# We remove any rows which have the rownames of removeSeq
if (length(intersect(removeSeq, rowName_vec)) > 0) { rowName_vec <- rowName_vec[-unlist(sapply(removeSeq,function(x){ which(rowName_vec == x) })) ] }

# We will not be keeping all the gene files in a reduced version of the wholeAln, this is because we want to remove all of those which
# do not carry a minimum sequence length within a minimum number of sequences.  We will store, for each gene represented by a file (in allFiles)
# if it is kept for the final reduced alingnment, and the reason why it was removed.  Lastly, we will also store what positions are polymorphic 
# in a growing matrix that we can use as our smallest alignment.
geneReview <- matrix(NA,nrow=length(allFiles),ncol=3,dimnames=list(NULL,c("Locus_Tag","Kept","wholeAln_polySites")))
#polyAln <- matrix(NA,nrow=nrow(tmpAln),ncol=1,dimnames=list(rownames(tmpAln),NULL)) 
# Since we will be cbinding many potrions to become polyAln we will store each portion seprately and then assemble them in one 
# action, this will save time writting the object to memoery once rather than after each iteration of the for loop.
polyPiece <- list()
aCounter <- 1

for (thisFile in allFiles){
	# We read thisFile 
	tmpAln <- read.dna(paste(muscleDir,thisFile,sep=""),format="fasta",as.character=TRUE)
	# This is the row in pa14_genicDB which relates to thisFile
	seqName <- sub(fileSuffix,"",thisFile)
	
	# We remove any rows which have the rownames of removeSeq
	if (length(intersect(removeSeq,rownames(tmpAln))) > 0) { tmpAln <- tmpAln[-unlist(sapply(removeSeq,function(x){ which(rownames(tmpAln) == x) })),] }
		
	
	# Now we look at the number of non-gap characters in each row
	qualPerc <- sapply(1:nrow(tmpAln),function(x){
		# We look i characters other than {"a","c","g","t"} are returned in a table of the characters
		tmpTable <- table(tolower(tmpAln[x,]))
		return( sum(tmpTable[which(is.element(dimnames(tmpTable)[[1]],ntChars))])/sum(tmpTable) )
	})
	
	# Now we identify if there are sufficient sequences which have sufficient non-"gap" positions
	if (length(which(qualPerc >= (1-min_numSites)))/nrow(tmpAln) <= min_numSeqs){
		# if so we're done, we enter that this gene is not kept
		geneReview[aCounter, ] <- c(seqName,"No","NA")
		aCounter <- aCounter + 1
	} else {
		# This means that we will keep this gene 
		
		# Now we will review all columns in tmpAln and identify if they are polymorphic
		polySite <- unlist(sapply(1:ncol(tmpAln),function(x){
			# We create a table of this column
			tmpTable <- table(tmpAln[,x])
			# Easiest is if the length of tmpTable == 1, this is a monomorphic site
			if (length(tmpTable) == 1){
				return( NULL )
			} else {
				# We look for any characters that are not of normal quality expectations, remove them from tmpTable
				tmpTable <- tmpTable[which(is.element(dimnames(tmpTable)[[1]],ntChars))]
				# Now this column may have had no elements of ntChars in which case this is now a blank
				# object and we could also stop going forward, ALSO if our tmpTable is now of length 1, then we also can stop (as per previous).
				if (length(tmpTable) <= 1){
					return( NULL )
				# There must be a larger fraction sequences with good characters than our min_numSeqs
				} else if (sum(tmpTable)/nrow(tmpAln) < 0.8){
					return( NULL )
				} else {
					# Otherwise our tmpTable has identified that this is a polymorphic site of quality positions
					return( x )
				}
			} }))
		# If there are no polymorphic sites we have nothing else to perform....
		if (length(polySite) > 0){
			# If there are we will keep all those cols which are polymorphic sites, and add them to our growing matrix
			# Before we do this we need to ensure that we have the correct order of sequences such that it matches polyAln
			corrOrder <- sapply(rowName_vec,function(x){ which(rownames(tmpAln) == x) })
			# Since we use our wholeAln_index to help recall wholeAln position equivalents we will draw out colnames from this
			tmpPos <- which(wholeAln_index$Locus_Tag == seqName)
			indexPos <- tmpPos[unlist(sapply(polySite,function(x){ which(wholeAln_index$Seq_ntPos[tmpPos] == x) }))]
			# Now is any of the Ref_ntPos are NA, this means we've returned some insert positions (HIGHLY UNLIKELY) but we will
			# record that they are insert positions of Seq_ntPos (so we know the position in the sequence file)
			tmp_colNames <- sapply(indexPos,function(x){ paste(wholeAln_index$Gene[x],
				# Here we look if this is an insert position
				if (is.na(wholeAln_index$Ref_ntPos[x])){
					paste("Insert_Seq_ntPos",wholeAln_index$Seq_ntPos[x],sep="")
				} else {
					wholeAln_index$Ref_ntPos[x]
				}, sep="_") })
			# Now we make the addition
			polyPiece[[length(polyPiece)+1]] <- matrix(tmpAln[corrOrder,polySite],nrow=nrow(tmpAln),ncol=length(polySite),dimnames=list(c(rownames(tmpAln)[corrOrder]),c(tmp_colNames)))
			geneReview[aCounter, ] <- c(seqName,"Yes",paste(polySite,collapse="_"))
			aCounter <- aCounter + 1
		} else {
			# There were no polymorphic sites
			geneReview[aCounter, ] <- c(seqName,"Yes","None")
			aCounter <- aCounter + 1
		}
	}
}
# Now we tidy up and remove the first column of polyAln
#if (is.na(polyAln[1,1])){ polyAln <- polyAln[,-1] }
# Now we assemble all the pieces found within polyPiece
polyAln <- foreach(thisPiece = polyPiece, .combine="cbind") %dopar% { return( thisPiece ) }

# Now we can save information we've generated 
save(geneReview,file= paste(outDir,"geneReview.RData",sep=""))
save(polyAln,removeSeq,polyName,file= paste(outDir,polyName,".RData",sep=""))
write.fasta(lapply(1:nrow(polyAln),function(x){ polyAln[x,] }),rownames(polyAln),paste(outDir,polyName,".fas",sep=""))

###################################################################################################
################################## END OF BODY ####################################################
###################################################################################################

q(save="no")

