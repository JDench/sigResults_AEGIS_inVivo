# After a run that uses BayesTraits to estimate pairs of sites under correlated evolution we have tended to a standard format for our output
# which we will collect using this script and then generate some post analysis from the collected information.

# This updated version is robust to if this was a single run or there were multiple runs held within a job directory.

########## DEPENDENCIES ##########
# This is an updated and more robust version of post analysis for output from Bayes Traits run's, 
# it comes after: PAFixedPostAn_foreachApproach.v.2.r
# This is dependent on output from cGene.aln.Clean.Index.v.1.r, MasterFixedSite.v.7.r, and MasterAlign_SearchNBinarise.v.14.r or equiv.

#### NOTE: This is not properly set-up to compiled equivalent data sets run under different parameters, it is when we've farmed out parts of a job into subsets.
####		Use CompareSets scripts, using the output of this script, in order to compare runs as described not being possible here.

# This removes any hidden previous session info
system("rm .RData")
# This also removes any objects
rm(list=ls())

# Load our libraries
library(ape)
library(seqinr)
library(doMC)
library(foreach)

# This causes all our data.frames to not use factor as the main type for strings loaded
options(stringsAsFactors=FALSE)

#######################################################################
########################### Initial Objects ###########################
#######################################################################

########### File and system location settings #############
# Basic working directory wherein all other information can be found
baseDir <- "/Users/Jdench/Desktop/"
# This is where our output data is found  ### CHANGE THIS ###
workDir <- paste(baseDir,"PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8_retry/",sep="")
# This is where we are holding a database of supporting information
dbDir <- paste(baseDir,"PARuns/Sequences/PA14DBase/",sep="")
# Results can be stored here to avoid clutter
resultDir <- paste(workDir,"postAnal/",sep="")

# This is the location of the blast databases and program 
blastLoc <- paste(baseDir,"Blast/ncbi-blast-2.2.30+/bin/",sep="")
blastDb <- paste(baseDir,dbDir,"blastDB_pa14/",sep="")


###### Constants ######
# Here we register cores for our foreach processes
registerDoMC(cores=4)
# Here is the name of the alignment file that should have been passed through out process  ### CHANGE THIS ###
alnName <- "polyAln"
# This is the name of the object our run should output and we want to load, and then the assembled object name
runOutput_object <- "pairAMat_adj_sorted"
compiled_object <- "compiledRun_data"
badCompile_object <- "badCompilations"
# This is our cutoff point for significance in terms of q-value
ourAccept <- 1e-5
# These are the extensions of output files that our run could have made
oFile <- c(paste(runOutput_object,".RData",sep=""),"BanalFixed.txt","BanalRun.txt")
# This is a value of the maximum number of repeated found siteID matches that exist we will investigate
maxSiteID <- 50
# This is a vector of the column names from our dBase object that we want to return in our gather_dbInfo function
gatherInfo_Genecols	<- c("Locus_Tag","Gene","Product","Type","Strand","COG", "Gene_Ontology")
gatherInfo_Sitecols	<- c("alignSite","Int_ntPos","Codon","Codon_ntPos")
# This is a logical to know if we've used the "fixedSite" approach (where we use sites of interest against the alignment)
###### CHANGE THESE ######
analysis_fixedSite <- TRUE
if(analysis_fixedSite){
	# And then further, the second indicates how we should find the int_ntPos and gene information.  It may use the posNames[alnPosition]
	# or geneName + ntPos inherited from the folder.  Otherwise I may have updated the actualy output to hold this information directly
	where_intPos <- "outputDirect"
}
				

########### Data Object Loading - Make changes if output names differ - ##########
# We look for the expected data objects that ought to have been created, the first being the polyFrame and reportFrame objects
if (length(list.files(path=workDir,pattern="_polyFrame.RData",recursive=FALSE)) == 1){
	load(list.files(path=workDir,pattern="_polyFrame.RData",recursive=FALSE,full.names=TRUE))
} else {
	print(paste(" Could not find a unique/existing polyFrame.RData object in ",workDir,sep=""))
}
# Now we look if there is an index file, otherwise we will expect that relevant information can be obtained from the alignment file
# Regardless we create an object called posNames which will hold the name of the position which maps to each traitMat location
useIndex <- FALSE
if (file.exists(paste(workDir,alnName,"_index.RData",sep=""))){
	load(paste(workDir,alnName,"_index.RData",sep=""))
	useIndex <- TRUE
	# Now we create posNames
	assign("posNames", unlist(sapply(1:nrow(eval(as.name(paste(alnName,"_index",sep="")))), function(x){ 
		paste(eval(as.name(paste(alnName,"_index",sep="")))$Gene[x], eval(as.name(paste(alnName,"_index",sep="")))$Ref_ntPos[x],sep="_") })),
		pos=".GlobalEnv")
} else {
	# This means our information should be held in the colnames of our alignment
	if (file.exists(paste(workDir,alnName,".RData",sep=""))){
		load(paste(workDir,alnName,".RData",sep=""))
		# We take posNames from the colnames of the alignment
		assign("posNames",colnames(eval(as.name(alnName))),pos=".GlobalEnv")
		# we now clear up some memory since we don't need the whole alignment
		rm (list=(alnName) )
	} else {
		print (" Could not find the alignment file and there is no _index file cannot proceed, please review ")
		q(save="no")
	}
}

# Now we look for the processed traitMat and other files from the MasterAlign_SearchNBinarise scripts, includes: traitMat, builtMat, patternMat, numMono, corrOrder
if (length(list.files(path=workDir,pattern=paste(alnName,"_traitMat_",sep=""),recursive=FALSE)) == 1){
	load(list.files(path=workDir,pattern=paste(alnName,"_traitMat_",sep=""),recursive=FALSE,full.names=TRUE))
} else {
	print(paste(" Could not find a unqiue/existing _traitMat_ file to load in ",workDir,sep=""))
}

# This is a reference database which contains additional information with respect to genes
load(file=paste(dbDir,"pa14_genicDB.RData",sep=""))

# We look if there is a resultDir else we will create it
if (!file.exists(resultDir)){ system(paste("mkdir -p ",resultDir,sep="")) }
											
#######################################################################
#######################################################################
#######################################################################

#######################################################################
############################### Functions #############################
#######################################################################
# This function will take a long string separated by any value and return a vector of each separated item
# It takes as input the string to be separated and the separator
deString <- function(sTring, sEp){
	# This will hold all the positions where the separator is found
	sepPos <-regexpr(sEp, sTring)
	# In the event there is no instance of sEp then sepPos should return as value -1, otherwise we need to substring
	if (sepPos != -1){
		# This is the vector that we will use to destring the input, we initalise it with the first character
		tmpVec <- substr(sTring, 1, sepPos[1]-1)
		
		# If sepPos has no more than 2 values we don't need to loop through it
		if (length(sepPos) > 2){
			for (s in 1:length(sepPos)-1){
				tmpVec <- c(tmpVec, substr(sTring, sepPos[s]+1, sepPos[s+1]-1))
			}
		}
		# Now we add the last item from the list 
		tmpVec <- substr(sTring, sepPos[length(sepPos)]+1, nchar(sTring))
	} else {
		# This is done to be certain that the format of the vector is consistent throughout
		tmpVec <- sTring
	}
	# Now we return our vector 	
	return(tmpVec)
}

# This will identify if all the objects in a vector fed in as v, are identical
all_identical <- function(v) all(sapply( as.list(v[-1]),FUN=function(z) {identical(z, v[1])}))


# This is a function to find rows where the contents of two columns appear together in alternate rows
find_duplicateRows <- function(inFrame, col1, col2){
	# We create a list of all the unions of col1 and col2 within inFrame
	allUnions <- lapply(1:nrow(inFrame), function(i){ return( union(inFrame[i,col1],inFrame[i,col2])[order(union(inFrame[i,col1],inFrame[i,col2]))] ) })
	return( unique(unlist(lapply(1:(nrow(inFrame)-1),function(thisRow){
		# we create an ordered union of thisRow's col1 and col2 elements
		tmpUnion <- union(inFrame[thisRow,col1],inFrame[thisRow,col2])[order(union(inFrame[thisRow,col1],inFrame[thisRow,col2]))]
		# Now we look for any of the ordered allUnions are identical to the ordered tmpUnion, this would mean it is an equivRow
		return( unique(unlist(lapply((thisRow +1):length(allUnions),function(thisUnion){ 
					if(all(allUnions[[thisUnion]] == tmpUnion)){ return( thisUnion ) } }))) )
	}))) )
}

# This is a function that will take a <gene>_<Ref_ntPos> to gather information from a reference database and return this
gather_dbInfo <- function(inString){
	# We first verify that a string with an "_" pattern has been sent through, otherwise we return NULL
	if (!grepl("_",inString)){ 
		return( NULL )
	} else {
		# Sanity check
		#print(inString)
		# we separate the gene and site information
		tmpGene <- substr(inString,1,regexpr("_",inString)-1)
		tmpGene <- data.frame(pa14_genicDB[which(pa14_genicDB$geneID == tmpGene), gatherInfo_Genecols])
		
		tmpSite <- as.numeric(substr(inString,regexpr("_",inString)+1, nchar(inString)))
		# Now we use the tmpSite information to inform on the nucleotide, codon and amino acid information, if tmpSite is NA this was an insert....
		# But realistically I would be surpirsed if this ever occured in our identified correlated positions (though not strictly impossible)
		if (is.na(tmpSite)){
			tmpSite <- data.frame(which(posNames == inString),"INSERT","INSERT","INSERT")
		# Now if this is a protein coding gene we will acquire some other information
		# Thus it must be tagged as gene and must be a multiple of 3 in length
		} else if (grepl("gene", tmpGene$Type) && nchar(pa14_genicDB$Seq_nt[which(pa14_genicDB$Locus_Tag == tmpGene$Locus_Tag)]) %% 3 == 0){
			# The codon can be found by simply dividing the internal nucleotide position (Int_ntPos) by 3 and rounding up
			# As for the codon's internal position we use modulo and recall that values of x %% 3 are the third positions
			tmpSite <- data.frame(which(posNames == inString),tmpSite,
				ceiling(tmpSite/3),
				if (tmpSite %% 3 == 0){	3 }else{ tmpSite %% 3 } )	
		# Otherwise it is intergenic, tRNA or RNA in which case the codon and aa information can be disregarded as NA
		} else {
			tmpSite <- data.frame(which(posNames == inString),"NA","NA","NA")
		}
		colnames(tmpSite) <- gatherInfo_Sitecols
		
		# Now we can return this information
		return( data.frame("siteRef"=inString,tmpSite,tmpGene) )
	}
}

##############################################################################################################################################
######################################################### BEGIN BODY OF SCRIPT ###############################################################
##############################################################################################################################################
# We first identify all the oFile[1] files which exist, as they indicate the sources of usable information 
allComplete_runFiles <- list.files(workDir,pattern=oFile[1],recursive=TRUE)
# This will just inform us on the failed runs by using our oFile type 2 and 3
allFailed_runFiles <- c(list.files(workDir,pattern=oFile[2],recursive=TRUE),list.files(workDir,pattern=oFile[3],recursive=TRUE))

if (!file.exists(paste(resultDir, compiled_object,".RData",sep=""))){
	# Now as our first step we load all the allComplete_runFiles, create a data frame of all their parts and run an FDR correction
	# NOTE: there may only be one completed run in which case this is strictly a waste of computational time so we check that
	if (length(allComplete_runFiles) > 1){
		assign(compiled_object, foreach(thisFile = allComplete_runFiles, .combine="rbind") %dopar% {
			# Bug checking for the foreach looping
			#debugFrame <- NULL
			#for (thisFile in allComplete_runFiles){		
	
			# We create a fake env in order to have a parent frame for loading this data, we replace any punctuation in this to be underscores
			# This step avoids any irregular characters or expressions, such as filepath slashes.
			envString <- gsub("[[:punct:]]","_",paste(thisFile,"_Env",sep=""))
			assign(envString, new.env(), pos = -1 ) #envir = parent.frame() )
			# Now we load thisFile and add all pieces together, BUT we need a new environemnt for each run to work in and load this into...
			load(paste(workDir,thisFile,sep=""), envir= eval(as.name(envString)) )
			# We review if the siteID_2 column of our output contains identical values of #2, if so that means we used a fixed site
			# and to get the positional information we must trust to the sub directory information held in thisFile: 
			# the name should be structured as:  <gene>/Site<alnName_position>/<runOutput_object>.RData
			if ( analysis_fixedSite ){
				if ( where_intPos == "alnPosition" ){
					# So we will update the infromation in this column to use the gene and alnName position information
					# regardless of useIndex we can use that Site information in posNames, we run a sanity check though to verify the gene matches
					tmpGene <- substr(thisFile,1,regexpr("/",thisFile)-1)
					tmpSite <- substr(thisFile,regexpr("/Site",thisFile)+5,gregexpr("/",thisFile)[[1]][2]-1)
					# Sanity check
					if ( tmpGene == substr(posNames[as.numeric(tmpSite)],1,nchar(tmpGene)) ){
						#rbind(debugFrame, cbind( eval(as.name(runOutput_object)),list("siteRef_2"=rep(posNames[as.numeric(tmpSite)],nrow(eval(as.name(runOutput_object))))) ) )
						return(	cbind( eval( as.name(runOutput_object), envir = eval(as.name(envString)) ),list("siteRef_2"=rep(posNames[as.numeric(tmpSite)],nrow(eval( as.name(runOutput_object), envir = eval(as.name(envString)) )))) ) )
					} else {
						print(" We had a problem confirming the position as tmpGene did not match the posNames gene ")
					}
				} else if( where_intPos == "folderDirect" ) {
					# So we will update the infromation in this column to use the gene and nt_intPos information inherited from the folder
					tmpGene <- substr(thisFile,1,regexpr("/",thisFile)-1)
					tmpSite <- substr(thisFile,regexpr("/Site",thisFile)+5,gregexpr("/",thisFile)[[1]][2]-1)
					return(	cbind( eval( as.name(runOutput_object), envir = eval(as.name(envString)) ),list("siteRef_2"=rep(paste(tmpGene,tmpSite,sep="_"),nrow(eval( as.name(runOutput_object), envir = eval(as.name(envString)) )))) ) )
				} else if( where_intPos == "outputDirect" ) {
					# We use the siteID_2 information directly, as this ought to be the information.
					return(	cbind( eval( as.name(runOutput_object), envir = eval(as.name(envString)) ),list("siteRef_2"=eval( as.name(runOutput_object), envir = eval(as.name(envString)))$siteID_2) ) )
				}
			} else {
				#rbind(debugFrame, cbind( eval(as.name(runOutput_object)),list("siteRef_2"=rep(posNames[as.numeric(tmpSite)],nrow(eval(as.name(runOutput_object))))) ) )
				return( eval( as.name(runOutput_object), envir = eval(as.name(envString)) ) ) 		
			}
			#} #Bug checking for the foreach loop, end of for statement	
		}, pos=".GlobalEnv")
		
		# if there was more than one output file read in then we want to adjust the fdr value and so will take care of that in this loop
		# This is to reset the p-value, since through saving we may have lost some of the information calculated by R
			# From experience this tmpPval is not required.... there are no rounding problems observed....
		tmpPval <- unlist(sapply(eval(as.name(compiled_object))$test_stat,function(x){ 1-pchisq(x,4) }))
		tmpPfdr <- p.adjust(tmpPval, "BH")
		# This will take the re-calculated p-value and get our FDR value using a Benjamani Hocheburg FDR correction
		assign(compiled_object, cbind(eval(as.name(compiled_object))[,setdiff(colnames(eval(as.name(compiled_object))),c("pvalue", "pvalue_fdr"))], data.frame("pValue"=tmpPval,"pValue_FDR"=tmpPfdr)), pos=".GlobalEnv")
	} else {
		# We simply load the one file and assign it to the same object after adjusting the column names for pvalue and pvalue_fdr
		tmpObject <- eval(as.name(load(paste(workDir, allComplete_runFiles,sep=""))))
		colnames(tmpObject)[which(colnames(tmpObject) == "pvalue")] <- "pValue"
		colnames(tmpObject)[which(colnames(tmpObject) == "pvalue_fdr")] <- "pValue_FDR"
		assign(compiled_object, tmpObject, pos=".GlobalEnv")
	}
	
	# Now we add siteRef_1 and _2 columns, but verify if we need to add the siteRef_2 column, else we only add the first
	# Note: the is.numeric call is to make certain this information can be interpreted as numeric, it ought to but if we don't get expected output at least this
	# should prevent a crash here, and we can inspect....
	if (!is.element("siteRef_2",colnames(eval(as.name(compiled_object))))){
		assign(compiled_object,
		cbind(eval(as.name(compiled_object)),
		list("siteRef_1"=sapply(eval(as.name(compiled_object))$siteID_1,function(x){ if (is.numeric(x)){ return( paste(posNames[siteID[[x]]],collapse="::") ) } else { return( x) } })),
		list("siteRef_2"=sapply(eval(as.name(compiled_object))$siteID_2,function(x){ if (is.numeric(x)){ return( paste(posNames[siteID[[x]]],collapse="::" ) ) } else { return( x) } })) ),
		pos=".GlobalEnv")
	} else {	
		assign(compiled_object,
		cbind(eval(as.name(compiled_object)),
		list("siteRef_1"=sapply(eval(as.name(compiled_object))$siteID_1,function(x){ if (is.numeric(x)){ return( paste(posNames[siteID[[x]]],collapse="::") ) } else { return( x) } })) ),
		pos=".GlobalEnv")
	}
	# Now we also add the length of the siteID values, this can be used later as a means of reviewing if a site of interest is generated as a result
	# of an ambiguous site pattern.  For siteRef_2, we use this rather than the siteID information since it is more reliably accurate and of a standard format
	
	# We create a vector, using all the unique siteRef_2 instances and then search for their length data, this is done separately since 
	# it is a time consuming step and we want to run this more efficiently
	tmp_lenSite2 <- vapply(unique(eval(as.name(compiled_object))$siteRef_2), USE.NAMES=TRUE,FUN.VALUE = vector(mode="numeric",length=1),
					function(x){
						polyPoint <- which(posNames == x)
						found = FALSE; tmp_refSite <- 1
						while(!found && tmp_refSite < length(siteID)){
							found = is.element(polyPoint,siteID[[tmp_refSite]])
							if(!found){ tmp_refSite <- tmp_refSite + 1 }
						}			
						if(found){ return( length(siteID[[tmp_refSite]]) ) } else { return( NA ) }
					})
	# Now we take our found siteID_2 lengths and assign then based on the names values we can associate.
	assign(compiled_object, cbind(eval(as.name(compiled_object)),
							list("len_siteID_1"=sapply(eval(as.name(compiled_object))$siteID_1, function(x){ length(siteID[[x]]) })), 
							list("len_siteID_2"=sapply(eval(as.name(compiled_object))$siteRef_2, function(x){ 
								return( tmp_lenSite2[which(names(tmp_lenSite2) == x)] ) })) ), 
							pos=".GlobalEnv") 						
	# Here we look for, and remove if present, any rows where the siteRef_ information has come back as ""
	# We've previously found this was occuring before implementing upstream removal of badGenes from our alignments, 
	# cGene.aln.Clean.Index.v.3.r or greater should not face this problem.
	badCheck <- unique(which(eval(as.name(compiled_object))$siteRef_1 == ""),which(eval(as.name(compiled_object))$siteRef_2 == ""))
	# We now assign a badCompile object 
	assign(badCompile_object, eval(as.name(compiled_object))[badCheck,], pos=".GlobalEnv")
	# And if there were any badCheck's we will update our compiled_object
	if (length(badCheck) > 0){  assign(compiled_object, eval(as.name(compiled_object))[-badCheck,], pos=".GlobalEnv")  }
	
	# Now let's save the compiled output
	save(list=(c(compiled_object, badCompile_object)),file=paste(resultDir, compiled_object,".RData",sep=""))
} else {
	load(paste(resultDir, compiled_object,".RData",sep=""))
}

# Now we extract only those parts that are deemed significant, and we use the pa14_genicDB reference object to start compiling information
# And we will then grab information from pa14_genicDB which informs on the sites found.

sigPairs <- eval(as.name(compiled_object))[which(eval(as.name(compiled_object))$pValue_FDR <= ourAccept), ]

# Now we define the unique significant pairs, first we get a list, then for each list entry we find which row has the smallest pValue_FDR
duplicateRows <- find_duplicateRows(sigPairs,"siteRef_1","siteRef_2")
removeRows <- if (length(duplicateRows) == 0){ 
					NULL 
				} else { 
					unlist(sapply(1:length(duplicateRows),function(x){
						# We look if this list entry is NULL, if not we need to review which to keep
						if (length(duplicateRows[[x]]) > 0){
							tmpSet <- c(x,duplicateRows[[x]])
							# now we return the non-minimum pValue_FDR rows as reported by sigPairs row x or those duplicates.
							return( tmpSet[-which(sigPairs$pValue_FDR[tmpSet] == min(sigPairs$pValue_FDR[tmpSet]))[1]] )
						} else {
							return( NULL )
						} })) 
				}
unique_sigPairs <- if (length(removeRows) > 0){ sigPairs[-removeRows,] } else { sigPairs }

# Now we may have multiple siteID's for either site, in which case we need to expand our returned information to accomodate all of these
# For thie reason we will create the sitesInfo as the ref genomic database such that each row relates to one siteRef element, we separate
# pairs but use the siteRef column as a reference to allow us to grab information based on the output held within sigPairs.

# we create a vector of all the unique siteRef information we will be working with.  We get the information from both siteRef_1 and _2
# Our first task is to identify those rows which have siteRef_ information from a site pattern relating to more than one site, this is found by a "::"
uniqueSites <- unique(unlist(unique(c(sapply(unique(sigPairs$siteRef_1),function(x){ strsplit(x,"::")[[1]] }),
						sapply(unique(sigPairs$siteRef_2),function(x){ strsplit(x,"::")[[1]] }) ))))
# We use the foreach mostly for the combine function but to save some time if there happends to be a large number of unique sites.
sitesInfo <- foreach(thisSite = uniqueSites, .combine="rbind") %dopar% {
	# We first address the siteRef_1 information, and use the "::" separator to identify the number of parts, if there is no instance this call will function
	return( gather_dbInfo(thisSite) )
}

# Now we reduce sitesInfo to unique siteRef instances
sitesInfo <- sitesInfo[unlist(sapply(unique(sitesInfo$siteRef),function(x){ grep(x,sitesInfo$siteRef)[1] })),]


#### This step is meant to check if something was written properly, but it is not properly coded, not to mention there is no intialisation step
save(sigPairs, unique_sigPairs, sitesInfo,file=paste(resultDir,"sitesInfo_",ourAccept,".RData",sep=""))

##############################################################################################################################################
######################################################### END BODY OF SCRIPT #################################################################
##############################################################################################################################################

q(save="no")
