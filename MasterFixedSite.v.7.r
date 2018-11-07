# This will be the master file that combines all portions previously written and considered
# Condensed as a single run file uses script previously used in:
# Author: Jonathan Dench - ControlFileWriter.v.3. ,eSiteGennInsert.v.1.5.,
# Author: Jean-Claude N., Stephane Aris-Brossou, Jonathan Dench - 1stTrials-MergedeSitesv.2.5.r

# Author for this file: Jonathan Dench, based from MasterNoBoot.v.5. this is prepped for a validation run
# vs. the previous pur. vs. pyr. trait matrix information.  Resultant from initial runs having poor recovery
# and specificity.

# This program is designed to detect epistasis at the SNP level from nucleotide sequences
# It considers an outgroup as being sequence(s) which have not experienced the selection
# of interest that may be creating epistasis, roots a tree with these and then removes them
# from BayesTraits analysis, not removing them would allow us to look for ALL epistasis that
# may be in effect within the sequences. (simply remove the drop.tip portions to perform this)
# NOTE if drop.tip is removed must update which objects are being written to outputs!

########## DEPENDENCIES ##########
# This is dependent directly upon MasterAlign_SearchNBinarise.v.14.r which is itself dependent upon:
# cGene.aln.Clean.Index.v.1.r and cGene_aln.Blast.v.12.r

# load library to read alignments
library(ape)
library(seqinr)
library(doMC)
library(foreach)

###################################################################################################
################# BEGIN OF PREDEFINITIONS #########################################################
###################################################################################################
# This section includes those objects/variables that the user should predefine before running the script

############# System Definitions ##############
# which computer is this script intended to be run on
wSystem <- 0 	# 1: run on SGE cluster, #2: run on HPCVL, any other value local machine
# define how many cores we will request with doPar
registerDoMC(cores=4)
# This is the root file path dependent on the location that a run is submitted.  Also wcAdj adjusts the wc value to report on data.
if (wSystem == 1){
    bDir <- "/export/home/jonathan/"
    phaseLoc <- "/export/home/jonathan/phase2.0/"
    bayesLoc <- "/export/apps/"
    wcAdj <- 1
}else if (wSystem == 2){
    bDir <- "/home/hpc3058/jonathan/"
    phaseLoc <- "/home/hpc3058/bin/phase2.0/"
    bayesLoc <- "/home/hpc3058/bin/BayesTraitsV2/"
    wcAdj <- 2
}else{
    bDir <- "/Users/Jdench/Desktop/"
    phaseLoc <- "/Users/Jdench/Desktop/Jonathan_Dench/3rd_Party_Software/phase-2.0-src/"
    bayesLoc <- "/Users/Jdench/Desktop/Jonathan_Dench/3rd_Party_Software/BayesTraitsV2-Beta-OSX/"
    wcAdj <- 1
}

# This is the working directory which should share the same filepath across all machines
wDir <- "PARuns/Sequences/PA_FixedSite/newTrials/polyAln_0.8_testTrial/"
# This sets the true working directory where a particular run will take place, NOTE this directory must exist
xDir <- "gyrA/"

# This is the actual working directory where we will find our current job
trueWD <- paste(bDir,wDir,sep="")
# These are the names of the tree files that will be used
treeFiles <- c("wholeAln_reduced.tree","wholeAln_reduced.nexus")

############# End System Definitions ##########

############# Preload Data ####################
# This loads our traitMat, builtID, simptMat, numMono, siteID and corrOrder information
# This is produced by MasterAlign_SearchNBinarise.v.9.r or later
load(paste(trueWD,"polyAln_traitMat_5perc_wholeAlnTree.RData",sep=""))
# This will load a vector of all the sites within this gene that this job script will handle
load(paste(bDir,wDir,xDir,"geneSites.RData",sep=""))
# This will load a list of all the site patterns for each fixed site we are observing, index of site must match list index
load(paste(bDir,wDir,xDir,"sitePatterns.RData",sep=""))

# Now we load the J. Dettman tree information we will later use
tTree <- read.nexus(paste(trueWD,treeFiles[2],sep=""))

### Sanity check for the corrOrder and our tree files ###

# Now we perform a sanity check, if the corrOrder does not cause the rownames of traitMat (or _names from an index file) to match the tTree$tip.label (species names)
# then we have had a problem in our data generation and we should not proceed further
if (length(rownames(traitMat)) == 0){
	# Then we look for an _index file in the trueWD, if it exists we load it and look for an _names object to have been loaded to ls()
	tmpFiles <- list.files(path = trueWD, pattern = "_index.RData")
	# If there is more one instance we cannot handle this properly without manual effort.
	if (length(tmpFiles) == 1){
		load(paste(trueWD,tmpFiles,sep=""))
		if (length(which(grepl("_names",ls()))) == 1){
			# Then we can use a unique _names object that ought to have just been loaded to act as the names in our sanity check
			if (mean(eval(as.name(ls()[grepl("_names",ls())]))[corrOrder] == tTree$tip.label) == 1){
				# This means that all the _names, in the correct order, match our tTree$tip.label
				print(" Sanity check complete, all is well with corrOrder, tTree$tip.label, and our alignments ")
			} else {
				print(" Sanity check complete with failure exiting R ")
				q(save="no")
			}
		} else {
			print(" Could not find a _names object in order to perform sanity check ")
		}
	} else if (length(tmpFiles) > 1){
		print(" Could not identify a unique _index.RData file to perform sanity check ")
	} else {
		print(" Could not find an _index.RData file to confirm a sanity check ")
	}
} else {
	# This means the rownames of traitMat can be used to compare to the species names in tTree$tip.label, corrOrder should be able to make these match
	if (mean(rownames(traitMat)[corrOrder] == tTree$tip.label) == 1){
		# This means that all the _names, in the correct order, match our tTree$tip.label
		print(" Sanity check complete, all is well with corrOrder, tTree$tip.label, and our alignments ")
	} else {
		print(" Sanity check complete with failure exiting R ")
		q(save="no")
	}
}
	


# Since we will never be performing site directed mutagenesis on multiple sites, we will simply exclude all siteID's
# for which there are more than 4 identical patterns, this should speed up calculations slightly
#validID <- which(sapply(siteID,function(x){ length(x)}) <= 4)  ## REMOVED FOR NOW AS SITES OF INTEREST MAY LAY HERE
	
# Now we loop this process for each site job in this gene.
for (thisSite in 1:length(geneSites)){
	curDir <- paste(bDir,wDir,xDir,"Site",geneSites[thisSite],"/",sep="")
	# If curDir already has a file "pairAMat_adj_sorted.RData" then we assume this run has been completed
	if (!file.exists(paste(curDir,"pairAMat_adj_sorted.RData",sep=""))){
		# We make a directory for this site's work, if it doesn't exist, if it does we clear it
		if(!file.exists(curDir))){ 
			system(paste("mkdir -p ",curDir,sep="")) 
		} else {
			# This means the run was started but not completed, we will restart it by deleting any txt intermediates that would be there.
			system(paste("cd ",curDir," ; rm *.txt",sep=""))
		}
		# This is done to unclass the read table, it ensures that what is returned is a simple vector.
		tTrait <- as.character(fixedList[[thisSite]])
		
		############# End Preload Data ################
		
		###################################################################################################
		################# END OF PREDEFINITIONS ###########################################################
		###################################################################################################
		
		
		###################################################################################################
		################# Begin Of Bayes Traits ###########################################################
		###################################################################################################	    
		
		# initialize object containing results
		pairAMat<- NULL
		# starts the clock
		ptm <- proc.time()
		# start the loops; only inner loop forks processes, since we farm our computations across multiple sites of interest we set the
		# i value to be the fixed site of interest, otherwise this should be a range of values in our siteID / builtMat objects loaded.
		for(i in paste(substr(xDir,1,nchar(xDir )-1),geneSites[thisSite],sep="_")){
		    pairAMat_tmp <- foreach(j = 1:length(siteID),.combine='rbind') %dopar% {
		        #print(paste("Now doing pair of sites ",j," / ",i,sep=""))
		        # Column 1 is going to be our fixed site
		        col1 <- tTrait
		        # Here we have traitMat sites as the second site against which we compare
		        col2 <- builtMat[,j]
		        
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
		            write.table(trait_data, paste(curDir,"traits.",i,".",j,".txt",sep=""),quote=F,row.names=F)
		            
		            # run BayesTraits under the dependent and indep models
		            ########### NOTE: The Input files InputBT_ML_indep.txt, InputBT_ML_dep.txt must be placed into each wd
		            system(paste(bayesLoc, "BayesTraits ",trueWD,treeFiles[2]," ",curDir,"traits.",i,".",j,".txt < ", trueWD,"InputBT_ML_indep.txt > ",curDir,"traits_res_indep.",i,".",j,".txt",sep=""))
		            system(paste(bayesLoc, "BayesTraits ",trueWD,treeFiles[2]," ",curDir,"traits.",i,".",j,".txt < ", trueWD,"InputBT_ML_dep.txt > ",curDir,"traits_res_dep.",i,".",j,".txt",sep=""))
		            
		            # extract likelihood values, by counting the # of lines in the traits_res_indep.i.j.txt file
		            system(paste("wc -l ",curDir,"traits_res_indep.",i,".",j,".txt > ",curDir,"wc",i,".",j,".txt",sep=""))
		            # This reads the wc file produce which should have simply the number of lines in the file
		            wc <- read.table(paste(curDir,"wc",i,".",j,".txt",sep=""))            
		            # Since the last line in the "traits" file is the valuable output from ByesTraits, we skip wc[1] -1 lines, to read the last.
		            res_indep <- read.table(paste(curDir,"traits_res_indep.",i,".",j,".txt",sep=""),header=FALSE,fill=TRUE,skip=as.numeric(wc[1]-wcAdj))
		            # Repeat now for the dependent model
		            system(paste("wc -l ",curDir,"traits_res_dep.",i,".",j,".txt > ",curDir,"wc",i,".",j,".txt",sep=""))
		            wc <- read.table(paste(curDir,"wc",i,".",j,".txt",sep=""))           
		            res_dep <- read.table(paste(curDir,"traits_res_dep.",i,".",j,".txt",sep=""),header=FALSE,fill=TRUE,skip=as.numeric(wc[1]-wcAdj))
		            # This adjustment is required only beacuse on the HPCVL server the traits_res_indep & _dep files contain an additional line
		            # after results section leading to a table read in with 2 rows, yet we only need to 1 row's value
		            
		            res_dep <- res_dep[1,]
		            res_indep <- res_indep[1,]
		            # cleanup tmp files releasing them from memory
		            unlink(paste(curDir,"traits_res_indep.",i,".",j,".txt",sep=""))
		            unlink(paste(curDir,"traits_res_dep.",i,".",j,".txt",sep=""))
		            unlink(paste(curDir,"traits.",i,".",j,".txt",sep=""))
		            unlink(paste(curDir,"traits.",i,".",j,".txt.log.txt",sep=""))
		            unlink(paste(curDir,"wc",i,".",j,".txt",sep=""))
		            unlink(paste(curDir,"traits.",i,".",j,".txt.log.txt.Schedule.txt",sep=""))
		            
		            # Since we expect our Bayes values to be negative this tests to assure us we are within expectations otherwise ignore
		            # Note these likelihood values are the log-likelihoods of the models.
		            if((res_indep[,2] < 0) & (res_dep[,2] < 0)){
		                # As this is a two sided test we multiply our test stat by two
		                test_stat <- 2 * abs(res_indep[,2] - res_dep[,2])
		                # this should be the last instruction of the foreach loop
		                # This caculatees the significance of our LRT and other information in a vector
		                return(data.frame("siteID_1"=j,
		                					"siteID_2"=i,
		                					"lkl_indep"=res_indep[,2],
		                					"rootP_00_indep"=res_indep[,ncol(res_indep)-3],
		                					"rootP_01_indep"=res_indep[,ncol(res_indep)-2],
		                					"rootP_10_indep"=res_indep[,ncol(res_indep)-1],
		                					"rootP_11_indep"=res_indep[,ncol(res_indep)],
		                					"lkl_dep"=res_dep[,2],
		                					"rootP_00_dep"=res_dep[,ncol(res_dep)-3],
		                					"rootP_01_dep"=res_dep[,ncol(res_dep)-2],
		                					"rootP_10_dep"=res_dep[,ncol(res_dep)-1],
		                					"rootP_11_dep"=res_dep[,ncol(res_dep)],
		                					"test_stat"= test_stat,
		                					"pvalue"=1-pchisq(test_stat,4), stringsAsFactors=FALSE) )
		            }
		        }
		    }
		    # merge objects into pairAMatto build a memory of each bayesTraits run
		    pairAMat<- pairAMat_tmp
		}
		# Stop the clock
		time_par1 <- proc.time() - ptm
		pairAMat<- na.omit(pairAMat)
		pairAMat<- data.frame(pairAMat)
		
		# FDR calculation, add actual site positions and sort by P-value
		adj_p <- p.adjust(pairAMat$pvalue, "BH")
		pairAMat_adj <- cbind(pairAMat, adj_p)
		
		colnames(pairAMat_adj) <- c(colnames(pairAMat), "pvalue_fdr")
		pairAMat_adj_sorted <- pairAMat_adj[order(pairAMat_adj$pvalue),]
		
		# write results to file
		#write.csv(pairAMat_adj_sorted, file=paste(curDir,"pairAMat_adj_sorted.csv",sep=""), quote=F, row.names=F) # relocated and set row.names to "F"
		save(pairAMat_adj_sorted,file=paste(curDir,"pairAMat_adj_sorted.RData",sep=""))
		
	} # This closes out the conditional of doing this run if the output file hasn't already been created.
}

q(save="no")
