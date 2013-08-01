cat("############################################################################\n")
cat("##                                                                        ##\n")
cat("## QUASI - Quality Assessment and Statistical Inference                   ##\n")
cat("##                                                                        ##\n")
cat("##                                                                        ##\n")
cat("## This file is part of the QUASI Pipeline. For further information about ##\n")
cat("## the functionality of the pipeline visit www.mybioinformatics.de        ##\n")
cat("##                                                                        ##\n")
cat("############################################################################\n")

cat("\n\nPipeline Version: 0.87\n\n")


# Defines a function to test wether the given package is installed or not
is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
}

# Check wether the required packages exist. If not, try to install them from
# the Bioconductor site.
packages <- c("edgeR", "DESeq", "baySeq", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "limma", "snow", "multicore", "tcltk2")
check <- is.installed(packages)

if(any(!check)){
    source("http://bioconductor.org/biocLite.R")
    cat("Installing missing packages (", packages[which(!check)] ,")...\n")
    biocLite(packages[which(!check)])
}

lapply(packages, function(x){ 
    cat("Loading ", x ," - Version: ", as.character(packageVersion(x)) ,"\n")
    suppressPackageStartupMessages(library(x, character.only=T, quietly = T, verbose = F)) 
    })

require(tcltk2)

# Presets the significance level to 5%
pval <- 0.05

# Booleans to check which packages have been used and which have not
bool_deseq <- FALSE
bool_edger <- FALSE
bool_bayseq <- FALSE

# To keep track of FastQ files, save them in this vector
fastqFiles <- c()

# Build the table containing the annotation information (e.g GeneBank ID, Gene ID, Gene Name) # 01.08.2013
if(missing(annotation_db)){
	cat("\n\nPre-building the database for the annotation of the result files later on...")
	total_ACCNUM <- rbind(toTable(org.Hs.egACCNUM), toTable(org.Mm.egACCNUM), toTable(org.Rn.egACCNUM))
	total_GENENAME <- rbind(toTable(org.Hs.egGENENAME), toTable(org.Mm.egGENENAME), toTable(org.Rn.egGENENAME))
	annotation_db <- merge(total_ACCNUM, total_GENENAME, by="gene_id", all.x=TRUE)
	cat("finished.\n")
}


configure <- function(){
    
    ## create DESeq countData object from the count table provided in matrix.txt
    fastqFolder <<- tclvalue(tkchooseDirectory(title="Choose folder of input fastq files:"))
        
    if(fastqFolder == ""){
        cat("You did not choose a folder! Please run configure() again.\n")
        return
    }
    
    fastqFiles <<- list.files(fastqFolder, pattern = ".faq$|.fq$|.fastq$", full.names=TRUE)
    
	# TODO : If no matrix.txt file has been found ask wether one already exists. If not try to create matrix.txt by running "count"
#    if(length(list.files(fastqFolder, pattern = "matrix")) != 0){
#        choice <- readline("Files containing the keyword 'matrix' have been found. Do you want to proceed with one of those? (y/n)")
#
#        switch(tolower(choice),
#               y={ file <- tclvalue(tkgetOpenFile(filetypes = "{{TXT} {.txt}} {{All files} *}", title="Select your count matrix file:")) },
#               n={ cat("Generating new count matrix file (matrix.txt).\n");
#                   code <- system(paste("count ", fastqFiles));
#                   if(code != 0) {cat("The count script encountered a problem! This usually implies inconsistencies in the FastQ file! If not, please send us a bug report.\n"); return();}
#                   file <- "matrix.txt" },
#                { cat("Input not recognized. Please call configure() again.\n"); return(); })
    
#        if(file == ""){ cat("You did not choose a file! Please run configure() again.\n"); return(); }
#    }
#    else{
#        code <- system(paste("count ", fastqFiles))
#        if(code != 0) {cat("The count script encountered a problem! This usually implies inconsistencies in the FastQ file! If not, please send us a bug report.\n"); return();}
#        file <- "matrix.txt"
#    }
	file <- tclvalue(tkgetOpenFile(filetypes = "{{TXT} {.txt}} {{All files} *}", title="Select your count matrix file:"))
    ct <<- read.delim(file, header=TRUE, stringsAsFactors=TRUE)
    rownames(ct) <<- ct$name
    ct <<- ct[,-1]
    
    ## save the total number of ids for later
    initial_number_ids <<- nrow(ct)    
    

    ## the user must input the experimental conditions for DESeq, baySeq and edgeR
	print("Header line of the count matrix file:")
	print(colnames(ct))
    conds <<- readline("Input conditions, in the same order as in the count table file, separated by spaces: ")
    conds <<- strsplit(conds," ")
    conds <<- unlist(conds)

    ## sanity check: is the number of entered conditions equal to the number of columns in the file
    if(length(conds) != ncol(ct))
        stop("The number of entered conditions does not match the number of columns in the count matrix file!")    
    
    
    ## correlating each sample with each other - only when there is more than one file
	if(length(fastqFiles) > 1)
		correlate_samples()    
    
    
    ## enter the conditions that are going to be used for pairwise differential expression analysis
    cat("Input conditions for differential expression analysis separated by spaces. Available conditions: \n")
    print(unique(conds))
    pairwise_list <<- character()
    pairwise_list <<- append(pairwise_list,readline("Enter conditions: "),length(pairwise_list))
	if(length(unique(conds)) > 2){
	    while(readline("Do you want to enter more conditions (y/n)? ") != "n"){
	        pairwise_list <<- append(pairwise_list,readline("Enter conditions: "),length(pairwise_list))
	        if(ncol(t(unlist(strsplit(pairwise_list[length(pairwise_list)]," ")))) > 2){
	            cat("Wrong entry! Please enter only two conditions separated by a single whitespace (e.g: A B)!\n")
	            cat("Removing last entry from list!\n")
	            pairwise_list <<- pairwise_list[-length(pairwise_list)]
	        }
	        if(ncol(t(unlist(strsplit(pairwise_list[length(pairwise_list)]," ")))) < 2){
	            cat("Wrong entry! Please enter two conditions separated by a single whitespace (e.g: A B)!\n")
	            cat("Removing last entry from list!\n")
	            pairwise_list <<- pairwise_list[-length(pairwise_list)]
	        }
	        if(!all(unlist(strsplit(pairwise_list," ")) %in% conds)){
	            cat("Wrong entry! Please enter conditions that exist in your experiment!\n")
	            cat("Removing last entry from list!\n")
	            pairwise_list <<- pairwise_list[-length(pairwise_list)]
	        }
	    }
	}
    
}
    
edgeR <- function(){
    
	bool_edger <<- TRUE    
	output <- "edger"
    pdf(paste(output,paste(format(Sys.time(), "%Y-%m-%d"),".pdf",sep=""),sep="-"))
    par(bg="white")
    
    ## create edgeR countData object
    cat("edgeR: Creating edgeR countData object\n")
    dge <<- new("DGEList")
    dge$counts <<- ct
    dge$samples <<- data.frame(group = conds, lib.size = colSums(ct))
    cat("edgeR: Calculating normalization factors\n")
    dge <<- calcNormFactors(dge)
    cat("edgeR: Calculating common dispersion\n")
    dge <<- estimateCommonDisp(dge)
    cat("edgeR: Calculating tagwise dispersion\n")
    dge <<- estimateTagwiseDisp(dge)
    
    ## performing edgeR pairwise differential expression analysis
    cat("edgeR: Testing for differential expression and writing output ( depending on your data this may take a while )\n")
    for(i in pairwise_list){
        res <- unlist(strsplit(i," "))
        string <- paste(paste("edger_",res[1],sep=""), res[2], sep = "")
        
        ## performing exact test to calculate p-values
        cat("edgeR: Performing Exact-Fisher test to calculate p-values (",res[1]," vs ",res[2],")\n")
        assign(string,exactTest(dge, pair = c(res[1],res[2])), envir = globalenv())
        
        assign("edger_com_libsize", dge$common.lib.size, envir = globalenv())
        
        ## how many ids have a p.value below set cutoff
        topX <- sum(p.adjust(get(string)$table$PValue,"BH") < pval) # 07.12.2011 - Use FDR instead of unadjusted p-values
                                                                    # 11.06.2012 - syntax has changed; $p.value -> $PValue
        
        ## save names of ids with p.value below cutoff
        de.tags <- rownames(topTags(get(string), n = topX)$table)
        
        plotDE.edger(dge, de.tags, c(res[1],res[2]), pval)
		
        
        ## generate output files
		cat("edgeR: Annotating the hitlist...\n")
		tmp <- topTags(get(string),n=nrow(get(string)$table))$table
        tmp <- cbind(id=rownames(tmp),tmp)
        tmp <- cbind(tmp,accession=substr(tmp[,1],as.integer(attr(regexpr('([0123456789]+)_',rownames(tmp),ignore.case=TRUE),"match.length"))+1,nchar(rownames(tmp))))
        tmp <- merge(tmp, annotation_db, by = "accession", all.x=TRUE)		
		
		assign(string, tmp, envir = globalenv()) # 01.08.2013 - FDR should obviously be saved with the rest of the results
        
        write.table(tmp,file=paste(string,"_all.txt",sep=""),sep="\t",row.names=FALSE)
        tmp <- tmp[tmp$FDR < .05,]	# 27.07.2012 - syntax has changed; adj.P.Val -> FDR
	    tmp <- tmp[which(!is.na(tmp$id)),]
        write.table(tmp,file=paste(string,".txt",sep=""),sep="\t",row.names=FALSE)
		
    }
    
    dev.off()
    
}


deseq <- function(){
    
	bool_deseq <<- TRUE    
	output <- "deseq"
    pdf(paste(output,paste(format(Sys.time(), "%Y-%m-%d"),".pdf",sep=""),sep="-"))
    par(bg="white")
    
    ## create an DESeq countData object and estimate sizeFactors and dispersion coefficients
    cat("DESeq: Creating DESeq countData object\n")
    cds <<- newCountDataSet(ct, conds)
    cat("DESeq: Calculating normalization factors\n")
    cds <<- estimateSizeFactors(cds)
    
    if(packageVersion("DESeq") < "1.5.1"){
    	cat("DESeq: Calculating variance functions\n")
    	cds <<- estimateVarianceFunctions(cds)
    
    	scvPlot(cds, ylim = c(0,2))
    
    	for(i in unique(conds)){
        	residualsEcdfPlot( cds, i )
    	}
    
    	for(i in unique(conds)){
        	diag <- varianceFitDiagnostics( cds, i )
        	if(!bool_edger)
                plotdiag(diag, i)
            else
                plotdiag2(diag, i)
    	}
    }
    else{
        cat("DESeq: Warning! You are using a DESeq Version more recent than 1.5.1\n")
        cat("                Skipping deprecated diagnostic plots.\n")
        cat("DESeq: Estimating dispersion coefficients\n")
        tryCatch(
            cds <<- estimateDispersions(cds,fitType="parametric",method="per-condition"),
            warning = function(x){
                cat("DESeq: Couldn't estimate dispersion using fitType = parametric. Switching to local fit!\n")
                cds <<- estimateDispersions(cds,fitType="local",method="per-condition")
            }
        )
        
        for(i in unique(conds)){
            plotDispEsts(cds, i)
        }
    }
    
    ## performing DESeq pairwise differential expression analysis
    for(i in pairwise_list){
        res <- unlist(strsplit(i," "))
        string <- paste(paste("deseq_",res[1],sep=""),res[2],sep="")
        
        cat("DESeq: Testing for differential expression of conditions",res[1],"vs",res[2],"and writing output ( depending on your data this may take a while )\n")
        
        ## performing exact test to calculate p-values
        assign(string,nbinomTest(cds,res[1],res[2]),envir = globalenv())
        
        plotDE.deseq( get(string), c(res[1],res[2]), pval )
        
        tmp <- get(string)[order(get(string)$pval),]
        
        cat("DESeq: Annotating the hitlist...\n")
        tmp <- cbind(tmp,accession=substr(tmp[,1],as.integer(attr(regexpr('([0123456789]+)_',tmp[,1],ignore.case=TRUE),"match.length"))+1,nchar(tmp[,1])))
		tmp <- merge(tmp, annotation_db, by = "accession", all.x=TRUE)
        
        assign( string, tmp, envir = globalenv());
        
        #tmp <- na.omit(tmp) ## 02.12.2011 Don't discard ids that have no annotation information
        write.table(tmp,file=paste(string,"_all.txt",sep=""),sep="\t",row.names=FALSE)		
        #tmp <- na.omit(tmp[tmp$padj < .05,]) ## 02.12.2011 Don't discard ids that have no annotation information
        tmp <- tmp[tmp$padj < .05,]
        tmp <- tmp[which(!is.na(tmp$id)),]
        write.table(tmp,file=paste(string,".txt",sep=""),sep="\t",row.names=FALSE)
    }
    
    if(packageVersion("DESeq") < "1.5.1"){
        cdsBlind <- estimateVarianceFunctions( cds, method = "blind" )
    }
    else{
        cdsBlind <- estimateDispersions( cds, method = "blind" )
    }
    
    vsd <- getVarianceStabilizedData( cdsBlind )
    dists <- dist(t(vsd))
    heatmap(as.matrix(dists),symm=TRUE)
    
    dev.off()
    
}

    
bayseq <- function(){

	bool_bayseq <<- TRUE    
	output <- "bayseq"
    pdf(paste(output,paste(format(Sys.time(), "%Y-%m-%d"),".pdf",sep=""),sep="-"))
    par(bg="white")
    
    ## sorting columns of count table and condition vector by condition
    ct <<- ct[,c(order(conds))]
    conds <<- conds[order(conds)]
    
    ## performing BaySeq pairwise differential expression analysis ( this takes way longer! )
    cat("BaySeq: Testing for differential expression and writing output ( this will take very long! )\n")
    cl <- makeCluster(multicore:::detectCores(), "SOCK")
    for(i in pairwise_list){
        res <- unlist(strsplit(i," "))
        string <- paste(paste("bayseq_",res[1],sep=""),res[2],sep="")

        ## choose pair of conditions given in res vector
        ct_bayseq <- counts(cds)[,conditions(cds) %in% res]
        conds_bayseq <- conds[conds %in% res]

        ## baySeq requires the count table to be of matrix type
        ct_bayseq <- as.matrix(ct_bayseq)
        
        ## how many times do the first and second condition occur
        first <- rle(conds_bayseq)[[1]][1]
        second <- rle(conds_bayseq)[[1]][2]
        
        ## create libsizes vector
        libsizes <- getLibsizes(data = ct_bayseq, estimationType="quantile")
        
        ## create replicates
        replicates <- c(rep(1,first),rep(2,second))
        
        ## create groups list
        groups <- list(NDE = c(rep(1,(first+second))),
                        DE = replicates)
        
        ## create baySeq countData object for further analysis
        CD <- new("countData", data = ct_bayseq, replicates = replicates, libsizes = libsizes, groups = groups)
        
        ## Negative Binomial Approach
        obj_cdp.nbml <- paste(string,"CDP.NBML",sep="_")
        obj_cdpost.nbml <- paste(string,"CDPost.NBML",sep="_")
        
        assign(obj_cdp.nbml, getPriors.NB(CD, samplesize = (0.9*initial_number_ids), estimation = "QL", cl = cl),envir=globalenv())
        #assign(obj_cdp.nbml, getPriors.NB(CD, samplesize = (1000), estimation = "QL", cl = NULL),envir=globalenv())
        assign(obj_cdpost.nbml, getLikelihoods.NB(get(obj_cdp.nbml), pET = 'BIC', cl = cl),envir=globalenv())
        
        plotDE.bayseq(get(obj_cdpost.nbml), samplesA = (1:first), samplesB = ((first+1):(first+second)), cond = c(res[1], res[2]), pval)
        
    }
    stopCluster(cl)
    
    dev.off()
    
}


hitOverlap <- function(){
    
	if(!bool_deseq | !bool_edger | !bool_bayseq){
		stop("You must first run all three packages to create the Venn diagram!")
	}    
	output <- "venn"
    pdf(paste(output,paste(format(Sys.time(), "%Y-%m-%d"),".pdf",sep=""),sep="-"))
    par(bg="white")
    
    for(i in pairwise_list){
    
        res <- unlist(strsplit(i," "))
        obj_bayseq <- paste(paste(paste("bayseq_",res[1],sep=""),res[2],sep=""),"_CDPost.NBML",sep="")
        obj_deseq <- paste(paste("deseq_",res[1],sep=""),res[2],sep="")
        obj_edger <- paste(paste("edger_",res[1],sep=""),res[2],sep="")
        
        ## retrieve ids with topCounts and order by id name
        bayseq <- topCounts(get(obj_bayseq), group = 2, number = initial_number_ids)
        bayseq <- bayseq[order(rownames(bayseq)),]
        deseq <- get(obj_deseq)[order(get(obj_deseq)$id),]
        #edger <- get(obj_edger)$table[order(rownames(get(obj_edger)$table)),]
        tmp <- topTags(get(obj_edger),n=nrow(get(obj_edger)$table))$table
        edger <- tmp[order(rownames(tmp)),]
        
        ## bind p.value, pval and FDR column together
        table_both <- cbind(bayseq=bayseq$FDR,deseq=deseq$padj,edger=edger$FDR)
        table_both <- as.data.frame(table_both)
        
        ## determine which id has p.value/FDR/pval less then 5%
        ovl_bayseq <- (table_both$bayseq < pval)
        ovl_deseq <- (table_both$deseq < pval)
        ovl_edger <- (table_both$edger < pval)
        
        ## calculate the overlap of the two hit lists and plot the resulting vennDiagram
        overlap <- cbind(ovl_bayseq,ovl_deseq,ovl_edger)
        a <- vennCounts(overlap)

	## produce overlap list with ids, found by all three packages
	rownames(overlap) <- deseq$id
	all_three <- as.data.frame(names(which(overlap[,1] & overlap[,2] & overlap[,3])),optional=T)
	colnames(all_three) <- "id"
	out <- merge(merge(merge(all_three,deseq,by="id"),edger,by.x="id",by.y="row.names"),bayseq,by.x="id",by.y="row.names")
	out <- out[,c(1,8,12,20,21)]
	colnames(out) <- c("ID","DESeq_FDR","edgeR_FDR","baySeq_Likelihood","baySeq_FDR")
	write.table(out,paste(paste(paste("overlap_",res[1],sep=""),res[2],sep=""),".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
        
        vennDiagram(
            a,
            names = c("BaySeq","DESeq","edgeR"),
            main = paste(paste("Significant IDs predicted by DESeq/BaySeq/edgeR of conditions",res[1],sep=" "),res[2],sep=" ")
        )
    }
    
    dev.off()
    
}


quality <- function( ... ){
	
	file.list <- list.files(fastqFolder, pattern = ".base_dist.txt|.boxplotdata.txt|.length_dist.txt|.phred_dist.txt")
    if(length(file.list) != 0){
        cat("Quality assessment files have been found. If you proceed those files will be overwritten!\n")
        choice <- readline("Proceed? (y/n): ")
        
        switch(tolower(choice),
               y={ cat("Proceeding with quality assessment...\n");
					file.list <- list.files(fastqFolder, pattern = ".faq$|.fq$|.fastq$", full.names=TRUE)
					check <- lapply(file.list, function(x) system(paste("qa ", x)))
    
					if(any(check != 0)){
						cat("The following files returned a non-zero value (this usually indicates an error):\n")
						cat(file.list[which(check != 0)],"\n")
					}
				},
               n={ cat("Not creating a new quality assessment.\n"); },
               { cat("Input not recognized. Please call quality() again.\n"); return(); })
    }
    
    output <- "quality-report"
    pdf(paste(output,paste(format(Sys.time(), "%Y-%m-%d"),".pdf",sep=""),sep="-"))
	par(bg="white")
	
	read_base_dist()
	read_length_dist()
	read_phred_dist()
	read_boxplotdata()

	dev.off()
}


read_base_dist <- function(){
	
	file.list <- list.files(path = "./", pattern = "_base_dist.txt", full.names = TRUE)
	
	if(length(file.list) == 0){
		cat("Warning: no _base_dist.txt files have been found!\n")
	}
	else{
		for(file in file.list){
			matrix <- read.delim(file, header = FALSE)
            filename <- unlist(strsplit(file,"_base_dist.txt"))
           	plot(x = 1:nrow(matrix), xlim = c(1,nrow(matrix)), ylim = c(0,100), type = "n", xlab = "Cycle", ylab = "Content [%]", main = paste("Base content per cycle of file\n", filename, sep = ""))
            lines(matrix[,1], col = "red")      # A
            lines(matrix[,2], col = "orange")   # T
            lines(matrix[,3], col = "green")    # G
            lines(matrix[,4], col = "blue")     # C
            lines(matrix[,5], col = 5)          # N
            #lines(matrix[,6], col = 6)          # GC
            grid()
            #legend("topright", legend = c("A","T","G","C","N","GC"), col = c("red","orange","green","blue",5,6), border = "white", lwd = 2)
            legend("topright", legend = c("A","T","G","C","N"), col = c("red","orange","green","blue",5), border = "white", lwd = 2)
		}
	}
}


read_length_dist <- function(){
	
	file.list <- list.files(path = "./", pattern = "_length_dist.txt", full.names = TRUE)

	if(length(file.list) == 0){
		cat("Warning: no _length_dist.txt files have been found!\n")
	}
	else{
		for(file in file.list){		
			distribution <- scan(file)
			
			filename <- unlist(strsplit(file, "_length_dist.txt"))
            #plot(distribution, pch = 20, cex = .5, xlab = "Length", ylab = "# of reads", main = paste("Sequence length distribution of file\n", filename, sep=""))
			
			plot(distribution/sum(distribution)*100, type="h", lwd=5, col=ifelse(distribution > 0 ,"red","white"), 
			xlab = "Length", ylab = "Percentage of reads", main = paste("Sequence length distribution of file\n", filename, sep=""))
			
			#lines(distribution, col = "red", lty = 2)
            grid()
		}
	}
}


read_phred_dist <- function(){
	
	file.list <- list.files(path = "./", pattern = "_phred_dist.txt", full.names = TRUE)

	if(length(file.list) == 0){
		cat("Warning: no _phred_dist.txt files have been found!\n")
	}
	else{
		for(file in file.list){
			distribution <- scan(file)
            
            mat <- matrix(distribution, ncol = 4)

			filename <- unlist(strsplit(file, "_phred_dist.txt"))
            plot(c(0,45),c(0,max(mat)+3), type="n", pch = 20, cex= .5, xlab = "Quality", ylab = "Percentage of bases", main = paste("Quality-per-base distribution of file\n", filename, sep=""))
			lines(x=0:45,mat[1:46,1], col = "red", lty = 2)
            lines(x=0:45,mat[1:46,2], col = "orange", lty = 2)
            lines(x=0:45,mat[1:46,3], col = "green", lty = 2)
            lines(x=0:45,mat[1:46,4], col = "blue", lty = 2)
            
            grid()
            
            legend("topleft", legend = c("A","T","G","C"), col = c("red","orange","green","blue"), border = "white", lwd = 1)
		}
	}
}


read_boxplotdata <- function(){

	file.list <- list.files(path = "./", pattern = "_boxplotdata.txt", full.names = TRUE)

	if(length(file.list) == 0){
		cat("Warning: no _boxplotdata.txt files have been found!\n")
	}
	else{
		for(file in file.list){
            bin2boxplot(file)
		}
	}
}
