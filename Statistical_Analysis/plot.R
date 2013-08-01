###########################################################################
##
## QADE - Quality Assessment, Alignment, Differential Expression analysis
##
##			by Mark Onyango
##
## This file is part of the QADE Pipeline. For further information about 
## the functionality of the pipeline visit www.mybioinformatics.de
##
###########################################################################


## August 16th 2011:    - finalized plotDE.edger + bugfixes

## August 14th 2011:    - plotDE.edgeR added to plot edgeR DE data
##                      - plotDE.bayseq added to plot baySeq DE Data

## August 13th 2011:    - plotDispEsts not double logarithmic anymore



suppressPackageStartupMessages(library(hexbin))

baseVarFunc <- function(cds, cond){
	rvf <- rawVarFunc(cds,cond)
	sf <- sizeFactors(cds)[ conditions(cds) == cond ]
	xim <- sum(1/sf) / length(sf)
	function(q) rvf(q) + xim*q
}

xscale.components.log <- function( ... ){
	res = xscale.components.default( ... )
	res$bottom$labels$labels = do.call(expression,
		lapply(res$bottom$labels$at, function(a)
			substitute(10^b, list(b=a) ) ) )
	res
}

yscale.components.log <- function( ... ){
	res = yscale.components.default( ... )
	res$left$labels$labels = do.call(expression,
		lapply(res$left$labels$at, function(a)
			substitute(10^b, list(b=a) ) ) )
	res
}


plotdiag <- function(conditionobj,condition){
xg <- seq(log10(1/3),log10(max(conditionobj$baseMean)),length.out=1000)
print(xyplot(
	I( baseVar/baseMean^2 ) ~ baseMean,
	conditionobj[ conditionobj$baseMean>0,],
	panel=function( ... ) {
		panel.hexbinplot( ... )
		llines( xg, adjustScvForBias( baseVarFunc(cds,as.character(condition))(10^xg)/(10^xg)^2,
		attr(rawVarFunc(cds,as.character(condition))(NA), "size" )),col="orange",lwd=2)

		for(sf in sizeFactors(cds)[conditions(cds)==as.character(condition)])
			llines(xg, 1/(10^xg*sf),col="purple",lwd=1.5)
	},
	
	trans=function(x) x^(1/4),
	inv=function(x) x^4,
	xbins=80,
	scales=list(
		x=list(log=TRUE,axs="i"),
		y=list(log=FALSE,axs="i",tick.number=8 ) ),
	xlab="mean", ylab="squared coefficient of variation", main=paste("Diagnostics for condition",as.character(condition)),
	xscale.components = xscale.components.log,
) )
}


plotdiag2 <- function(conditionobj,condition){
xg <- seq(log10(1/3),log10(max(conditionobj$baseMean)),length.out=1000)
print(xyplot(
	I( baseVar/baseMean^2 ) ~ baseMean,
	conditionobj[ conditionobj$baseMean>0,],
	panel=function( ... ) {
		panel.hexbinplot( ... )
		llines( xg, adjustScvForBias( baseVarFunc(cds,as.character(condition))(10^xg)/(10^xg)^2,
		attr(rawVarFunc(cds,as.character(condition))(NA), "size" )),col="orange",lwd=2) # DESeq
	
		llines( xg, ( 10^xg + dge$common.dispersion * (10^xg)^2 ) / (10^xg)^2, col="orange",lwd=2,lty="dashed") # edgeR

		for(sf in sizeFactors(cds)[conditions(cds)==as.character(condition)])
			llines(xg, 1/(10^xg*sf),col="purple",lwd=1.5) # shot-noise  ->  Poisson
	},
	
	trans=function(x) x^(1/4),
	inv=function(x) x^4,
	xbins=80,
	scales=list(
		x=list(log=TRUE,axs="i"),
		y=list(log=FALSE,axs="i",tick.number=8 ) ),
	xlab="mean", ylab="squared coefficient of variation", main=paste("Diagnostics for condition",as.character(condition)),
	xscale.components = xscale.components.log,
) )
}


plotdiag3 <- function(conditionobj, condition){
xg <- seq(log10(1/3),log10(max(conditionobj$baseMean)),length.out=1000)
print(xyplot(
	baseVar ~ baseMean,
	conditionobj[conditionobj$baseMean > 0,],
	panel = function( ... ){
		panel.hexbinplot( ... )
		llines(xg, log10(baseVarFunc(cds,condition)(10^xg) ), col="orange", lwd=2) # DESeq
		llines(xg, log10(10^xg + dge$common.dispersion * (10^xg)^2 ), col="orange", lwd=2, lty="dashed") # edgeR

		for(sf in sizeFactors(cds)[conditions(cds) == as.character(condition)])
			llines(xg, log10(10^xg*sf), col="purple", lwd=1.5) # shot-noise  ->  Poisson
		},

	trans = function(x) x^(1/6),
	inv = function(x) x^6,
	xbins = 80,
	scales = list(
		x = list( log = TRUE, axs = "i"),
		y = list( log = TRUE, axs = "i", tick.number = 8) ),
	xlab = "mean", ylab = "variance",
	xscale.components = xscale.components.log,
	yscale.components = yscale.components.log
))
}


plotscv <- function(obj,cond){
	plot(obj$baseMean,(obj$baseVar/obj$baseMean^2),
		pch=20,cex=.5,log="x",
		xlab="BaseMean",
		ylab="Squared coefficient of variation (SCV)",
		main=paste("Diagnostics for condition",cond,sep=" ")
	)
	xg <- seq(log(1/3),log(max(obj$baseMean)), length.out=1000)
	lines(exp(xg), adjustScvForBias(baseVarFunc(cds,cond)(exp(xg))/(exp(xg))^2 , 
	attr(rawVarFunc(cds,cond)(NA), "size")),col="orange",lwd=2)	
}

plotDispEsts <- function( cds, cond )
{
	plot(
		rowMeans( counts( cds, normalized=TRUE ) ),
		cds@fitInfo[[cond]]$perGeneDispEsts,
		pch = '.', log = "x",
		xlab = "baseMean", ylab = "Dispersion",
		main = paste("Diagnostics for condition: ", cond, sep = "")
	)
	xg <- 10^seq( -.5, 5, length.out=300 )
	lines( xg, cds@fitInfo[[cond]]$dispFun( xg ), col="red" )
}

plotDE.deseq <- function( res, cond, pval = .05 )
{
	title <- sprintf( "DESeq: Significant (FDR: %.2f) FoldChange of conditions %s & %s", pval, cond[1], cond[2] )
	
    Adata <- res$baseMeanA
    Bdata <- res$baseMeanB
    
    Azeros <- which(Adata == 0)
    Bzeros <- which(Bdata == 0)
    
    infRatio <- ceiling(max(abs((log2(Bdata) - log2(Adata))[-union(Azeros, Bzeros)])))
    
    M <- log2(Bdata) - log2(Adata)
    
    M[Azeros] <- -infRatio - 2
    M[Bzeros] <- infRatio + 2
    
    A <- (Adata + Bdata)/2
    
    A[Azeros] <- Bdata[Azeros]
    A[Bzeros] <- Adata[Bzeros]
    
    plot(
        x = A,
        y = M, 
        ylim = c(-infRatio - 3, infRatio + 3), axes=FALSE,
		log = "x", pch = 20, cex = .3,
		col = ifelse( res$padj <= pval, "red", densCols(A) ),
		xlab = "baseMean", ylab = "log2-FoldChange",
		main = title
	)
    axis(side = 1)
    maxis <- pretty((-infRatio + 1):(infRatio - 1), min.n = 3, n = length(axTicks(side = 2)))
    maxis <- maxis[maxis < infRatio & maxis > -infRatio]
    axis(side = 2, at = c(-infRatio - 1, maxis, infRatio + 1), labels = c(-Inf, maxis, Inf))
    abline(h = c(-1, 1) * (1 + infRatio), col = "orange", lty = 3)
}

plotDE.edger <- function( obj, de.tags, cond, pval )
{
	if (!inherits(obj, what = "DGEList")) 
        stop("variable 'obj' must be of or descend from class 'DGEList'")
    
    title <- sprintf( "edgeR: Significant (FDR: %.2f) FoldChange of conditions %s & %s", pval, cond[1], cond[2] )
    
    group <- as.factor(obj$samples$group)
    if("pseudo.alt" %in% names(obj)){ # 01.08.2013 - syntax has obviously changed, adapting but keeping it backwards compatible
		x <- obj$pseudo.alt[,group == cond[1]]
		y <- obj$pseudo.alt[,group == cond[2]]
	}
	else{
			x <- obj$pseudo.counts[,group == cond[1]]
			y <- obj$pseudo.counts[,group == cond[2]]
	}


    Adata <- rowMeans(x)
    Bdata <- rowMeans(y)
    
    #Adata <- obj$conc$conc.group[,1]*obj$common.lib.size
    #Bdata <- obj$conc$conc.group[,2]*obj$common.lib.size
    
    Azeros <- which(Adata < 1)
    Bzeros <- which(Bdata < 1)
    
    infRatio <- ceiling(max(abs((log2(Bdata) - log2(Adata))[-union(Azeros, Bzeros)])))
    
    ## calculate the FoldChanges (log x - log y = log x/y)
    M <- log2(Bdata) - log2(Adata)
    
    ## set those ids with inf FoldChange to infRatio +/- 2
    M[Azeros] <- -infRatio - 2
    M[Bzeros] <- infRatio + 2
    
    A <- (Adata + Bdata)/2
    
    A <- sapply(A, function(x) if(x < 1){x = 0} else{x = x})
    
    plot(
    x = A,
    y = M,
    pch = 20, axes = FALSE, cex = .3,
    col = ifelse( names(A) %in% de.tags, "red", densCols(A) ),
    xlab = "log-baseMean", ylab = "log-FoldChange", log = "x",
    main = title
	)
    axis(side = 1)
    maxis <- pretty((-infRatio + 1):(infRatio - 1), min.n = 3, n = length(axTicks(side = 2)))
    maxis <- maxis[maxis < infRatio & maxis > -infRatio]
    axis(side = 2, at = c(-infRatio - 1, maxis, infRatio + 1), labels = c(-Inf, maxis, Inf))
    abline(h = c(-1, 1) * (1 + infRatio), col = "orange", lty = 3)
}

plotDE.bayseq <- function( cD, samplesA, samplesB, cond, pval = .05 )
{
    if (!inherits(cD, what = "countData")) 
        stop("variable 'cD' must be of or descend from class 'countData'")
    
    title <- sprintf( "baySeq: Significant (FDR: %.2f) FoldChange of conditions %s & %s", pval, cond[1], cond[2] )
    
    Adata <- cD@data[, samplesA]
    Bdata <- cD@data[, samplesB]
    
    selTags <- order(cD@posteriors[ , 2], decreasing = TRUE)
    topTags <- data.frame(cD@data[selTags, , drop = FALSE], FDR = cumsum(1 - exp(cD@posteriors[selTags, 2]))/1:nrow(cD@data))
    
    merged <- merge(cD@data, topTags, by = "row.names", sort = FALSE)
    
    ## normalize the count data
    Adata <- Adata/cD@libsizes[samplesA] * mean(cD@libsizes[c(samplesA, samplesB)])
    Bdata <- Bdata/cD@libsizes[samplesB] * mean(cD@libsizes[c(samplesA, samplesB)])
    
    if (nrow(cD@seglens) > 0){
        if (ncol(cD@seglens) == 1) {
            Adata <- Adata/cD@seglens[, 1]
            Bdata <- Bdata/cD@seglens[, 1]
        }
        else {
            Adata <- Adata/cD@seglens[, samplesA]
            Bdata <- Bdata/cD@seglens[, samplesB]
        }
    }
    
    ## combine the count data from the corresponding columns and calculate the median
    Adata <- colSums(t(cD@data[, samplesA]))/length(samplesA)
    Bdata <- colSums(t(cD@data[, samplesB]))/length(samplesB)
    
    ## filter ids that have zero count and save their row number in the following vector
    Azeros <- which(Adata == 0)
    Bzeros <- which(Bdata == 0)
    
    ## calculate the maximum log2 FoldChange that occurs and set anything above/below to +/-inf
    infRatio <- ceiling(max(abs((log2(Bdata) - log2(Adata))[-union(Azeros, Bzeros)])))
    
    ## calculate the FoldChanges (log x - log y = log x/y)
    M <- log2(Bdata) - log2(Adata)
    
    ## set those ids with inf FoldChange to infRatio +/- 2
    M[Azeros] <- -infRatio - 2
    M[Bzeros] <- infRatio + 2
    
    A <- (Adata + Bdata)/2
    
    A[Azeros] <- Bdata[Azeros]
    A[Bzeros] <- Adata[Bzeros]
    
    plot(
        y = M, 
        x = A, 
        ylim = c(-infRatio - 3, infRatio + 3), axes = FALSE, 
        ylab = "M", xlab = "A",log="x", main = title,
        pch = 20, cex = .3, 
        col = ifelse(merged$FDR < .05, "red", densCols(A))
    )
    axis(side = 1)
    maxis <- pretty((-infRatio + 1):(infRatio - 1), min.n = 3, n = length(axTicks(side = 2)))
    maxis <- maxis[maxis < infRatio & maxis > -infRatio]
    axis(side = 2, at = c(-infRatio - 1, maxis, infRatio + 1), labels = c(-Inf, maxis, Inf))
    abline(h = c(-1, 1) * (1 + infRatio), col = "orange", lty = 3)
    
    
}


########## Draw Pearson correlation scatterplots ##########

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.up <- function(x, y)
{
    points(x, y, pch = '.')
    abline(coef = c(0,1), col = "red")
}

correlate_samples <- function()
{
    output <- "sample_correlation"
    pdf(paste(output,paste(format(Sys.time(), "%Y-%m-%d"),".pdf",sep=""),sep="-"))
    par(bg="white")
    
    dta <- ct
    colnames(dta) <- conds
    
    pairs(dta, lower.panel = panel.cor, upper.panel = panel.up)
    
    dev.off()
}


bin2boxplot <- function(file){
    
    # This function will return TRUE if the value presented to it is an integer
    # and FALSE if not.
    check.integer <- function(N){
        !length(grep("[^[:digit:]]", format(N, scientific = FALSE)))
    }
    
    # Vectorizing the function check.integer enables us to check multiple numbers at once.
    vcheck.integer <- Vectorize(check.integer)
    
    # Read in the Boxplot data from the file
    t <- read.delim(file, header=FALSE, sep=" ")
    
    # Remove the last column as it contains only whitespaces
    t <- t[,-ncol(t)]
    
	# Create a matrix containing the cumulative sum up to each quality score at each cycle.
    # As the quality scores are compressed with the counting sum algorithm into an array,
    # the cumulative sum is needed to later identify the quality score. The score can be
    # identified as the cumsum can be used as an index for the compressed array.
    z <- apply(t,2,cumsum)
	
	# Get the total sum in each cycle as vector
    summe <- z[nrow(z),]
    
    bin2percentile <- function(p, xnz, notzero){
        j <- 1
        # If the value of the percentile is no integer the we use ceiling(percentile)
        # as the quasi-index in the cumsum vector.
        if(!check.integer(p)){
            
            # If the percentile is less than the first value in the non-zero cumsum vector
            # the location of the percentile is automatically known
            if(ceiling(p) < xnz[1]){
                # The value of the percentile's location is given in the notzero vector
                val <- notzero[1]
            }
            
            # As long as the location-value of the cumsum vector is less or equal to the percentile
            # iterate through the cumsum vector location values
            while(xnz[j] <= ceiling(p)){ 
                val <- notzero[j]
                j <- j+1
                
                # If the iteration has exceeded the length of the cumsum vector
                # obviously the last quality score corresponds to the percentile.
                if(j == length(xnz)+1){
                    j <- length(xnz)
                    cat("b")
                    break
                }
            }
            # Assign the value of the x-th percentile to the i-th cycle column of the result matrix.
            return(val)
        }
        else{			
            if(ceiling(p) < xnz[1]){
                val <- notzero[1]
            }
            while(xnz[j] <= p){ 
                if(xnz[j+1]-1 > xnz[j]){
                    val <- notzero[j]
                }
                else{
                    val <- (notzero[j]+notzero[j+1])/2
                }
                j <- j+1
            }
            return(val)
        }
    }
    
    # Determine the quantiles (vector)
    p10 <- summe * 0.10	
	p25 <- summe * 0.25
	p50 <- summe * 0.50
	p75 <- summe * 0.75
	p90 <- summe * 0.90
    
    # This matrix will contain the quantile values at each cycle
    result <- matrix(ncol=ncol(z),nrow=5)
    
	# Fill the result matrix one cycle at a time...
    for(i in 1:ncol(z)){
	    # Take the i-th cycle column from the z-matrix.
        x <- z[,i]
        
        # Which quality scores have zero counts int the t-matrix?
		notzero <- which(t[,i] != 0)
        
        # Subset the i-th column from the z-matrix.
        # Only non-zero elements are now contained.
        # Location vector as the values virtually correspond
        # to the location of the associated quality scores.
		xnz <- x[notzero]
		        
        result[1,i] <- bin2percentile(p10[i], xnz, notzero)
        result[2,i] <- bin2percentile(p25[i], xnz, notzero)
        result[3,i] <- bin2percentile(p50[i], xnz, notzero)
        result[4,i] <- bin2percentile(p75[i], xnz, notzero)
        result[5,i] <- bin2percentile(p90[i], xnz, notzero) 
	}

    
    xleft <- seq(0.5,ncol(result)-0.5,1)
    xright <- xleft+1
    ybottom <- result[2,]
    ytop <- result[4,]
    
    par(xpd=T, mar=par()$mar+c(0,0,0,2))
    filename <- unlist(strsplit(file, "_boxplotdata.txt"))
    plot(c(0,ncol(result)+1), c(0,50), type="n", xaxs="i", yaxs="i", xaxt="n",yaxt="n",
        xlab="Cycle", ylab="Quality")
    title(paste("Per cycle quality/error probability of file\n",filename, sep=""))
    axis(1, at = seq(0,ncol(result),ifelse(ncol(result) < 30,2,10)))
    axis(2, at = seq(0,50,5))
    axis(4, at = c(10,20,30,40), labels = c("10%","1%","0.1%","0.01%"))
    mtext("Error Probability",side=4,line=3)
    
    # Draw green, yellow and red rectangle
    rect(0.05, 30, ncol(result)+1, 49.95, col="lightgreen", border="transparent")
    rect(0.05, 20, ncol(result)+1, 30, col="yellow2", border="transparent")
    rect(0.05, 0, ncol(result)+1, 20, col="tomato1", border="transparent")
    
    # Draw actual boxes
    rect(xleft=xleft,xright=xright,ybottom=ybottom,ytop=ytop,col="azure2")
    
    # Draw Median
    segments(x0=xleft+0.1, x1=xright-0.1, y0=result[3,],y1=result[3,],col="red2",lwd=2)
    
    #Draw Whiskers
    segments(x0=xleft+0.5, x1=xleft+0.5, y0=result[1,],y1=result[2,],lty=2)
    segments(x0=xleft+0.5, x1=xleft+0.5, y0=result[5,],y1=result[4,],lty=2)
    segments(x0=xleft+0.2, x1=xright-0.2, y0=result[5,],y1=result[5,])
    segments(x0=xleft+0.2, x1=xright-0.2, y0=result[1,],y1=result[1,])
    
    par(mar=c(5, 4, 4, 2) + 0.1)
    
}
