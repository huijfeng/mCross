suppressMessages(library(MotIV, quietly=T))
suppressMessages(library(motifStack, quietly=T))
suppressMessages(library(getopt, quietly=T))
Sys.setenv(R_GSCMD="\"C:\\Program Files\\gs\\gs9.20\\bin\\gswin64.exe\"")
#####Database and Scores#####
## cluster motif using stamp github version output pairwise alignment scores
plotRPBmotiftree <- function(score.file, tf.file, plot=T, plot.format="pdf", dir){

  pcm.mat <- readPWMfile(tf.file)
  pw.score <- read.table(score.file, stringsAsFactors = F)

  pw.score <- as.dist(pw.score)
  creat_pcmobj <- function(x, name){
    motif <- new("pcm", mat=x, name=name)
    motif
  }

  pfm.mat <- list()
  for (i in 1:length(pcm.mat)){
    pfm.mat[[names(pcm.mat)[i]]] <- creat_pcmobj(pcm.mat[[i]] , names(pcm.mat)[i])
  }
  ## trim motif by Information content 0.2
  # pfm.mat <- lapply(pfm.mat, function(x) trimMotif(x, t=0.2))
  pfm.mat <- lapply(pfm.mat, pcm2pfm)

  # d <- motifDistances(lapply(pfm.mat, pfm2pwm), DBscores=jaspar.scores, cc="KL")
  ## using smu and cc as measurement
  d <- pw.score
  hc <- motifHclust(d, method="single")

  pfm.mat <- pfm.mat[hc$order]
  pfm.mat <- DNAmotifAlignment(pfm.mat, revcomp=rep(FALSE, length(pfm.mat)))

  if(plot){
    file.name <- sapply(strsplit(basename(tf.file),"\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))
    n <- length(pfm.mat)
    treewidth=1/6

    if(plot.format=="pdf"){
      width <- 7
      height <- (n/4)*7
      pdf(paste(dir, file.name, ".pdf", sep=""), width=width, height=height)
    }
    if(plot.format=="png"){
      width <- 480
      height <- (n/4)*480
      png(paste(dir, file.name, ".png", sep=""), width=width, height=height)
    }

    opar<-par(mar=c(2,0,0,0), mfrow=par("mfrow"))
    layout(matrix(c(rep(1,n),rep(2:(n+1),ceiling(1/treewidth)-1)),nrow=n,ncol=ceiling(1/treewidth)))
    plot(as.dendrogram(hc), main=NULL, horiz=T)
    # axis(side=1, at=seq(max(hc$height),0,by=-0.1))
    for (i in 1:(length(hc$height)-3)){
      abline(v=hc$height[i], lty=3, col=2)
    }
    par(mar=c(3.5,3.5,1.5,0.5))
    assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
    for (i in 1:(n - 1)) {
      plot(pfm.mat[[n - i + 1]], xlab = NA)
    }
    plot(pfm.mat[[1]], xlab=NA)
    rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
    par(opar)
    dev.off()
  }
}

plotRPBmotiftreeStamptree <- function(tree.file, tf.file, chidx.file, plot=T, plot.format="pdf", dir, treecut=0.05){

  pcm.mat <- readPWMfile(tf.file)

  tree <- read.tree(tree.file)

  chidx <- read.table(chidx.file, header=T)


  creat_pcmobj <- function(x, name){
    motif <- new("pcm", mat=x, name=name)
    motif
  }

  pfm.mat <- list()
  for (i in 1:length(pcm.mat)){
    pfm.mat[[names(pcm.mat)[i]]] <- creat_pcmobj(pcm.mat[[i]] , names(pcm.mat)[i])
  }
  ## trim motif by Information content 0.2
  # pfm.mat <- lapply(pfm.mat, function(x) trimMotif(x, t=0.2))
  pfm.mat <- lapply(pfm.mat, pcm2pfm)

  # d <- motifDistances(lapply(pfm.mat, pfm2pwm), DBscores=jaspar.scores, cc="KL")
  ## using smu and cc as measurement
  # d <- pw.score
  # hc <- motifHclust(d, method="single")
  tree.order <- match(tree$tip.label, names(pcm.mat))
  pfm.mat <- pfm.mat[tree.order]
  pfm.mat <- DNAmotifAlignment(pfm.mat, revcomp=rep(FALSE, length(pfm.mat)))

  if(plot){
    file.name <- sapply(strsplit(basename(tf.file),"\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))
    n <- length(pfm.mat)
    treewidth=1/4

    if(plot.format=="pdf"){
      width <- 7
      height <- (n/4)*7
      pdf(paste(dir, file.name, ".pdf", sep=""), width=width, height=height)
    }
    if(plot.format=="png"){
      width <- 480
      height <- (n/4)*480
      png(paste(dir, file.name, ".png", sep=""), width=width, height=height)
    }

    opar<-par(mar=c(2,0,0,0), mfrow=par("mfrow"))
    layout(matrix(c(rep(1,n),rep(2:(n+1),ceiling(1/treewidth)-1)),nrow=n,ncol=ceiling(1/treewidth)))

    plot(tree, show.tip.label=F)
    if(tree$root.edge!=0){
      axisPhylo()
      abline(v=getcuthc(tree, chidx, treecut)$hc, col="red", lty=3)
    }
    ##add tree cut lines need a function to decide n


    par(mar=c(3.5,3.5,1.5,0.5))
    assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
    for (i in 1:(n - 1)) {
      plot(pfm.mat[[n - i + 1]], xlab = NA)
    }
    plot(pfm.mat[[1]], xlab=NA)
    rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
    par(opar)
    dev.off()
  }
}


getcuthc <- function(tree, chidx, treecut){
  res <- list()
  if(dim(chidx)[1]<=2 || all(diff(chidx$C.H_Metric)>0) || all(diff(chidx$C.H_Metric)<0)){
    if(tree$root.edge>=treecut){
      res$hc <- tree$root.edge - treecut
      ## first class
      res$class <- 1
    }else{
      res$hc <- 0
      res$class <- 2
    }

  }else{
    k <- chidx$NumClust[which.min(chidx$C.H_Metric)]
    if(!is.null(k)){
      res$hc <- treek2hc(tree,k = k)$slice
    }
    res$class <- 3
  }
  res
}

treek2hc <- function(tree, k=NULL, h=NULL){
  tree.res <- list()
  # tree.hc <- as.hclust.phylo(tree)
  if(!is.null(k)){
    nodehc <- as.numeric(names(table(nodeHeights(tree))))
    num2hc <- list()
    num2hc[[1]] <- c(nodehc[1], nodehc[1])
    cutnum <- length(nodehc)
    if(cutnum==2){
      if(Ntip(tree)==2){
        num2hc[[2]]<- c(nodehc[1], nodehc[2])
      }else if(Ntip(tree)>2){
        i=2
        while(i>1 && i<=Ntip(tree)){
          num2hc[[i]]<- c(nodehc[2], nodehc[2])
          i <- i+1
        }
      }
    }else if(cutnum>=3){
      i=2
      while(i>1 && i < cutnum){
        num2hc[[i]] <- c(nodehc[i-1], nodehc[i])
        i <- i+1
      }
      j=cutnum
      while(j>=cutnum && j<=Ntip(tree)){
        num2hc[[j]] <- c(nodehc[cutnum], nodehc[cutnum])
        j <- j+1
      }

    }
    slice <- mean(num2hc[[k]])
    tree.res$slice <- slice
    # tree.res$group <- cutree(tree.hc, k)
  }else if(!is.null(h)){
    tree.res$slice <- h
    # tree.res$group <- cutree(tree.hc, h)
  }

  tree.res
}

creat_pcmobj <- function(x, name){
  motif <- new("pcm", mat=x, name=name)
  motif
}

summaryTree <- function(tree, score.file, hc, class, name){
  res <- list()
  res$RBPName <- name
  score <- read.table(score.file, header = T, sep="\t", stringsAsFactors = F)
  score <- score[!duplicated(score), ]
  mem.score <- list()
  mem.consensus <- list()

  for (i in 1:dim(score)[1]){
    mem.score[[score$name[i]]] <- score$score[i]
    if( grepl(",", score$consensus[i])){
      mem.consensus[[score$name[i]]] <- unlist(strsplit(score$consensus[i], split = ","))[1]
     }else{
      mem.consensus[[score$name[i]]] <- score$consensus[i]
    }
  }
  getrep <- function(x, mem.score){
   return(x$tip.label[which.max(as.numeric(unlist(mem.score[x$tip.label])))])
  }
  if(class==1 || class==3){
    clusters <- treeSlice(tree, slice=hc, trivial = T)
    res$ClusterNumber <- length(clusters)
    mem.list <- unlist(lapply(clusters, function(x) return(paste(x$tip.label, collapse = ","))))
    res$ClusterMember <- paste(mem.list, collapse = ";")
    rep.list <- unlist(lapply(clusters,function(x) getrep(x, mem.score) ))
    seq.list <- unlist(mem.consensus[rep.list])
    res$Representative <- paste(rep.list, collapse = ";")
    res$RepresentativeConsensus <- paste(seq.list, collapse = ";")
  }else if(class==2){
    res$ClusterNumber <- 1
    res$ClusterMember <- paste(tree$tip.label, collapse = ",")
    res$Representative <- tree$tip.label[which.max(as.numeric(unlist(mem.score[tree$tip.label])))]
    res$RepresentativeConsensus <- mem.consensus[[res$Representative]]
  }
  res$ClusterMethod <- class
  res
}

#argument mask: 0 - no argument, 1 - required argument, 2 - optional argument
optionSpec = matrix(c(
    'inputDir',    'i',    1, "character",
    'basename',   'n',    1, "character",
    'outDir',  'o',    1, "character",
    'plotFormat', 'p',  2, "character",
    'help'   ,   'h',    0, "logical"
    ), byrow=TRUE, ncol=4);

opt = getopt(optionSpec);

if ( !is.null(opt$help) |is.null(opt$inputDir) | is.null(opt$outDir)| is.null(opt$basename))
{
    #cat(getopt(optionSpec, usage=TRUE));
    cat (
        'plot motif clustering results from Stamp\n',
        'Usage: Rscript ', get_Rscript_filename(),"\n",
        '[required]\n',
        ' -i, --inputDir     Input folder\n',
        ' -n, --basename     Basename from Stamp\n',
        ' -o, --outDir       Plot file folder\n',
        '[options]\n',
        ' -p, --plotFormat   Plot format[pdf|png]\n',
        ' -h, --help         Print usage\n'
    );
    q(status=1);
}

if(!is.null(opt$inputDir) & !is.null(opt$basename)){
  transfac.file <- paste(opt$inputDir, "/", opt$basename,".transfac" sep="")
  tree.file <- paste(opt$inputDir, "/", opt$basename,".tree" sep="")
  chidx.file <- paste(opt$inputDir, "/", opt$basename, ".CHidx.txt", sep="")
}
if(!is.null(opt$outDir)){
  outdir <- opt$outDir
}
if(!is.null(opt$outDir)){
  plotformat <- opt$plotFormat
}else{
  plotformat <- "pdf"
}
plotRPBmotiftreeStamptree(tree.file, transfac.file, chidx.file,  plot.format = plotformat, dir=outdir)
