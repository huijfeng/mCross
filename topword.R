library (limma);

#get the list of top word showing asymmetric enrichment

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2)
{
    n <- length (args);
    print (args);
    stop ("topword.R <in.w7.zscore.txt> <out.dir>\n");
}

infile <- args[1];
outdir <- args[2];

RBP <- read.table(infile, sep="\t", header=T);

RBP.data <- as.matrix(RBP[,-1]);
RBP.names <- colnames(RBP)[-1];

colnames(RBP.data) <- RBP.names;
rownames(RBP.data) <- as.character(RBP[,1]);

dir.create(outdir);

word.sig <- matrix(data=NA, nrow=length(RBP.names), ncol=3);
colnames(word.sig) <- c("q_0.01", "q_0.05", "q_0.1");

for (i in 1:length(RBP.names))
{
	cat(RBP.names[i], "...\n");
	z <- RBP.data[,i];
	names(z) <- RBP[,1];

	z.sort <- sort(z, decreasing=T);
	ntry <- length(which(z.sort>0));
	q <- array (data=NA, dim=ntry);
	names(q) <- names(z.sort)[1:ntry];
	
	for (j in 1:ntry)
	{
		threshold <- z.sort[j];
		k <- length(which (z.sort <= -threshold));
		q[j] <- min(k/j, 1);
	}
	
	#make q monotonic
	for (j in ntry:1)
	{
		q[j-1] <- min(q[j-1],q[j]);
	}
	
	
	word.sig[i,1] <- length(q[q<0.01]);
	word.sig[i,2] <- length(q[q<0.05]);
	word.sig[i,3] <- length(q[q<0.1]);
	
	q.top <- q[q<0.1];
	write.table(q.top, paste(outdir, "/top.", RBP.names[i], ".txt", sep=""), quote=F, col.names=F);
}

write.table(cbind(RBP.names, word.sig), paste(outdir, "/w7.sig.summary.txt", sep=""), sep="\t", row.names=F, quote=F);



