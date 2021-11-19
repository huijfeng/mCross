library (gplots);

#plot the frequency of top words around crosslink sites

args <- commandArgs(trailingOnly=TRUE)

if (length(args) <3)
{
    n <- length (args);
    print (args);
    stop ("plot_topword_xl.R <in.profile.txt> <out.profile.pdf> <RBP.name>\n");
}

infile <- args[1];
outfile <- args[2];
RBP.name <- args[3];
w=7;
h=7;
if(length(args) > 3)
{
	w=as.integer(args[4]);
}
if(length(args) > 4)
{
	h = as.integer(args[5]);
}


profile <- read.table (infile, sep="\t", header=T);
profile.data <- as.matrix(profile[, -1]);
rownames (profile.data) <- profile[,1];
profile.data <- t(profile.data);
	
quantile.range <- quantile(profile.data, probs = seq(0, 1, 0.01));
if (quantile.range["90%"] > quantile.range["10%"])
{
	palette.breaks <- seq(quantile.range["10%"], quantile.range["90%"], 0.1);
}else
{
	palette.breaks <- seq(quantile.range["0%"], quantile.range["100%"], 0.1);	
}
	# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
color.palette  <- colorRampPalette(c("#99d594", "#ffffbf", "#fc8d59"))(length(palette.breaks) - 1);

	
pdf(file=outfile, w, h);
lmat = rbind(c(0,3),c(2,1),c(0,4));
lwid = c(0.5,4)
lhei = c(0.5,4,0.5)
heatmap.2(profile.data, main=RBP.name, cexRow=1, col=color.palette, dendrogram='none', Colv=NA, Rowv=NA, trace="none", density.info="none", margins=c(8,8), lmat=lmat, lwid=lwid, lhei=lhei);
dev.off();




