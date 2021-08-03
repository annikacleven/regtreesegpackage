#' Plotting Whole Genome
#'
#' This plots whole genome data with options to smooth, color, and other customizations
#'
#' @param png_filename A .png filename that the plot will be exported to
#' @param chr The chromosome the datapoint is a part of
#' @param x The data point value that is being plotted
#' @param s The smoothing value
#' @param segmentnumber The number of segments in the tree, the value here will be printed under the title of the plot.
#' @param sample.name The name of the sample
#' @param ylab The y axis title
#' @param up.y The y axis upper limit
#' @param lo.y The y axis lower limit
#' @param odd.col The color of the odd numbered chromosomes
#' @param even.col The color of the even numbered chromosomes
#' @param x.col The color of the x chromsome
#' @param y.col The color of the y chromosome
#' @param plot.bg.color The color of the plot background
#' @param med.line.col The color of the median line
#' @return A .png file with the whole genome data
#' @examples
#'   WholeGenome.Plot("test.png", chr = data$Chr ,x = data$log2r, up.y=5, lo.y = -5)
#' @author Annika Cleven
#' @export


WholeGenome.Plot <- function(png_filename, chr, x, s=NA, segmentnumber = NA,
                             sample.name="sample", ylab = "Allele Imbalance", BAF = F,
                             up.y=3, lo.y = -3, plot.log = "",
                             odd.col='magenta', even.col='dodgerblue', x.col = 'blueviolet', y.col = "darkmagenta", plot.cex=3, plot.pch = ".", plot.bg.color = "white",
                             is.chr.labeled=TRUE, is.txt.output=TRUE, output.digit=3,
                             odd.y.pos=-10,even.y.pos=-38,label.cex=1, med.line.col = "red")
{
  if(is.character(chr[1])){
    chr <- gsub("chr", "", chr)
    chr <- gsub("Chr", "", chr)
  }
  chrN <- chr;  chrN[chr == "chrX"] <- "23"; chrN[chr == "chrY"] <- "24"
  chrN <- as.numeric(chrN)

  ## order data by chr numbers
  chr <- chr[order(chrN)]
  x <- x[order(chrN)]
  s <- s[order(chrN)]
  chrN <- chrN[order(chrN)]

  chrN <- chrN

  color_vec <- rep("black",length(chr))
  color_vec[which(chrN %% 2 != 0)] <- odd.col
  color_vec[which(chrN %% 2 == 0)] <- even.col
  color_vec[which(chrN == 23)] <- x.col
  color_vec[which(chr == 24)] <- y.col


  if(BAF){
    nx <- -1*x+1
    ns <- -1*s+1
    lo.y <- 0
    up.y <- 1
  }

  png(png_filename,width=3.6e3,height=1.1e3,res=210)
  if(is.null(up.y)){up.y <- quantile(x,0.995,na.rm=TRUE)*2.5}

  plot(x,cex=plot.cex, pch = plot.pch, bg=plot.bg.color,
       col = color_vec, ylim = c(lo.y, up.y),
       xlab = 'Bin index',ylab = ylab,
       main=paste('profile of',sample.name), log = plot.log)

  mtext(side = 3, paste(segmentnumber))

  if(BAF){
    points(nx,cex=plot.cex, pch = plot.pch, col = color_vec)
  }

  ## add overflow points back
  ##
  ox <- rep(NA, length(x))
  #ox[x <= lo.y] <- -3; ox[x >= up.y] <- 3
  opch <- ifelse(ox == -3, 6, ifelse(ox == 3, 2, NA))

  if(!all(is.na(ox))){
    points(ox, cex=1, pch = opch, col = color_vec)
  }

  if(!all(is.na(s))){
    lines(s)
    if(BAF){
      lines(ns)
    }
  }

  if(BAF){
    abline(h=0.5, lwd=3, col=med.line.col, lty=3)
  } else {
    abline(h=median(x, na.rm = TRUE),lwd=3,col=med.line.col,lty=3)
  }

  chr.Str <- c( seq(1,22), "X", "Y")
  if(is.chr.labeled)
  {
    ## add chr labels and separation lines
    temp <- data.frame(idx = 1:length(x), chr = chrN)
    median_pos <- aggregate(idx ~ chrN, data = temp, median)
    max_pos <- aggregate(idx ~ chrN, data = temp, max)
    tmp.y.pos <- rep(c(odd.y.pos, even.y.pos), length(unique(chr)))[1:length(unique(chr))]
    text(x=median_pos$idx,y=-tmp.y.pos,labels= unique(chr), cex=label.cex)
    text(x=median_pos$idx,y=tmp.y.pos,labels= unique(chr), cex=label.cex)
    abline( v =  max_pos$idx[-nrow(max_pos)], lwd = 3, lty = 2, col = "grey")
  }

  dev.off()
}

