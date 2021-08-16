#'Segment and Visualize a Whole Genome
#'
#' Creates a plot of the Whole Genome segmented, a dataframe including all the segments,
#' and a dataframe including the predictions from the regression tree using the optimal cp value for each chromosome
#'
#' @param df A dataframe with columns of Start.Pos, log2r, and Chr columns.  The Chr column should have format like "chr1", "chr21", "chrY".
#' @param png_filename A .png filename
#' @param upper.y.lim The upper limit for the y value of the plot
#' @param lower.y.lim The lower limit for the y value of the plot
#' @param cpvalue  Specified constant cp value for the regression tree to use instead of the optimal cp value
#' @return A .png file with the whole genome plot and segmentation,
#' a list containing a dataframe with all of the segmentation data(segments),
#' and a dataframe with the predictions from the regression tree(regtreepred), and a dataframe with the cp used
#' for each chromosome (cpdf)
#' @examples
#' example <- seg.genome(datafr127, png_filename = "datafr.png",upper.y.lim = 5, lower.y.lim = -5)
#' example$regtreepred
#' example$segments
#' @author Annika Cleven
#' @export


seg.genome <- function(df, png_filename, upper.y.lim = 5, lower.y.lim = -5, cpvalue = NA, conserve = FALSE){
  upper.y.lim <- upper.y.lim
  lower.y.lim <- lower.y.lim
  ifelse(is.null(c(levels(df$Chr))), names <- c(unique(df$Chr)), names <- c(levels(df$Chr)))

  names2 <- c()
  for(i in names){
    subset <- df%>%
      dplyr::filter(Chr == i)
    obs <- NROW(subset$Chr)
    if(obs>0){
      names2 <- c(names2, i)
    }
  }
  emptydf <- data.frame(matrix(ncol = ncol(df) + 1, nrow = 0))
  full_pred <- emptydf
  cplist <- c()
  for(i in names2){
    subset <- df%>%
      dplyr::filter(Chr == i)
    #if (cpvalue == NA){
    if (!is.na(cpvalue)){
    CP <- cpvalue}
    else if (conserve == TRUE){CP <- cpopt(subset, conserve = TRUE)}
    else {CP <- cpopt(subset)}



    model1<- rpart::rpart(log2r~Start.Pos, subset,
                   control=rpart.control(cp = CP))

    firstrow <- subset[1,]
    #print(firstrow)
    # if((firstrow$Chr == "chr1")){full_pred <- subset%>%
    #   modelr::add_predictions(model1)}

    # else{
      pred_added <- subset%>%
      modelr::add_predictions(model1)
    full_pred <- rbind(full_pred, pred_added)
    # }
    cplist <- c(cplist, CP)

  }

  full_pred <- full_pred %>%
    dplyr::mutate(chrom = gsub("chr", "", Chr),
           chrN = as.numeric(ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrom))))

  full_pred$chrN <- as.numeric(as.character(full_pred$chrN))
  full_pred <- full_pred[order(full_pred$chrN),]


  cpdf <- data.frame(names2, cplist)

  counter <- 1
  df_split <- split(df, df$Chr)
  segments <- NULL
  for(chrid in names2)
  {
    #print(paste('segment chr', chrid))
    df_chr <- df_split[[chrid]]
    df_tmp <- data.frame(Start.Pos = 1:NROW(df_chr), log2r = df_chr$log2r)
    CP2 <- cpdf[counter, 2]
    #print(CP2)
    counter <- counter +1
    # fit the regression tree model
    model<- rpart::rpart(log2r~Start.Pos, df_tmp,
                  control=rpart.control(cp = CP2, maxdepth = 10))

    if(!is.null(model$splits))
    {
      splitPoints = model$splits[, 4][order(model$splits[, 4])]
      startIDs = c(1, ceiling(splitPoints))
      startPos <- df_chr$Start.Pos[startIDs]
      endIDs = c(floor(splitPoints), nrow(df_tmp))
      endPos <-  df_chr$Start.Pos[endIDs]
      meanlog2ratio = sapply(1:length(startIDs), function(x)  mean(df_chr$log2r[startIDs[x]:endIDs[x]],na.rm = T))
      segments_chr <- data.frame(chr = chrid, start = startPos, end = endPos, meanlog2ratio = meanlog2ratio, location  = 0.5*(startPos+endPos), widths = endPos - startPos +1)


      segments <- rbind(segments, segments_chr)
    }else{
      segments_chr <- data.frame( chr = chrid,
                                  start = df_chr$Start.Pos[1], end = df_chr$Start.Pos[nrow(df_chr)],
                                  meanlog2ratio = mean(df_chr$log2r), location  = 0.5*(startPos+endPos), widths = endPos - startPos +1)
    }
  }
  segments <- segments %>%
    dplyr::mutate(chrom = gsub("chr", "", chr),
           chrN = ifelse(chr == "chrX", 23, ifelse(chr == "chrY", 24, chrid)),
           chrnumber = as.numeric(ifelse(chr == "chrX", 23, ifelse(chr == "chrY", 24, chrom))))

  segments <- segments%>%
    dplyr::select(chr, chrnumber, start, end, meanlog2ratio, location, widths)

  segments$chrnumber <- as.numeric(as.character(segments$chrnumber))
  segments <- segments[order(segments$chrnumber),]


  numsegments <- paste("Number of Segments:",toString(nrow(segments)))
  numsegments2 <- numsegments
  png_filename <- png_filename


  WholeGenome.Plot(png_filename, chr = full_pred$chrN ,x = full_pred$log2r, s = full_pred$pred, segmentnumber = numsegments, up.y=upper.y.lim, lo.y = lower.y.lim, is.chr.labeled=FALSE, is.txt.output=TRUE, output.digit=3, sample.name=deparse(substitute(df)))

  listOfDataframe = list(
    "regtreepred" = full_pred,
    "segments" = segments,
    "cpdf" = cpdf
  )

  return(listOfDataframe)

}


