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

  #Finding chromosome names
  names2 <- c()
  for(i in names){
    subset <- df%>%
      dplyr::filter(Chr == i)
    obs <- NROW(subset$Chr)
    if(obs>0){
      names2 <- c(names2, i)
    }
  }

  #Fitting each a regression tree for each chromosome
  emptydf <- data.frame(matrix(ncol = ncol(df) + 1, nrow = 0))
  full_pred <- emptydf
  cplist <- c()
  for(i in names2){
    subset <- df%>%
      dplyr::filter(Chr == i)
    if (!is.na(cpvalue)){
      CP <- cpvalue}
    else if (conserve == TRUE){CP <- cpopt(subset, conserve = TRUE)}
    else {CP <- cpopt(subset)}

    model1<- rpart::rpart(log2r~Start.Pos, subset,
                          control=rpart.control(cp = CP))

    firstrow <- subset[1,]

    pred_added <- subset%>%
      modelr::add_predictions(model1)
    full_pred <- rbind(full_pred, pred_added)
    cplist <- c(cplist, CP)

  }

  full_pred <- full_pred %>%
    dplyr::mutate(chrom = gsub("chr", "", Chr),
                  chrN = as.numeric(ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrom))))

  full_pred$chrN <- as.numeric(as.character(full_pred$chrN))
  full_pred <- full_pred[order(full_pred$chrN),]


  cpdf <- data.frame(names2, cplist)

  segmentdf <- full_pred[1,]

  #Finding each segment

  lengthlist <- seq(2,length(full_pred$pred), by = 1)
  for(i in lengthlist){
    if(full_pred[i,"pred"] != full_pred[i-1,"pred"]){
      segmentdf <- rbind(segmentdf, full_pred[i-1,])
      segmentdf <- rbind(segmentdf, full_pred[i,])}

  }
  segmentdf <- rbind(segmentdf, full_pred[nrow(full_pred),])

  segmentdf <- segmentdf %>%
    dplyr::select(Chr, Start.Pos, chrN, pred)%>%
    dplyr::mutate(row_num <- seq.int(nrow(segmentdf)))


  end.df <- segmentdf %>% dplyr::filter(row_number() %% 2 == 0)
  start.df <- segmentdf %>% dplyr::filter(row_number() %% 2 == 1)

  StartIDs = start.df$Start.Pos
  EndIDs = end.df$Start.Pos
  Chr = start.df$Chr
  chrN = start.df$chrN
  log2ratio = start.df$pred
  location = .5 * ( StartIDs + EndIDs)
  width = EndIDs - StartIDs +1

  IDs <- data.frame(Chr, chrN, StartIDs, EndIDs, log2ratio, location, width)
  cpdf <- data.frame(names2, cplist)


  listOfDataframe = list(
    "regtreepred" = full_pred,
    "segments" = IDs,
    "cpdf" = cpdf
  )
  #Calculating number of segments
  numsegments <- paste("Number of Segments:",toString(nrow(IDs)))
  #Plotting the segmentation
  WholeGenome.Plot(png_filename, chr = full_pred$chrN, s = full_pred$pred, x = full_pred$log2r, segmentnumber = numsegments, sample.name=deparse(substitute(df)))

  return(listOfDataframe)


}


