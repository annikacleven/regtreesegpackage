#'Segment and Visualize a Whole Genome
#'
#'Creates a plot of the Whole Genome segmented, a dataframe including all the segments,
#'and a dataframe including the predictions from the regression tree using the optimal cp value for each chromosome
#'and iterative regression trees.
#'
#' @param df A dataframe with columns of Start.Pos, log2r, and Chr columns.  The Chr column should have format like "chr1", "chr21", "chrY".
#' @param png_filename A .png filename
#' @param upper.y.lim The upper limit for the y value of the plot
#' @param lower.y.lim The lower limit for the y value of the plot
#' @param conserve This will use the conservative cpopt method
#' @param cpvalue Specified constant cp value for the regression tree to use instead of the optimal cp value
#'
#' @return A .png file with the whole genome plot and segmentation data,
#'  a list containing a dataframe with all of the segmentation data (segments),
#'  a dataframe with the predictions from the regression tree (regtreepreds),
#'  and a dataframe with the cp values used in the first iteration (cpdf)
#' @examples
#' example <- iterseg.genome(s131L001, "example.png")
#' example$regtreepred
#' example$segments
#' @author Annika Cleven
#' @export



iterseg.genome <- function(df,png_filename, cpvalue = NA, conserve = FALSE, upper.y.lim = 5, lower.y.lim = 5){
  upper.y.lim <- upper.y.lim
  lower.y.lim <- lower.y.lim
  png_filename <- png_filename

  ###Finding Chromosome Names
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

  #Fitting Regression tree for each chromosome
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

    #Optimal model with log2r and Start.Pos

    model1<- rpart::rpart(log2r~Start.Pos, subset,
                          control=rpart.control(cp = CP))

    #Creating full_pred df
    firstrow <- subset[1,]
    # if((firstrow$Chr == "chr1")){full_pred <- subset%>%
    #   modelr::add_predictions(model1)}

    #else{
      pred_added <- subset%>%
      modelr::add_predictions(model1)
    full_pred <- rbind(full_pred, pred_added)
    #}
    cplist <- c(cplist, CP)
  }
  #cplist <- c(cplist, CP)

  #Calculating the Error between the model and the log2r
  #Creating chrN so that chromosomes can be sorted numerically
  full_pred <- full_pred %>%
    dplyr::mutate(chrom = gsub("chr", "", Chr),
                  chrN = as.numeric(ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrom))),
                  error1 = log2r - pred)

  #One iteration was done by code above now there are numiterations -1 left
  iterations = 2

  full_pred <- as.data.frame(full_pred)

  ####Other iterations
  counter = 1
  pred = 2

  for (i in c(1:iterations)) {
    counterstr <- as.character(counter)
    errorcol <- paste("error",counterstr,sep="")

    df_split <- split(full_pred, full_pred$Chr)
    segments <- NULL
    #Finding the chromosome names left
    for(chrid in names2){
      df_chr <- df_split[[chrid]]


      if (!is.na(cpvalue)){
        CP <- cpvalue}
      else if (conserve == TRUE){CP <- cpopt(subset, conserve = TRUE)}
      else {CP <- cpopt(subset)}


      # fit the regression tree model

      model<- rpart::rpart(df_chr[, ncol(df_chr)]~Start.Pos, df_chr,
                           control=rpart.control(cp = CP))

      predstr <- toString(pred)
      predcol <- paste("pred",predstr, sep = "")

      if((chrid == "chr1")){
        df_chr[[predcol]] <- predict(model, df_chr)
        full_pred <- df_chr

      }

      else {
        df_chr[[predcol]] <- predict(model, df_chr)
        pred_added <- df_chr
        full_pred <- rbind(full_pred, pred_added)
      }
    }

    counter = counter + 1
    counterstr <- as.character(counter)
    errorcol <- paste("error",counterstr,sep="")

    full_pred[[errorcol]] <- ((full_pred[, ncol(full_pred)-1] - full_pred[, ncol(full_pred)]))

    pred = pred + 1

  }

  full_pred$chrN <- as.numeric(as.character(full_pred$chrN))
  full_pred <- full_pred[order(full_pred$chrN),]


  full_pred <- full_pred %>%
    dplyr::mutate(totalpred = pred + pred2 + pred3)

  segmentdf <- full_pred[1,]

  #Finding each segment

  lengthlist <- seq(2,length(full_pred$totalpred), by = 1)
  for(i in lengthlist){
    if(full_pred[i,"totalpred"] != full_pred[i-1,"totalpred"]){
      segmentdf <- rbind(segmentdf, full_pred[i-1,])
      segmentdf <- rbind(segmentdf, full_pred[i,])}

  }
  segmentdf <- rbind(segmentdf, full_pred[nrow(full_pred),])

  segmentdf <- segmentdf %>%
    dplyr::select(Chr, Start.Pos, chrN, totalpred)%>%
    dplyr::mutate(row_num <- seq.int(nrow(segmentdf)))


  end.df <- segmentdf %>% dplyr::filter(row_number() %% 2 == 0)
  start.df <- segmentdf %>% dplyr::filter(row_number() %% 2 == 1)

  StartIDs = start.df$Start.Pos
  EndIDs = end.df$Start.Pos
  Chr = start.df$Chr
  chrN = start.df$chrN
  log2ratio = start.df$totalpred
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
  WholeGenome.Plot(png_filename, chr = full_pred$chrN, s = full_pred$totalpred, x = full_pred$log2r, segmentnumber = numsegments, sample.name=deparse(substitute(df)))

  return(listOfDataframe)

}
