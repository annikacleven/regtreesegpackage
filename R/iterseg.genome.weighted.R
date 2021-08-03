#'Segment and Visualize a Whole Genome
#'
#' Creates a plot of the Whole Genome segmented, a dataframe including all the segments,
#' and a dataframe including the predictions from the regression tree using the optimal cp value for each chromosome,
#' iterative regression trees,and weighting to incentitze the regression to catch spikes.
#'
#' @param df A dataframe with columns of Start.Pos, log2r, and Chr columns. The Chr column should have format like "chr1", "chr21", "chrY".
#' @param png_filename A .png filename
#' @param upper.y.lim The upper limit for the y value of the plot
#' @param lower.y.lim The lower limit for the y value of the plot
#' @param cpvalue Specify a constant cp value for the regression tree to use instead of the optimal cp value
#'
#' @return A .png file with the whole genome plot and segmentation data,
#'  a list containing a dataframe with all of the segmentation data (segments),
#'  and a dataframe with the predictions from the regression tree (regtreepred)
#' @examples
#' example <- iterseg.genome.weighted(datafr127, png_filename = "datafr.png")
#' example$regtreepred
#' example$segments
#' @author Annika Cleven
#' @export

iterseg.genome.weighted <- function(df,png_filename, numiterations = 3, cpvalue = NA, upper.y.lim = 5, lower.y.lim = 5 ){

  upper.y.lim <- upper.y.lim
  lower.y.lim <- lower.y.lim
  png_filename <- png_filename

  ###Finding names

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

  cplist <- c()
  for(i in names2){
    subset <- df%>%
      dplyr::filter(Chr == i)

    if(is.na(cpvalue)){CP <- cpopt(subset)
    #print(CP)
    }

    else{CP <- cpvalue
    #print(CP)
    }

    #Optimal model with log2r and Start.Pos

    mad.diff = mad(diff(subset$log2r))

    if (mad.diff != 0){
      subset <- subset%>%
        dplyr::mutate(weight = ifelse(abs(log2r) <= mad.diff, 1, abs(log2r)/mad.diff))}

    else{subset <- subset%>%
      dplyr::mutate(weight = 1)}

    model1<- rpart::rpart(log2r~Start.Pos, subset, weights = subset$weight,
                          control=rpart.control(cp = CP))

    #Creating full_pred df
    if((subset$Chr == "chr1")){full_pred <- subset%>%
      modelr::add_predictions(model1)}

    else{pred_added <- subset%>%
      modelr::add_predictions(model1)
    full_pred <- rbind(full_pred, pred_added)}
  }

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
  cplist <- c()
  for (i in c(1:iterations)) {
    counterstr <- as.character(counter)
    errorcol <- paste("error",counterstr,sep="")

    df_split <- split(full_pred, full_pred$Chr)
    segments <- NULL

    for(chrid in names2){
      df_chr <- df_split[[chrid]]


      if(is.na(cpvalue)){CP <- cpopt(subset)
      #print(CP)
      }

      else{CP <- cpvalue
      #print(CP)
      }

      cplist <- c(cplist, CP)
      # fit the regression tree model

      mad.diff = mad(diff(df_chr[, ncol(df_chr)]))

      if(mad.diff != 0){

        df_chr <- df_chr%>%
          dplyr::mutate(weight = ifelse(df_chr[, ncol(df_chr)] <= mad.diff, 1, abs(df_chr[, ncol(df_chr)])/mad.diff))}

      else{df_chr%>%
          dplyr::mutate(weight = 1)}


      model<- rpart::rpart(df_chr[, ncol(df_chr)]~Start.Pos, df_chr, weights = df_chr$weight,
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

  lengthlist <- seq(2,length(full_pred$totalpred), by = 1)
  for(i in lengthlist){
    if(full_pred[i,ncol(full_pred)] != full_pred[i-1,ncol(full_pred)]){
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
  weightedavglog2ratio = start.df$totalpred
  location = .5 * ( StartIDs + EndIDs)
  width = EndIDs - StartIDs +1

  IDs <- data.frame(Chr, chrN, StartIDs, EndIDs, weightedavglog2ratio, location, width)


  listOfDataframe = list(
    "regtreepred" = full_pred,
    "segments" = IDs
  )

  numsegments <- paste("Number of Segments:",toString(nrow(IDs)))

  WholeGenome.Plot(png_filename, chr = full_pred$chrN, s = full_pred$totalpred, x = full_pred$log2r, segmentnumber = numsegments,
                   sample.name=deparse(substitute(df)), up.y=upper.y.lim, lo.y = lower.y.lim)

  return(listOfDataframe)


}
