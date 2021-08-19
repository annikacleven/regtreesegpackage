#'Segment and Visualize a Single Chromosome
#'
#' Creates a plot of the chromosome of interest segmented, a dataframe including all the segments,
#' and a dataframe including the predictions from the regression tree using the optimal cp value for the chromosome,
#' iterative regression trees, and weighting to incentitze the regression to catch spikes.
#'
#' @param df A dataframe with columns of Start.Pos, log2r, and Chr columns.The Chr column should have format like "chr1", "chr21", "chrY".
#' @param chromid The chromosome of interest in the form "chr3" or "chr21"
#' @param cpvalue Specify a constant cp value for the regression tree to use instead of the optimal cp value
#' @param conserve  This will use the conservative cpopt method
#' @return A list containing a dataframe with all of the segmentation data (segments),
#' and a dataframe with the predictions from the regression tree (regtreepred), a dataframe with the cp values from the
#' first iteration (cpdf) a plot of ,the chromosome data and final predictions (chrplot),
#' a list 5 plots of each step in iteration (plots)
#' @examples
#' example <- iterseg.chr.weighted(s131L001, "chrX")
#' example$regtreepred
#' example$segments
#' example$chrplots
#' example$plots$iter1
#' @author Annika Cleven
#' @export


iterseg.chr.weighted <- function(df, chromid, cpvalue = NA, conserve = FALSE){
  ifelse(is.null(c(levels(df$Chr))), names <- c(unique(df$Chr)), names <- c(levels(df$Chr)))

  #Getting names
  names2 <- c()
  for(i in names){
    subset <- df%>%
      dplyr::filter(Chr == i)
    obs <- NROW(subset$Chr)
    if(obs>0){
      names2 <- c(names2, i)
    }
  }

  #Fitting the regression tree and getting the regression tree model
  cplist <- c()
  emptydf <- data.frame(matrix(ncol = ncol(df) + 1, nrow = 0))
  full_pred <- emptydf
  for(i in names2){
    subset <- df%>%
      dplyr::filter(Chr == i)


    if (!is.na(cpvalue)){
      CP <- cpvalue}
    else if (conserve == TRUE){CP <- cpopt(subset, conserve = TRUE)}
    else {CP <- cpopt(subset)}


    #Optimal model with log2r and Start.Pos

    mad.diff = mad(diff(subset$log2r))

    if (mad.diff != 0){
      subset <- subset%>%
        dplyr::mutate(weight = ifelse(abs(log2r) <= mad.diff, 1, abs(log2r)/mad.diff))}

    else{subset <- subset%>%
      dplyr::mutate(weight = 1)}


    model1<- rpart::rpart(log2r~Start.Pos, subset, weights = subset$weight,
                          control=rpart.control(cp = CP))


    pred_added <- subset%>%
      modelr::add_predictions(model1)

    full_pred <- rbind(full_pred, pred_added)
    cplist <- c(cplist, CP)

  }



  #Calculating the Error between the model and the log2r
  #Creating chrN so that chromosomes can be sorted numerically
  full_pred <- full_pred %>%
    dplyr::mutate(chrom = gsub("chr", "", Chr),
                  chrN = as.numeric(ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrom))),
                  error1 = log2r - pred)

  #One iteration was done by code above now there are 2 left
  iterations = 2
  full_pred <- as.data.frame(full_pred)

  ####Other iterations
  counter = 1
  pred = 2

  for (i in c(1:iterations)) {
    counterstr <- as.character(counter)
    errorcol <- paste("error",counterstr,sep="")

    df_split <- split(full_pred, full_pred$Chr)
    dflength<- length(df_split[[1]])
    segments <- NULL


    emptydf2 <-  data.frame(matrix(ncol = dflength, nrow = 0))
    full_pred <- emptydf2
    for(chrid in names2){
      df_chr <- df_split[[chrid]]


      if (!is.na(cpvalue)){
        CP <- cpvalue}
      else if (conserve == TRUE){CP <- cpopt(subset, conserve = TRUE)}
      else {
        CP <- cpopt(df_chr)}


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



      df_chr[[predcol]] <- predict(model, df_chr)
      pred_added <- df_chr
      full_pred <- rbind(full_pred, pred_added)


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

  full_pred <- full_pred %>%
    filter(Chr == chromid)

  #Getting segments

  segmentdf <- full_pred[1,]

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

  segments <- data.frame(Chr, chrN, StartIDs, EndIDs, log2ratio, location, width)


  ##VISUALIZATIONS


  df = full_pred

  a <- df %>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = log2r), size = .5, color = "blue")+
    ggplot2::geom_step(ggplot2::aes(Start.Pos, pred), color = "red", size = 1)+
    ggplot2::labs(title = "Predictions After One Iteration", x = "Bin Index", y = "CNV Profile")


  b <- df%>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = error1), size = .5, color = "blue")+
    ggplot2::geom_step(ggplot2::aes(Start.Pos, pred2), color = "red", size = 1)+
    ggplot2::labs(title = "Residual Error Predictions", x = "Bin Index", y = "CNV Profile")


  c <- df%>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = log2r), size = .5, color = "blue")+
    ggplot2::geom_step(ggplot2::aes(Start.Pos, pred + pred2), color = "red", size = 1)+
    ggplot2::labs(title = "Predictions After Two Iterations", x = "Bin Index", y = "CNV Profile")


  d <- df%>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = error2), size = .5, color = "blue")+
    ggplot2::geom_step(ggplot2::aes(Start.Pos,pred3), color = "red", size = 1)+
    ggplot2::labs(title = "Residual Error Predictions", x= "Bin Index", y = "CNV Profile")


  e <- df%>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = log2r), size = .5, color = "blue")+
    ggplot2::geom_step(ggplot2::aes(Start.Pos, pred + pred2 + pred3), color = "red", size = 1)+
    ggplot2::labs(title = paste("Chromid:", chromid, ", Number of segments:", nrow(segments)),
                  subtitle = "Predictions After Three Iterations", x = "Bin Index", y = "CNV Profile")

  plotlist = list("iter1" = a, "error1" = b, "iter2" = c, "error2" = d, "iter3"= e )

  cpdf <- data.frame(names2, cplist)

  listOfDataframe = list(
    "regtreepred" = full_pred,
    "segments" = segments,
    "cpdf" = cpdf,
    "chrplot" = e,
    "plots" = plotlist
  )

  return(listOfDataframe)

}
