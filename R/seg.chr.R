#'Segment and Visualize a Chromosome
#'
#' Creates a plot of chromosome segmented, a dataframe including all the segments in the chromosome,
#' and a dataframe including the predictions from the regression tree using the optimal cp value for the chromosome.
#'
#' @param df A dataframe with columns of Start.Pos, log2r, and Chr columns. The Chr column should have format like "chr1", "chr21", "chrY".
#' @param chromid The chromosome id, for example "chr5"
#' @param cpvalue Specified constant cp value for the regression tree to use instead of the optimal cp value
#' @param conserve This will use the conservative cpopt method
#' @return A list containing a dataframe with all of the segmentation data (segments),
#' and a dataframe with the predictions from the regression tree (regtreepred),
#' a plot of the segmentation and chromosome data, and a dataframe with the optimal cp used
#' for each chromosome on the first round of iteration(cpdf)
#' @examples
#' example <- seg.chr(datafr, "chr5")
#' example$regtreepred
#' example$segments
#' example$plot
#' @author Annika Cleven
#' @export
#'
#'


seg.chr <- function(df, chromid, cpvalue = NA, conserve = FALSE){
  ifelse(is.null(c(levels(df$Chr))), names <- c(unique(df$Chr)), names <- c(levels(df$Chr)))


#Finding Names
  names2 <- c()
  for(i in names){
    subset <- df%>%
      dplyr::filter(Chr == i)
    obs <- NROW(subset$Chr)
    if(obs>0){
      names2 <- c(names2, i)
    }
  }

  #Fitting the regression tree and getting predictions
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

    model1<- rpart::rpart(log2r~Start.Pos, subset,
                          control=rpart.control(cp = CP))



    pred_added <- subset%>%
      modelr::add_predictions(model1)

    full_pred <- rbind(full_pred, pred_added)
    cplist <- c(cplist, CP)}


  full_pred <- full_pred %>%
    dplyr::mutate(chrom = gsub("chr", "", Chr),
                  chrN = as.numeric(ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrom))))

  full_pred$chrN <- as.numeric(as.character(full_pred$chrN))
  full_pred <- full_pred[order(full_pred$chrN),]
  full_pred <- full_pred %>%
    filter(Chr == chromid)


  #######
  #Getting segments
  segmentdf <- full_pred[1,]


  for(i in (2:length(full_pred$pred))){
    if(full_pred[i,"pred"] != full_pred[i-1,"pred"]){
      segmentdf <- rbind(segmentdf, full_pred[i-1,])
      segmentdf <- rbind(segmentdf, full_pred[i,])}
  }


  segmentdf <- rbind(segmentdf, full_pred[nrow(full_pred),])


  segmentdf <- segmentdf %>%
    dplyr::select(Chr, Start.Pos, chrN, pred)

  end.df <- segmentdf %>% dplyr::filter(row_number() %% 2 == 0)

  start.df <- segmentdf %>% dplyr::filter(row_number() %% 2 == 1)


  StartIDs = start.df$Start.Pos
  EndIDs = end.df$Start.Pos
  Chr = start.df$Chr
  chrN = start.df$chrN
  log2ratio = start.df$pred
  location = .5 * ( StartIDs + EndIDs)
  width = EndIDs - StartIDs +1

  segments <- data.frame(Chr, chrN, StartIDs, EndIDs, log2ratio, location, width)

  #Tracking the cp values
  cpdf <- data.frame(names2, cplist)

  #Plotting
  chrplot <- full_pred%>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = log2r), color = "blue", size = 1)+
    ggplot2::geom_step(ggplot2::aes(Start.Pos, pred), color = "red", size = 1)+
    ggplot2::labs(title = paste(deparse(substitute(df)), ":", chromid), x = "Bin index", y = "Allele Imbalance",
                  subtitle = paste("Number of Segments:", nrow(segments)))

  listOfDataframe = list(
    "regtreepred" = full_pred,
    "segments" = segments,
    "cpdf" = cpdf,
    "chrplot" = chrplot

  )

  return(listOfDataframe)
}
