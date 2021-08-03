#'Segment and Visualize a Chromosome
#'
#' Creates a plot of chromosome segmented, a dataframe including all the segments in the chromosome,
#' and a dataframe including the predictions from the regression tree using the optimal cp value for the chromosome.
#'
#' @param df A dataframe with columns of Start.Pos, log2r, and Chr columns. The Chr column should have format like "chr1", "chr21", "chrY".
#' @param chromid The chromosome id, for example "chr5"
#' @return A list containing a dataframe with all of the segmentation data (segments),
#' and a dataframe with the predictions from the regression tree (regtreepred),
#' and a plot of the segmentation and chromosome data
#' @examples
#' example <- seg.chr(datafr, "chr5")
#' example$regtreepred
#' example$segments
#' example$plot
#' @author Annika Cleven
#' @export

seg.chr <- function(df, chromid){
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
    if (is.na(cpvalue)){
    CP <- cpopt(subset)}
    else{CP <- cpvalue}

    model1<- rpart::rpart(log2r~Start.Pos, subset,
                          control=rpart.control(cp = CP))


    if((subset$Chr == "chr1")){full_pred <- subset%>%
      modelr::add_predictions(model1)}

    else{pred_added <- subset%>%
      modelr::add_predictions(model1)
    full_pred <- rbind(full_pred, pred_added)}
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
                         control=rpart.control(cp = CP2))

    if(!is.null(model$splits))
    {
      splitPoints = model$splits[, 4][order(model$splits[, 4])]
      startIDs = c(1, ceiling(splitPoints))
      startPos <- df_chr$Start.Pos[startIDs]
      endIDs = c(floor(splitPoints), nrow(df_tmp))
      endPos <-  df_chr$Start.Pos[endIDs]
      meanlog2ratio = sapply(1:length(startIDs), function(x)  mean(df_chr$log2r[startIDs[x]:endIDs[x]],na.rm = T))
      segments_chr <- data.frame(Chr = chrid, start = startPos, end = endPos, meanlog2ratio = meanlog2ratio, location  = 0.5*(startPos+endPos), widths = endPos - startPos +1)


      segments <- rbind(segments, segments_chr)
    }else{
      segments_chr <- data.frame( Chr = chrid,
                                  start = df_chr$Start.Pos[1], end = df_chr$Start.Pos[nrow(df_chr)],
                                  meanlog2ratio = mean(df_chr$log2r), location  = 0.5*(startPos+endPos), widths = endPos - startPos +1)
    }
  }
  segments <- segments %>%
    dplyr::mutate(chrom = gsub("chr", "", Chr),
                  chrN = ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrid)),
                  chrnumber = as.numeric(ifelse(Chr == "chrX", 23, ifelse(Chr == "chrY", 24, chrom))))

  segments <- segments%>%
    dplyr::select(Chr, chrnumber, start, end, meanlog2ratio, location, widths)

  segments$chrnumber <- as.numeric(as.character(segments$chrnumber))
  segments <- segments[order(segments$chrnumber),]

  chrsegments <- segments %>%
    dplyr::filter(Chr == chromid)

  chr_pred <- full_pred %>%
    dplyr::filter(Chr == chromid)

  plot <- chr_pred%>%
    ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = Start.Pos, y = log2r), color = "blue", size = 1)+
    ggplot2::geom_step(ggplot2::aes(Start.Pos, pred), color = "red", size = 1)+
    ggplot2::labs(title = paste(deparse(substitute(df)), ":", chromid), x = "Bin index", y = "Allele Imbalance",
                  subtitle = paste("Number of Segments:", nrow(chrsegments)))

  listOfDataframe = list(
    "regtreepred" = chr_pred,
    "segments" = chrsegments,
    "cpdf" = cpdf,
    "plot" = plot

  )

  return(listOfDataframe)
}
