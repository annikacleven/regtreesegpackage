#' Find the optimal complexity parameter value
#'
#' The function finds the complexity parameter that is most optimal for reducing cross validation error.
#'
#' @param df A dataframe with columns of Start.Pos and log2r
#' @param conserve The default is FALSE.  When TRUE the optimal complexity parameter has a minimum around .01 to reduce the chance of overfitting.
#' When FALSE the optimal complexity has a minimum around .001.
#' @return The optimal complexity parameter in numeric format
#' @examples
#'  my <- cpopt(df)
#'  print(my)
#' @author Annika Cleven
#' @export

cpopt <- function(df, conserve = FALSE){

  if(conserve == FALSE){
    model<-rpart::rpart(log2r~Start.Pos, data = df, control = rpart.control(cp = .01))
    df2 = as.data.frame(model[["frame"]][["complexity"]])

    nrow <- nrow(df2)
    smallcp <-df2[nrow, ]

    model1b <- rpart::rpart(log2r~Start.Pos, data = df, control = rpart.control(cp = smallcp))
    cp <- as.data.frame(model1b$cptable[,1])
    cp<- cp%>%
      dplyr::mutate(id = as.numeric(rownames(cp)),
             cp = model1b$cptable[,1])

    xerror <- as.data.frame(model1b$cptable[,4])
    xerror <- xerror %>%
      dplyr::mutate(id = as.numeric(rownames(cp)),
             xerror = model1b$cptable[,4])

    relerror <- as.data.frame(model1b$cptable[,3])
    relerror <- relerror %>%
      dplyr::mutate(id = as.numeric(rownames(cp)),
             relerror = model1b$cptable[,3])

    datafr <- full_join(cp, xerror)
    full_df <- full_join(datafr, relerror)

    id <- which.min(full_df$xerror)

    cpopt <- full_df$cp[id]

  }

  if(conserve == TRUE){
    model<-rpart::rpart(log2r~Start.Pos, data = df, control = rpart.control(cp = .01))
    cpVals <- model$cptable
    id <- which.min(cpVals[,4])
    cpopt2 <- cpVals[id,1]
    return(as.numeric(cpopt2))
  }


  return(as.numeric(cpopt))
}
