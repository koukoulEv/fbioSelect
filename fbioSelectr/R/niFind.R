#' Compute regression weights for repeated measures data
#'
#' Takes the data and subject ids and creates the independence weights.
#' @param data a dataframe including the response vector.
#' @param idLabel the name of the variable in the dataframe containing subjects ids.
#' @return Independence structure repeated measures weights.
#' @export
niFind = function(data, idLabel){

  freq = data.frame(table(data[, idLabel]))
  names(freq) = c(idLabel, "ni")
  newData = merge(data.frame(id = data[,idLabel]), freq, by = c(idLabel))
  return(newData$ni)
}
