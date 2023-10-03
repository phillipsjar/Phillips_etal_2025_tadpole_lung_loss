#function meant for cleaned up BT output of model with model string intact
RJ_model_testing = function(data, model = "dependent"){
  if(model == "dependent"){
    models <- unique(data$Model.string)
    model_matrix <- matrix(nrow = length(models), ncol = 3)
    colnames(model_matrix) <- c("model_string", "number_observed", "percent_of_posterior")
    rownames(model_matrix) <- 1:length(models)
    model_matrix[,1] <- models
    for(i in 1:length(models)){
      model_matrix[i,2] <- length(which(data$Model.string %in% model_matrix[i,1]))
    }
    runs <- dim(data)[1]
    
    model_matrix[,3] <- (as.numeric(model_matrix[,2])/runs)
    output <- matrix(NA,nrow = 2, ncol = 3)
    colnames(output) <- c("No regains", "No lentic regains", "No lentic losses")
    rownames(output) <- c("posterior odds", "BayesFactor")
    output[1,1] <- sum(as.numeric(model_matrix[grep("Z . . . . Z . .", model_matrix[,1]),3]))
    output[1,2] <- sum(as.numeric(model_matrix[grep(". . . . . Z . .", model_matrix[,1]),3]))
    output[1,3] <- sum(as.numeric(model_matrix[grep(". . . . . . . Z", model_matrix[,1]),3]))
    output[2,1] <- output[1,1]/(1-output[1,1])
    output[2,2] <- output[1,2]/(1-output[1,2])
    output[2,3] <- output[1,3]/(1-output[1,3])
    return(output)}
  if(model == "independent"){
    models <- unique(data$Model.string)
    model_matrix <- matrix(nrow = length(models), ncol = 3)
    colnames(model_matrix) <- c("model_string", "number_observed", "percent_of_posterior")
    rownames(model_matrix) <- 1:length(models)
    model_matrix[,1] <- models
    for(i in 1:length(models)){
      model_matrix[i,2] <- length(which(data$Model.string %in% model_matrix[i,1]))
    }
    runs <- dim(data)[1]
    model_matrix[,3] <- round((as.numeric(model_matrix[,2])/runs), 4)
    output <- matrix(NA,nrow = 2, ncol = 1)
    colnames(output) <- c("No regains")
    rownames(output) <- c("posterior odds", "BayesFactor")
    output[1,1] <- sum(as.numeric(model_matrix[grep(". . Z .", model_matrix[,1]),3]))
    output[2,1] <- output[1,1]/(1-output[1,1])
    
    return(output)}
  if(model == "six_state"){
    models <- unique(data$Model.string)
    model_matrix <- matrix(nrow = length(models), ncol = 3)
    colnames(model_matrix) <- c("model_string", "number_observed", "percent_of_posterior")
    rownames(model_matrix) <- 1:length(models)
    model_matrix[,1] <- models
    for(i in 1:length(models)){
      model_matrix[i,2] <- length(which(data$Model.string %in% model_matrix[i,1]))
    }
    runs <- dim(data)[1]
    
    model_matrix[,3] <- (as.numeric(model_matrix[,2])/runs)
    output <- matrix(NA,nrow = 2, ncol = 3)
    colnames(output) <- c("No regains", "No lentic regains", "No lentic losses")
    rownames(output) <- c("posterior odds", "BayesFactor")
    output[1,1] <- sum(as.numeric(model_matrix[grep(". Z . . Z . . . . . . Z . .", model_matrix[,1]),3]))
    output[1,2] <- sum(as.numeric(model_matrix[grep(". . . . . . . . . . . Z . .", model_matrix[,1]),3]))
    output[1,3] <- sum(as.numeric(model_matrix[grep(". . . . . . . . . . . . . Z", model_matrix[,1]),3]))
    output[2,1] <- output[1,1]/(1-output[1,1])
    output[2,2] <- output[1,2]/(1-output[1,2])
    output[2,3] <- output[1,3]/(1-output[1,3])
    return(output)}
}
