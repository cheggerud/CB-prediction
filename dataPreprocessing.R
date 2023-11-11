# this should be dataForPrediction.csv.csv
data <- read.csv(file = "data.csv", header = TRUE, sep = ",")
dataBeforeFactorizing = data
# This is only for data preprocessing
if (F) {
  data[,1:3] = NULL
  names(data)[c(1:2,5)] = c("lake", "Month","Target")
  data = data[,c(1:4,6:20,5)]
  data = data[,c(1,3,2,4:20)]
  data$Target = data$Target - 1
  ######################### Save the discretized data
  write.csv(data, file = "data.csv", row.names = FALSE)   
}