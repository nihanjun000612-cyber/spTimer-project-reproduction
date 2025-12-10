#setwd('D:/data')
carData=read.table("car.data", 
                      header = FALSE,
                      sep=',',
                      col.names = c("1", "2", "3", "4", 
                                    "5", "6", "7"),
                      stringsAsFactors = FALSE)
carData3=factor(carData[, 3], 
                       levels = c("2", "3", "4", "5more"),
                       ordered = TRUE)
carData4=factor(carData[, 4],
                       levels = c("2", "4", "more"),
                       ordered = TRUE)
ordCarData1=carData[order(carData[, 3], carData[, 4]), ]
ordCarData2=carData[order(carData[, 4], carData[, 3]), ]
day1Data <- data.frame(idNum = 1:10, 
                       measure = rnorm(10))
day2Data <- data.frame(idNum = 11:20, 
                       measure = rnorm(10))
stackedData=rbind(day1Data, day2Data)
sideBySide=cbind(day1Data, day2Data)