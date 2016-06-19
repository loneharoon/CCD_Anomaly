executeARIMA_withoutTemperature <- function(train_data,test_data){
  #this function uses univarite energy data only to forecast energy consumption
  library(forecast)
  numcores <- detectCores() - 1 #no. of cores used for computation
  lastday_traindata <- as.Date(last(index(train_data)),tz="Asia/Kolkata") # last day of training data
  noOfDays <-length(unique(as.Date(index(test_data),tz="Asia/Kolkata"))) # no. of predicting days
  forecastday <- lastday_traindata+1 # day for which forecast is needed
  for(i in 1:noOfDays){
    # create time series object, 48 is seasonality i.e, 48 observations contstitue one season
    tsobject <- ts(coredata(train_data$power), frequency = 144)
    arimaObject <- Arima(tsobject,order=c(4,1,2),seasonal=c(2,1,0),method="CSS")
    # extract actual power consumption and temperature on forecast day
    testday_data <- test_data[as.Date(index(test_data),tz="Asia/Kolkata") %in% forecastday]
    p2 <- forecast(arimaObject,h=144,level = 95)

    if (!exists('res_obj')){
      res_obj <- xts(data.frame(actual=testday_data$power,forecast=p2$mean,lower=p2$lower,upper=p2$upper),index(testday_data))
    }else{
      temp <- xts(data.frame(actual=testday_data$power,forecast=p2$mean,lower=p2$lower,upper=p2$upper),index(testday_data))
      res_obj <- rbind(res_obj,temp)
    }
    print(as.Date(index(head(testday_data,1))))
    train_data <- rbind(train_data,testday_data) # append actual data of forecast day with previous training data
    train_data <- train_data[-c(1:144),] # remove first day data from training data as extra day data is appended
    forecastday <- forecastday + 1 # index of next forecast day
  }
  return(res_obj)
}

find_arimaresult_status <- function(arima_result){
  # this function is used to check whether for n instances within an interval anomaly status was high
suppressWarnings(rm(park,info_dframe))
  anomalythreshold_len <- 4
pastanomaly_vec <- as.numeric()
for(i in 1:NROW(arima_result)){
  if(.indexmin(arima_result[i])%in%c(0,10,20,30,40))
  {
    if(exists("park"))
      park<-rbind(park,arima_result[i])
    else
      park <- arima_result[i]
  }else{
    park<-rbind(park,arima_result[i])
    decision_vector <- c(ifelse((coredata(park$actual) > coredata(park$upper)),1,0 ))
    decision_vector <- c(pastanomaly_vec,decision_vector) # May be in the past interval there were some anomaly continuing
    rm(park)
    length_list <- rle(decision_vector) # run length encoding
    is_anom <- any(length_list$lengths >= anomalythreshold_len & length_list$values == 1) # http://stackoverflow.com/a/37571854/3317829

    ## this part of code saves the status of current interval instances with anomaly
    cnt <- 0; lcs <-0 #last common sequence with status as 1
    last <- decision_vector[length(decision_vector)-cnt]
    print(index(arima_result[i]))
    while(!is.null(last) & last==1){
      lcs <- lcs + 1
      cnt <- cnt +1
      if(lcs >= anomalythreshold_len)
        break
      last <- decision_vector[length(decision_vector)-cnt]
      #browser()
    }
    pastanomaly_vec <- rep(1,lcs)
    if(length(pastanomaly_vec) >= anomalythreshold_len)
      pastanomaly_vec <- as.numeric() # this sequence is already detected in current interval no need to carry it forward
    ## this part ends here
    if(!exists("info_dframe")){
      info_dframe <- xts(is_anom,index(arima_result[i]))
    }else{
      temp <- xts(is_anom,index(arima_result[i]))
      info_dframe <- rbind(info_dframe,temp)
    }
  }

}
  return(info_dframe)
}

plot_arimaresults <- function(arima_result,fname){
  #this function prints predicted,actual and 2standard deviation values
  fname = paste0(strsplit(filename,"[.]")[[1]][1],"_WithARIMA.pdf")
  pdf(file = paste0("Dropbox/",fname),width = 14,height = 10)
  par(mfrow=c(4,3))
  len<- NROW(arima_result)
  par(mfrow=c(4,3))
  for(i in 1:len){
    if(.indexmin(arima_result[i])%in% c(0,10,20,30,40)){
      if(!exists("plotframe")){
        plotframe <- arima_result[i]
      }else{
        temp <- arima_result[i]
        plotframe <- rbind(plotframe,temp)
      }
    }else{
      plotframe <- rbind(plotframe,arima_result[i])
      plot(plotframe$actual/1000,ylim=c(0,4),main="")
      lines(plotframe$fit/1000,col="blue")
      # # lines(temp2/1000,col="blue")
      lines(plotframe$upper/1000,col="red",lt=3)
      rm(plotframe)
    }
  }
  dev.off()
}

metric_SMAPE <- function(object){
  numer <- sum(abs(object$fit-object$actual))
  denom <- sum(abs(object$fit)+abs(object$actual))
  smape <- numer/denom
  return (smape)
}
metric_MASE <- function(forecast_ob){
  # https://en.wikipedia.org/wiki/Mean_absolute_scaled_error 
  numerator <- abs(forecast_ob$actual - forecast_ob$fit)
  naive_error <- 0
  multiplier <- NROW(forecast_ob)/NROW(forecast_ob-1)
  for(i in 2:NROW(forecast_ob)){
    naive_difference <- abs (coredata(forecast_ob$actual[i]) - coredata(forecast_ob$actual[i-1]))
    naive_error <- naive_error + naive_difference
  }
  metric_val <- sum(numerator)/(multiplier*naive_error)
  return(metric_val)
}