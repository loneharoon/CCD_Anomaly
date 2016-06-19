
rm(list=ls())
library(data.table)
library(xts)
library(fpc)
library(clv)
library(TSclust)
library(ggplot2)
library(RColorBrewer)# to increase no. of colors
library(forecast)
#library(reshape2)# for melt and dcast

resample_data<- function(xts_datap,xminutes){
  #downsampled data
  ds_data <- period.apply(xts_datap,INDEX = endpoints(index(xts_datap), on = "minutes", k = xminutes ),FUN= mean)
  # align data to nearest time boundary
  align_data <- align.time(ds_data,xminutes*60) # aligning to x minutes
  # return(ds_data)
  return(align_data)
}

file <- "xxx.csv"
dir <- "path to directory containing csv"
df <-  read.csv(paste0(dir,file),header = TRUE)
cat(colnames(df))
df_xts <- xts(df[,-1],as.POSIXct(strftime(df$localminute,format = "%Y-%m-%d %H:%M:%S")))
df_xts_slice <- df_xts["2014-06-01/2014-08-30"]
seqs <- seq(from= index(first(df_xts_slice)), to = index(tail(df_xts_slice,1)), by ="1 min")
temp <- xts(1:length(seqs),seqs)
if(NROW(df_xts_slice)!=NROW(temp)){
  cat("INTERPOLATION DONE")
  cat("before interpolation", NROW(df_xts_slice))
  df_xts_slice <- na.approx(cbind(temp[,-1],df_xts_slice))
  cat("after interpolation", NROW(df_xts_slice))
}

data_10min <- resample_data(df_xts_slice,10) # sample to 10 mintues


data_visualize <- function(){
  # VISUALIZE SPECiFIC PORTION OF DATA
  #http://novyden.blogspot.in/2013/09/how-to-expand-color-palette-with-ggplot.html
  dframe <- data_10min["2014-08-9"]
  dframe <- data.frame(timeindex=index(dframe),coredata(dframe))
  # dframe$dataid <- NULL ; dframe$air1 <-NULL ; dframe$use<- NULL ; dframe$drye1 <- NULL
  df_long <- reshape2::melt(dframe,id.vars = "timeindex")
  colourCount = length(unique(df_long$variable))
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount) # brewer.pal(8, "Dark2") or brewer.pal(9, "Set1")
  g <- ggplot(df_long,aes(timeindex,value,col=variable,group=variable))
  g <- g + geom_line() + scale_colour_manual(values=getPalette)
  g
}

# INTERESTED PORTION OF DATA
use_data <- data_10min$use
names(use_data) <- "power"
train_data <- use_data["2014-06-01/2014-06-20"]
test_data <-  use_data["2014-06-21/2014-08-31"]

arima_module <-function(){
  arima_result <- executeARIMA_withoutTemperature(train_data,test_data)
  #names(arima_result)<- c("actual","fit","lower","upper")
  dropbox_dir <- "/Volumes/MacintoshHD2/Users/haroonr/Dropbox/"
  arima_rds <-readRDS("1105_arima_june_aug.rds",file = paste0(dropbox_dir,"1105/1105_arima_june_aug.rds"))
  names(arima_rds)<- c("actual","fit","lower","upper")
  arima_endresult <-find_arimaresult_status(arima_rds)
  names(arima_endresult) <- "val"
  arima_endresult[arima_endresult$val==TRUE,]
  plot_arimaresults(arima_result,file)

}

neural_module <- function(){
  neural_result <- executeNeuralnet_withoutTemperature(train_data,test_data)
  names(neural_result)<- c("actual","fit","lower","upper")
  #arima_endresult <-find_arimaresult_status(arima_result)
}

if((.indexday(last(train_data)) -.indexday(first(train_data)))%%2!=0)
  stop("Ensure training days are even in number")
firstob <- xts(coredata(train_data[1,]),index(train_data[1,])-10*60)
train_data <- rbind(firstob,train_data)
result <- compute_hourlyBaseline_end50min(train_data,test_data,file)
# saveRDS(result,file="115_ccd_june_aug.rds")
result_df <- data.frame(time=index(result[[1]]),val = coredata(result[[1]]))
result_df$label <- ifelse(result_df$val==TRUE,result_df$time,NA)
result_df$label <- as.POSIXct(result_df$label,origin = "1970-01-01")
result_df[result_df$val==TRUE,1:2]
ggplot(result_df,aes(time,val,label=label))+ geom_text(check_overlap = TRUE,angle=30)
cat(NROW(result_df))



compute_hourlyBaseline_end50min <- function(train_data,test_data,filename){
  # this function computes anomalies hourly, but decision is taken at xx:00:00 time of every hour
  fname = paste0(strsplit(filename,"[.]")[[1]][1],".pdf")
  pdf(fname,width = 14,height = 10)
  par(mfrow=c(4,3))

  looplength =  length(test_data)
  # looplength = 145 # No. of 10 min intervals in a day + 1
  suppressWarnings(rm(tempframe)) # remove existing copy if any exists
  pastanomaly_vec <- as.numeric()# this contains elements of last interval with anom status as 1
  anomalythreshold_len <- 4
  for(i in 1:looplength){
    a = test_data[i,]
    if(!.indexmin(a)%in% c(50))
    { #collect all the data lying within the intended interval except index of detection point
      # browser()
      if(exists("tempframe")){
        # tempframe : contains all the data collected witin an interval
        # df: contains all the historical data in the same time duaration as that of testing data
        tempframe <- rbind(tempframe,a)
        temp <- train_data[.indexhour(train_data) %in% .indexhour(a) & .indexmin(train_data) %in% .indexmin(a) ]
        df <- rbind(df,temp)
      }else{
        tempframe <- a
        df <- train_data[.indexhour(train_data) %in% .indexhour(a) & .indexmin(train_data) %in% .indexmin(a) ]
        
      }
    }else{
      tempframe <- rbind(tempframe,a) # add last observation of each interval
      temp <- train_data[.indexhour(train_data) %in% .indexhour(a) & .indexmin(train_data) %in% .indexmin(a) ]
      df <- rbind(df,temp)
      df_daywise <- split.xts(df,"days",k=1)
      if(length(unique(lapply(df_daywise,length))) != 1){
        stop("Error: there are missing observations on some days. Please Interpolate!")
        browser()
      }
      # browser()
      #df_daywise <- df_daywise[1:10]
      #browser()
      df_matrix <- do.call (rbind,  lapply(df_daywise,function(x)  xts(t(coredata(x$power)),as.Date(first(index(x)),tz="Asia/Kolkata"))) )
      if(dim(df_matrix)[1]%%2 !=0){
        # this constraint is used to ensure that clustering division into 3 clusters does not result in more deltetion
        cat("Problem at Stage A")
        browser()
        stop("Error: Please choose no. of historical days as even")
      }
      #dist_matrix <- dist(df_matrix, method="euclidean",diag=TRUE,upper=TRUE) # alternative is diss function in another package
      
      
      plotintervaldata <- function(){
        # this code is only used if you want to plot interval specific historical data
        loops <- dim(df_matrix)[1]
        plot(x = c(1:6),y = df_matrix[1,]/1000, ylim = c(0,4),type="l",main = "",ann = FALSE)
        for(i in 2:loops){
          lines(x = c(1:6),y = df_matrix[i,]/1000)
        }
      }
      
      
      df_matrix_subset <- get_selected_observations(df_matrix)
      #browser()
      plotselected_clustereddays <- function(){
        # function used only on need basis
        loops <- dim(df_matrix_subset)[1]
        plot(x = c(1:6),y = df_matrix_subset[1,]/1000, ylim = c(0,4),type="l",main = "",ann = FALSE)
        for(i in 2:loops){
          lines(x = c(1:6),y = df_matrix_subset[i,]/1000)
        }
      }
      
      
      hierchicalClustering <- function(){
        # currently this function is not being used. Instead we use mediod clustering
        clusters <- hclust(dist_matrix)
        clustercut <- cutree(clusters,3)
        #  cbind(df_matrix,clustercut)
        clus_op <- xts(as.vector(clustercut),index(df_matrix))
        df_matrix<- cbind(df_matrix,cluster=clus_op)
        table_count <- as.data.frame(table(df_matrix$cluster))
        group_size <- max(table_count$Freq)
        cluster_no <- as.numeric(table_count[table_count$Freq==group_size,]$Var1)
        db_subset <- subset(df_matrix,cluster==cluster_no)
      }
      
      withoutclustering <- function(){
        # this is used when we don't want any clustering of historical data
        temp1 <- apply(df_matrix_subset,2,mean)
        temp1<- xts(temp1,index(tempframe))
        temp_sd <- apply(df_matrix_subset,2,sd)
        temp_sd<- xts(temp_sd,index(tempframe))
        temp_sd_3 <- temp_sd*3
        temp_sd_2 <- temp_sd*2
        plot(tempframe/1000,ylim=c(0,4),main="")
        lines(temp1/1000,col="blue")
        # # lines(temp2/1000,col="blue")
        lines((temp1+temp_sd_3)/1000,col="red",lt=3)
        lines((temp1-temp_sd_3)/1000,col="red",lt=3)
        lines((temp1+temp_sd_2)/1000,col="black",lt=4)
        lines((temp1-temp_sd_2)/1000,col="black",lt=4)
      }
      
      #browser()
      
      baseline_val <- apply(df_matrix_subset,2,mean)
      temp1<- xts(baseline_val,index(tempframe)) # temp1 represents baseline value
      temp_sd <- apply(df_matrix_subset,2,sd)
      temp_sd<- xts(temp_sd,index(tempframe))
      temp_sd_3 <- temp_sd*3
      temp_sd_2 <- temp_sd*2
      plot(tempframe/1000,ylim=c(0,4),main="")
      lines(temp1/1000,col="blue")
      # # lines(temp2/1000,col="blue")
      lines((temp1+temp_sd_3)/1000,col="red",lt=3)
      lines((temp1-temp_sd_3)/1000,col="red",lt=3)
      lines((temp1+temp_sd_2)/1000,col="black",lt=4)
      lines((temp1-temp_sd_2)/1000,col="black",lt=4)
      decision_vector <- c(ifelse((coredata(tempframe) > coredata(temp1+temp_sd_2)),1,0 ))
      decision_vector <- c(pastanomaly_vec,decision_vector) # May be in the past interval there were some anomaly continuing
      #  decision_vector <- c(ifelse((coredata(tempframe) > coredata(temp1+temp_sd_2)) | (coredata(tempframe) < coredata(temp1-temp_sd_2)),1,0 ))
      length_list <- rle(decision_vector) # run length encoding
      is_anom <- any(length_list$lengths >= anomalythreshold_len & length_list$values == 1) # http://stackoverflow.com/a/37571854/3317829
      ## this part of code saves the status of current interval instances with anomaly
      cnt <- 0; lcs <-0 #last common sequence with status as 1
      last <- decision_vector[length(decision_vector)-cnt]
      while(last==1){
        lcs <- lcs + 1
        cnt <- cnt +1
        if(lcs >= anomalythreshold_len)
          break
        last <- decision_vector[length(decision_vector)-cnt]
      }
      pastanomaly_vec <- rep(1,lcs)
      if(length(pastanomaly_vec) >= anomalythreshold_len)
        pastanomaly_vec <- as.numeric() # this sequence is already detected in current interval no need to carry it forward
      
      ## this part ends here
      if(!exists("info_dframe")){
        info_dframe <- xts(is_anom,index(a))
        res_object <-xts(data.frame(actual=coredata(tempframe),fit=coredata(temp1),upper=coredata(temp1+temp_sd_2),lower=coredata(temp1-temp_sd_2)),index(tempframe))
        
      }else{
        temp <- xts(is_anom,index(a))
        info_dframe <- rbind(info_dframe,temp)
        temp2 <-  xts(data.frame(actual=coredata(tempframe),fit=coredata(temp1),upper=coredata(temp1+temp_sd_2),lower=coredata(temp1-temp_sd_2)),index(tempframe))
        res_object <- rbind(res_object,temp2)
        
      }
      # browser()
      if(!exists("dayframe")){
        # dayframe  contains all readings of a day
        dayframe <- tempframe
        if(exists("first_zerobservation")) {
          dayframe <- rbind(first_zerobservation,dayframe)
          rm(first_zerobservation)
        }
      }else{
        dayframe <- rbind(dayframe,tempframe)
      }
      if(.indexhour(a) %in% c(23) & .indexmin(a) %in% c(50)){
        # At the end of day, training data is updated,i.e., window is moved ahead
        cat("\nEntered day end phase");print(index(a))
        #browser()
        train_data <- rbind(train_data,dayframe)
        start_window <- index(first(train_data)) + 24 * 3600
        #end_window <- index(last(train_data))
        train_data <- train_data[index(train_data) >= start_window]
        #train_data <- window(train_data,first= start_window,end=end_window)
        rm(dayframe)
      }
      rm(tempframe)
    }
  }
  dev.off()
  #par(mfrow=c(1,1))
  return(list(statusdf=info_dframe,resultdf=res_object))
}

get_selected_observations  <- function(df_matrix){
  # this function takes input matrix of energy conusmption of n days and outputs selected days for
  # computing prediction intervals and mean/median consumption
  dist_matrix <- dist(df_matrix, method="euclidean",diag=TRUE,upper=TRUE)
  #dist_matrix <- diss(df_matrix,METHOD = "DTWARP")
  dismat <- dist_matrix # reduntant but ok
  pamk.result<-pamk(dismat,krange=1:3,criterion="ch",usepam=TRUE,diss=TRUE)
  nc <-  pamk.result$nc # no. of clusters
  #nc <- 3 # hardcode
  cat("clusters", nc)
  if(nc == 3){
    kmm <- pam(dismat,nc,diss=TRUE,trace.lev = 0) # K-Mediods
    #clusplot(kmm,main="",sub="",labels=2,xlab="",ylab="",xaxt="n",yaxt="n",cex=0.7,span=FALSE)
    db_subset <- get_db_subset3(df_matrix,kmm,nc)
    db_subset <- db_subset[,-dim(db_subset)[2]]
  }else if(nc == 2){
    kmm <- pam(dismat,nc,diss=TRUE,trace.lev = 0)
    #clusplot(kmm,main="",sub="",labels=2,xlab="",ylab="",xaxt="n",yaxt="n",cex=0.7,span=FALSE)
    db_subset <- get_db_subset2(df_matrix,kmm,nc)
    db_subset <- db_subset[,-dim(db_subset)[2]]
  }else{
    cat("\nONLY ONE CLUSTER")
    #stop("ERROR: this case has not been handled")
    # n is one and we don't cluster and take all points as our subset
    db_subset <- df_matrix
  }
  return(db_subset)
}

get_db_subset3 <- function(df_matrix,kmm,nc) {
  # When there are three clusters in data
  # IGNORE CODE REDUNDANCY
  dist_matrix <- dist(df_matrix, method="euclidean",diag=TRUE,upper=TRUE)
  #dist_matrix <- diss(df_matrix,METHOD = "DTWARP")
  dismat <- dist_matrix
  clus_op <- xts(as.vector(kmm$clustering),index(df_matrix)) # create compatible xts for next step
  df_matrix<- cbind(df_matrix,cluster=clus_op) # cbind cluster labels to data
  clus_count <- as.data.frame(table(df_matrix$cluster)) # count # of elements in each cluster
  unique_clus <- length(unique(clus_count$Freq)) #no. of unique clusters via length
  if (nc == unique_clus){ #case <4,3,2> when all are unique; removes last one
    cat("\n Case: 3D")
    group_size <- min(clus_count$Freq)
    del_cluster <- as.numeric(clus_count[clus_count$Freq==group_size,]$Var1)
    cat("\n Removed cluster",del_cluster)
    db_subset <- subset(df_matrix,cluster != del_cluster)
    return(db_subset)
  }else{
    min_ele <- min(clus_count$Freq)
    freq_minelement <- length(which(clus_count$Freq==min_ele))
    if(freq_minelement == 1){ #case <4,4,2>when there is a small cluster; removes last one
      cat("\n Case: 3B")
      group_size <- min(clus_count$Freq)
      del_cluster <- as.numeric(clus_count[clus_count$Freq==group_size,]$Var1)
      cat("\n Removed cluster",del_cluster)
      db_subset <- subset(df_matrix,cluster != del_cluster)
      return(db_subset)
    }else{  #case <4,3,3> or <8,1,1>
      #cls.scatt.data(data, clust, dist="euclidean")
      if(min_ele == 1){ # case <8,1,1>
        # remove  both of clusters containing single element
        cat("\n Case: 3C")
        group_size <- max(clus_count$Freq)
        retain_cluster <- as.numeric(clus_count[clus_count$Freq==group_size,]$Var1)
        cat("\n Retained cluster",retain_cluster)
        db_subset <- subset(df_matrix,cluster == retain_cluster)
        return(db_subset)
      }else{ #case <4,3,3>
        cat("\n Case: 3A")
        cls_ob <- cls.scatt.diss.mx(as.matrix(dismat), kmm$clustering)
        interdist_mat <- cls_ob$intercls.single
        clusters <- c(1,2,3) # names of all clusters
        group_size <- max(clus_count$Freq) # find index of max cluster in next step
        cluster_no <- as.numeric(clus_count[clus_count$Freq==group_size,]$Var1)
        equal_clusters <- setdiff(clusters,cluster_no) # find indexes of equal clusters
        dist_xy1 <- interdist_mat[cluster_no,equal_clusters[1]]
        dist_xy2 <- interdist_mat[cluster_no,equal_clusters[2]]
        del_cluster <- ifelse(dist_xy1 >= dist_xy2,equal_clusters[1],equal_clusters[2])
        cat("\n Removed cluster",del_cluster)
        db_subset <- subset(df_matrix,cluster != del_cluster)
        return(db_subset)
      }
    }
  }
}

get_db_subset2 <- function(df_matrix,kmm,nc) {
  ## When there are two clusters in data
  ## PLEASE IGNORE CODE REDUNDANCY
  cat ("\nCASE: 2 CLUSTERS")
  #dist_matrix <- diss(df_matrix,METHOD = "DTWARP")
  #dist_matrix <- dist(df_matrix, method="euclidean",diag=TRUE,upper=TRUE)
  #dismat <- dist_matrix
  clus_op <- xts(as.vector(kmm$clustering),index(df_matrix))
  df_matrix <- cbind(df_matrix,cluster=clus_op)
  clus_count <- as.data.frame(table(df_matrix$cluster))
  unique_clus <- length(unique(clus_count$Freq)) #no. of unique clusters via length
  if(clus_count$Freq[1] == clus_count$Freq[2]){ # case <5,5>
    # remove none of the elements
    cat("\nCASE: 2A")
    db_subset <- df_matrix
    return(db_subset)
  }else{
    group_size <- min(clus_count$Freq) # find cluster with min. no. of elem.
    if (group_size < (0.25 * dim(df_matrix)[1])){
      # delete cluster containing less than 1/4 of total observations
      cat("\nCASE: 2C")
      del_cluster <- as.numeric(clus_count[clus_count$Freq==group_size,]$Var1)
      cat("\nRemoved cluster",del_cluster)
      db_subset <- subset(df_matrix,cluster != del_cluster)
      return(db_subset)
    }else{
      cat("\nCASE: 2B")
      db_subset <- df_matrix
      return(db_subset)
    }
  }
}

plot_intervalConsumptionAcrosDays <- function(){
  # this function helps to visualize power consumption vertically,i.e, across days
  # df contains power consumption of n days
  
  df_daywise <- split.xts(dfs,"days",k=1)
  daydat <- df_daywise[1:24]
  suppressWarnings(rm(cont_dat))
  for(i in 1:length(daydat)){
    if(i==1){
      cont_dat <- daydat[[i]]
    }else{
      cont_dat <- rbind(cont_dat,daydat[[i]])
    }}
  cont_dat$day<-as.factor(strftime(index(cont_dat),"%m-%d"))# for day wise visualization
  cont_dat$xaxis <- as.numeric(strftime(index(cont_dat),"%M")) # for axis
  cont_dat$xaxis <- ifelse(cont_dat$xaxis==0,60,cont_dat$xaxis)# converting 0th minute to 60 for ploting purposes
  long_df<-data.frame(timestamp= index(cont_dat), coredata(cont_dat))
  ggplot(long_df,aes(xaxis,power/1000)) + facet_grid(day~.) + geom_line()+ labs(x="Time(Minutes)",y= "Power")
}
compute_hourlyBaseline_end00min <- function(train_data,test_data){
  # this function computes anomalies hourly, but decision is taken at xx:50:00 time of every hour
  par(mfrow=c(4,3))
  #pdf("subplot.pdf")
  if(.indexmin(test_data[1,]) %in% 0)
  { # I want to to take decision at 30 min  and 60 min. So for each range I will decide based on readings
    #  c(10,20,30) and c(40,50,0) respectively
    # cat("Starts with 0 min. which I don't like")
    first_zerobservation <- test_data[1,]
    test_data <- test_data[-1,] # remove this observation
  }
  looplength =  length(test_data)
  # looplength = 145 # No. of 10 min intervals in a day + 1
  suppressWarnings(rm(tempframe)) # remove existing copy if any exists
  
  for(i in 1:looplength){
    a = test_data[i,]
    if(!.indexmin(a)%in% c(0))
    { #collect all the data lying within the intended interval except index of detection point
      # browser()
      if(exists("tempframe")){
        # tempframe : contains all the data collected witin an interval
        # df: contains all the historical data in the same time duaration as that of testing data
        tempframe <- rbind(tempframe,a)
        temp <- train_data[.indexhour(train_data) %in% .indexhour(a) & .indexmin(train_data) %in% .indexmin(a) ]
        df <- rbind(df,temp)
      }else{
        tempframe <- a
        df <- train_data[.indexhour(train_data) %in% .indexhour(a) & .indexmin(train_data) %in% .indexmin(a) ]
        
      }
    }else{
      tempframe <- rbind(tempframe,a) # add last observation of each interval
      temp <- train_data[.indexhour(train_data) %in% .indexhour(a) & .indexmin(train_data) %in% .indexmin(a) ]
      df <- rbind(df,temp)
      df_daywise <- split.xts(df,"days",k=1)
      if(length(unique(lapply(df_daywise,length))) != 1){
        
        stop("Error: there are missing observations on some days. Please Interpolate!")}
      # browser()
      #df_daywise <- df_daywise[1:10]
      browser()
      df_matrix <- do.call (rbind,  lapply(df_daywise,function(x)  xts(t(coredata(x$power)),as.Date(first(index(x)),tz="Asia/Kolkata"))) )
      if(dim(df_matrix)[1]%%2 !=0){
        
        stop("Error: Please choose no. of historical days as even")
      }
      #dist_matrix <- dist(df_matrix, method="euclidean",diag=TRUE,upper=TRUE) # alternative is diss function in another package
      
      
      plotintervaldata <- function(){
        # this code is only used if you want to plot interval specific historical data
        loops <- dim(df_matrix)[1]
        plot(x = c(1:6),y = df_matrix[1,]/1000, ylim = c(0,4),type="l",main = "",ann = FALSE)
        for(i in 2:loops){
          lines(x = c(1:6),y = df_matrix[i,]/1000)
        }
      }
      
      
      df_matrix_subset <- get_selected_observations(df_matrix)
      #browser()
      plotselected_clustereddays <- function(){
        # function used only on need basis
        loops <- dim(df_matrix_subset)[1]
        plot(x = c(1:6),y = df_matrix_subset[1,]/1000, ylim = c(0,4),type="l",main = "",ann = FALSE)
        for(i in 2:loops){
          lines(x = c(1:6),y = df_matrix_subset[i,]/1000)
        }
      }
      
      
      hierchicalClustering <- function(){
        # currently this function is not being used. Instead we use mediod clustering
        clusters <- hclust(dist_matrix)
        clustercut <- cutree(clusters,3)
        #  cbind(df_matrix,clustercut)
        clus_op <- xts(as.vector(clustercut),index(df_matrix))
        df_matrix<- cbind(df_matrix,cluster=clus_op)
        table_count <- as.data.frame(table(df_matrix$cluster))
        group_size <- max(table_count$Freq)
        cluster_no <- as.numeric(table_count[table_count$Freq==group_size,]$Var1)
        db_subset <- subset(df_matrix,cluster==cluster_no)
      }
      
      withoutclustering <- function(){
        # this is used when we don't want any clustering of historical data
        temp1 <- apply(df_matrix_subset,2,mean)
        temp1<- xts(temp1,index(tempframe))
        temp_sd <- apply(df_matrix_subset,2,sd)
        temp_sd<- xts(temp_sd,index(tempframe))
        temp_sd_3 <- temp_sd*3
        temp_sd_2 <- temp_sd*2
        plot(tempframe/1000,ylim=c(0,4),main="")
        lines(temp1/1000,col="blue")
        # # lines(temp2/1000,col="blue")
        lines((temp1+temp_sd_3)/1000,col="red",lt=3)
        lines((temp1-temp_sd_3)/1000,col="red",lt=3)
        lines((temp1+temp_sd_2)/1000,col="black",lt=4)
        lines((temp1-temp_sd_2)/1000,col="black",lt=4)
      }
      
      baseline_val <- apply(df_matrix_subset,2,mean)
      temp1<- xts(baseline_val,index(tempframe))
      temp_sd <- apply(df_matrix_subset,2,sd)
      temp_sd<- xts(temp_sd,index(tempframe))
      temp_sd_3 <- temp_sd*3
      temp_sd_2 <- temp_sd*2
      plot(tempframe/1000,ylim=c(0,4),main="")
      lines(temp1/1000,col="blue")
      # # lines(temp2/1000,col="blue")
      lines((temp1+temp_sd_3)/1000,col="red",lt=3)
      lines((temp1-temp_sd_3)/1000,col="red",lt=3)
      lines((temp1+temp_sd_2)/1000,col="black",lt=4)
      lines((temp1-temp_sd_2)/1000,col="black",lt=4)
      decision_vector <- c(ifelse((coredata(tempframe) > coredata(temp1+temp_sd_2)) | (coredata(tempframe) < coredata(temp1-temp_sd_3)),1,0 ))
      length_list <- rle(decision_vector) # run length encoding
      is_anom <- any(length_list$lengths >=4 & length_list$values == 1) # http://stackoverflow.com/a/37571854/3317829
      if(!exists("info_dframe")){
        info_dframe <- xts(is_anom,index(a))
      }else{
        temp <- xts(is_anom,index(a))
        info_dframe <- rbind(info_dframe,temp)
      }
      # browser()
      if(!exists("dayframe")){
        # dayframe  contains all readings of a day
        dayframe <- tempframe
        if(exists("first_zerobservation")) {
          dayframe <- rbind(first_zerobservation,dayframe)
          rm(first_zerobservation)
        }
      }else{
        dayframe <- rbind(dayframe,tempframe)
      }
      if(.indexhour(a) %in% c(0) & .indexmin(a) %in% c(0)){
        # At the end of day, training data is updated,i.e., window is moved ahead
        #
        train_data <- rbind(train_data,dayframe)
        start_window <- index(first(train_data)) + 24 * 3600
        #end_window <- index(last(train_data))
        train_data <- train_data[index(train_data) >= start_window]
        #train_data <- window(train_data,first= start_window,end=end_window)
        rm(dayframe)
      }
      rm(tempframe)
    }
  }
  # dev.off()
  par(mfrow=c(1,1))
  return(info_dframe)
}

arima_extra <- function(){
  # STATIONARITY of series####
  #FIND ARIMA ORDER
  plot(use_data,type="l")
  datas <- ts(use_data,frequency=144)
  Acf(diff(datas))
  Pacf(diff(datas))
  
  # FIND DIFFERENCING
  dat=use_data
  tsob <-  ts(coredata(dat$power),frequency=144)
  stlob<- stl(tsob[,1],"periodic",robust = TRUE)
  plot(stlob)
  
  dat2 = diff(use_data)
  dat2[1,] <- dat2[2,]
  tsob <-  ts(coredata(dat2$power),frequency=144)
  stlob<- stl(tsob[,1],"periodic",robust = TRUE)
  plot(stlob)
}

