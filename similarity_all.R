############################ 1. library ################################


############################ 2. workspace and functions ##############################
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script.dir) # workspace 
loc<-".\\Model_Input\\"
print(script.dir)

#Sys.setlocale('LC_ALL','German')  #active this whe you have "Error in gsub("\"", qstring, col.names, fixed = TRUE)"
################## 4. Read in parameters ##########################
para<-read.csv(paste(script.dir, "Model_Param/parameters.csv", sep = "/"))   #read in parameter

if (ncol(para)!=4 |nrow(para)!=1 ){
  message("Please check Parameter input")} #check data structure

################## 5. Read in data set ##########################
##Read-in files
data_peak<-read.csv(paste(script.dir, "Model_Input/data_peak.csv", sep = "/"))   #read in parameter
# Sorting on reference column
data_peak_sorted <- data_peak[order(data_peak$reference), ]
#count ref samples
n_ref <- sum(data_peak_sorted$reference == 1, na.rm = TRUE) 
#count test samples
n_test<-nrow(data_peak_sorted)-n_ref 

################## 8. data pre-treat ##########################   
#remove first column
data_peak_data<-data_peak_sorted[,-(1:3)] 
#Median
data.median<-apply(data_peak_data[1:n_ref,],2,median) #can be changed to mean if needed
data_peak_data<-rbind(data.median,data_peak_data)
#filter small peaks
if(length(which(data.median<=para$area_median_filter))>0)
{data_peak_data_filter<-data_peak_data[,-which(data.median<=para$area_median_filter)]}else{
  data_peak_data_filter<-data_peak_data}

#numeric the dataframe
data_peak_data_filter<-as.data.frame(lapply(data_peak_data_filter, as.numeric)) 

n_col<-ncol(data_peak_data_filter) #column number
n_row<-nrow(data_peak_data_filter) #row number

################## 9. similarity ##########################    
################## 9.1. calculation cosine similarity and export########################## 
cc<-matrix(rep(0,n_row*n_row),n_row,n_row) # create a matrix with "0"s
log.data_peak_data_filter<-data_peak_data_filter
#log transfer needed, replace with "log.data_peak_data_filter<-log10(data_peak_data_filter)"  and unmute next line
#log.data_peak_data_filter[log.data_peak_data_filter=="-Inf"]<-0 #replace -infs

#n.row=1
for (j in 1:n_row)  #
  cc[1,j] = sum(t(log.data_peak_data_filter[1,])*log.data_peak_data_filter[j,])/sqrt((sum(log.data_peak_data_filter[1,]^2))*sum(log.data_peak_data_filter[j,]^2))

# Calculate the correlation for each pair of rows
for (j in 1:n_row) {
  cc[2, j] <- cor(as.numeric(log.data_peak_data_filter[1, ]), as.numeric(log.data_peak_data_filter[j, ]), method = "pearson")
}

#get result
cc.frame<-data.frame(cc[1:2,])
################## 9.2. calculation Euclidean distance ########################## 

Euclidean_distance<-vector()
for (j in 1:n_row) {
  # Calculate the Euclidean distance between row 1 and row j
  Euclidean_distance[j] <- sqrt(sum((data_peak_data_filter[j, ] - data_peak_data_filter[1, ])^2))
  }
cc.frame<-rbind(cc.frame,Euclidean_distance)
################## 9.3. calculation extent similarity ########################## 

Extent_similarity<-vector()
#n.row=1
for (j in 1:n_row){  #
  Extent_similarity[j] = 1-sum(sqrt((1-data_peak_data_filter[j,]/data_peak_data_filter[1,])^2))/n_col
    }
cc.frame<-rbind(cc.frame,Extent_similarity)
rownames(cc.frame)<-c("Cosine","Pearson","Euclidean_distance","Extent_similarity") #rename columns
colnames(cc.frame)<-c("median",data_peak_sorted$name)
################## 10. export########################## 
write.csv(cc.frame, file = "./Model_Output/similarity.csv") #export similarity
