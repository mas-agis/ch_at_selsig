##Summarize Ts/Tv results
#set main directory containing each chromosome folder for ts_tv ratio
setwd("D:/maulana/third_project/metadata")
#crating an empty dataframe for holding up all ts-tv values
combined = data.frame()
#looping over folder 1 to 29
for (i in 1:29) {
  #set the current directory to specific folder "i"
  setwd(paste0("D:/maulana/third_project/metadata/",i))
  #list files in the current folder
  files = list.files()
  #iterate over files in the current folder
  for (file_name in files) {
    #read each file to temp dataframe and skip the first 7 lines to get only ts and tv values
    temp = read.table(file_name, skip = 7)  
    #added chr number as a column
    temp$chr = i
    #added substring of file name containing breed name
    temp$breed = substr(file_name, 22, 50)  
    #row bind to the combined folder
    combined = rbind(combined, temp)
    #remove the temporary dataframe
    rm(temp)
  }
  
}

#sum up total ts and tv for each breed
library(dplyr)
#get breed names (substring of filename)
breed_name = unique(combined$breed)
#get number of breeds
length(breed_name)
#create empty dataframe to hold values
summary_ts_tv = data.frame()
#looping over each breed
for (name in breed_name) {
  #filter combined dataframe for each breed
  temp = filter(combined, breed == name)
  #grouping and sum the value for ts and tv to tempor dataframe
  tempor = temp  %>% group_by(V1) %>% summarise(total = sum(V2))
  #add breed name to tempor dataframe
  tempor$breed = name
  #row bind the dataframe
  summary_ts_tv = rbind(summary_ts_tv, tempor)
  #remove the temporary dataframe
  rm(temp)
  rm(tempor)
}

#write the summary_ts_tv as csv
write.csv(summary_ts_tv, "D:/maulana/third_project/metadata/summary_ts_tv.csv", row.names = FALSE)
