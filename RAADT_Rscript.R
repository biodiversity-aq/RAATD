#==============================================================
#
# biodiversity.aq data processing protocol for the RAATD data
#
# 2020-08-03
#==============================================================

# Setting up
#---------------------------------------------
#libraries
library(readxl)
library(stringr)
library(worrms)
library(tidyr)
library(dplyr)
library(mapview)
library(Orcs)
library(sp)
library(rgeos)

# file directory
loc <- "/Users/msweetlove/Royal Belgian Institute of Natural Sciences/Royal Belgian Institute of Natural Sciences/Anton Van de Putte - data_processing/03_processed/RAATD"
loc_files <- list.files(paste(loc, "02_interim", sep="/"));loc_files

# input file name :
fl <- "RAATD_metadata.csv"
datafl <- setdiff(loc_files, fl)

# read data
coreData <- read.csv(paste(loc, "02_interim", fl, sep="/"));head(coreData)

# Event Core
#---------------------------------------------
# https://tools.gbif.org/dwca-validator/extension.do?id=dwc:Event
head(coreData)
colnames(coreData)

### example just for PANGEA records
#coreData<-coreData[1:29,]

eventCore <- data.frame(eventID = paste("capture", coreData$individual_id,sep=":"),
                        parentEventID = "",
                        eventRemarks = "",
                        datasetID = coreData$dataset_identifier,
                        samplingProtocol = paste("capture for bio-logging", 
                                                 paste("device type= ", coreData$device_type, sep=""), 
                                                 paste("device id= ", coreData$device_id, sep=""),
                                                 sep=":"),
                        year = coreData$deployment_year,
                        month = coreData$deployment_month,
                        day = coreData$deployment_day,
                        eventDate = "", 
                        eventTime = coreData$deployment_time,
                        startDayOfYear = "",
                        endDayOfYear = "",
                        locality = coreData$deployment_site,
                        decimalLatitude=coreData$deployment_decimal_latitude,
                        decimalLongitude=coreData$deployment_decimal_longitude,
                        geodeticDatum = rep("WGS84", nrow(coreData)),
                        georeferenceProtocol = rep("ArgosLocation (http://www.argos-system.org/manual/3-location/34_location_classes.htm)", nrow(coreData)),
                        footprintWKT = "",
                        license = rep("CC BY 4.0", nrow(coreData)),
                        stringsAsFactors = FALSE)


#eventDate  as [YYYY]-[MM]-[DD]T[hh]:[mm]
eventCore$eventDate<-paste(eventCore$year, eventCore$month, eventCore$day, sep="-")
eventCore$eventDate<-gsub("-$", "", eventCore$eventDate)
eventCore$eventDate<-gsub("-$", "", eventCore$eventDate)

#eventRemarks to store keep_or_not info, comments, data contacts and contact emails
for(i in 1:nrow(coreData)){
  msg <- paste("individual_id= ", coreData[i,]$individual_id, sep="")
  if(tolower(coreData[i,]$keepornot)=="discard"){
    msg <- paste(msg, " |analysis remark= not used in analyses", sep="")
  }
  if(!is.na(coreData[i,]$comments) &&
     tolower(coreData[i,]$comments)!="" &
     coreData[i,]$comments != "NA"){
    msg <- paste(msg, " |comments= ", coreData[i,]$comments, sep="")
  }
  if(!is.na(coreData[i,]$data_contact) &&
     (coreData[i,]$data_contact)!="" &
     coreData[i,]$data_contact != "NA"){
    msg <- paste(msg, " |data_contact= ", coreData[i,]$data_contact, sep="")
  }
  if(!is.na(coreData[i,]$contact_email) &&
     tolower(coreData[i,]$contact_email)!="" &
     coreData[i,]$contact_email != "NA"){
    msg <- paste(msg, " |contact_email= ", coreData[i,]$contact_email, sep="")
  }
  eventCore[i,]$eventRemarks <- msg
}




# occurrence extension
#---------------------------------------------
# prep coredata
coreData$eventID <- eventCore$eventID

# make an occurerence extention for each file, pool at the end
occID_max <- 0
occurrenceExtension_complete <- data.frame()

# make a dictionary of the taxa and their urnID
taxid_key <- data.frame(taxname = unique(coreData$scientific_name), scientificNameID=NA)
for(nc in 1:nrow(taxid_key)){
  taxon<-as.character(taxid_key[nc,]$taxname)
  if(!taxon %in% c("", "NA", NA)){
    taxid <- tryCatch({
      tx <- worrms::wm_name2id(taxon)
    }, error = function(e){
      tx <- ""
    }
    ) 
    if(taxid != ""){
      taxid<-paste("urn:lsid:marinespecies.org:taxname:", taxid, sep="")
    }
    taxid_key[nc,]$scientificNameID <- taxid
  }
}

for(f in 1:length(datafl)){
  print(paste("1.   file", as.character(f)))
  data_f <- read.csv(paste(loc, "02_interim", datafl[f], sep="/"), header=TRUE, stringsAsFactors = FALSE)
  
  data_f <- data_f[data_f$individual_id %in% coreData$individual_id,]

  # dates
  dates <- str_split_fixed(data_f$date, "-| ", 4)
  data_f$year <- dates[,1]
  data_f$month <- dates[,2]
  data_f$day <- dates[,3]
  data_f$time <- dates[,4]
  data_f$date<- gsub(" ", "T", data_f$date)

  print(paste("2.   initiating occurrence extension"))
  
  occurrenceExtension <- data.frame(occurrenceID = paste(data_f$individual_id, "deployment", 
                                                         stringr::str_pad((occID_max+1):(occID_max+nrow(data_f)), 7, pad = "0"), 
                                                         sep=":"),
                                    datasetID = unname(unlist(sapply(data_f$individual_id, FUN = function(x){
                                      gsub(x,coreData[coreData$individual_id==x,]$dataset_identifier,x)}))),
                                    type = rep("deployment", nrow(data_f)),
                                    modified = rep(Sys.Date(), nrow(data_f)),
                                    basisOfRecord = rep("MachineObservation", nrow(data_f)),
                                    eventID = unname(unlist(sapply(data_f$individual_id, FUN = function(x){
                                      gsub(x,coreData[coreData$individual_id==x,]$eventID,x)}))),
                                    year = data_f$year,
                                    month = data_f$month,
                                    day = data_f$day,
                                    eventTime = data_f$time,
                                    eventDate = data_f$date,
                                    decimalLatitude = data_f$decimal_latitude,
                                    decimalLongitude = data_f$decimal_longitude,
                                    #fieldNotes = "",
                                    dynamicProperties = "",
                                    scientificName = unname(unlist(sapply(data_f$individual_id, FUN = function(x){
                                      gsub(x,coreData[coreData$individual_id==x,]$scientific_name,x, fixed=TRUE)}))),
                                    scientificNameID = "", 
                                    vernacularName = unname(unlist(sapply(data_f$individual_id, FUN = function(x){
                                      gsub(x,coreData[coreData$individual_id==x,]$common_name,x)}))),
                                    lifeStage = unname(unlist(sapply(data_f$individual_id, FUN = function(x){
                                      gsub(x,coreData[coreData$individual_id==x,]$age_class,x)}))),
                                    sex = unname(unlist(sapply(data_f$individual_id, FUN = function(x){
                                      gsub(x,coreData[coreData$individual_id==x,]$sex,x)}))),
                                    stringsAsFactors = FALSE)
  
  occID_max <- occID_max + nrow(data_f)
  
  #dynamicProperties
  print(paste("3.   dynamic properties"))
  msg <- rep("", nrow(data_f))
  idx <- !is.na(data_f$breeding_stage) & nzchar(data_f$breeding_stage) & data_f$breeding_stage != "NA"
  msg[idx] <- paste0("{\"breeding_stage\": ", data_f$breeding_stage[idx], "}")
  idx <- !is.na(data_f$location_quality) & nzchar(data_f$location_quality) & data_f$location_quality != "NA"
  msg[idx] <- paste0(msg[idx], paste0(", {\"location_quality\": ", data_f$location_quality[idx], "}"))
  msg <- gsub("^, ", "", msg)
  occurrenceExtension$dynamicProperties <- msg
  
  #scientificNameID
  print(paste("4.   scientificNameID"))
  occurrenceExtension$scientificNameID <- unname(unlist(sapply(as.character(occurrenceExtension$scientificName), FUN = function(x){
    gsub(x,taxid_key[taxid_key$taxname==x,]$scientificNameID,x)})))
  
  #eventDate  as [YYYY]-[MM]-[DD]T[hh]:[mm]
  #print(paste("5.   eventDate"))
  #occurrenceExtension$eventDate<-paste(occurrenceExtension$year, occurrenceExtension$month, occurrenceExtension$day, sep="-")
  #occurrenceExtension$eventDate<-gsub("-$", "", occurrenceExtension$eventDate)
  #occurrenceExtension$eventDate<-gsub("-$", "", occurrenceExtension$eventDate)
  
  occurrenceExtension_complete <- dplyr::bind_rows(occurrenceExtension_complete, occurrenceExtension)
}

occurrenceExtension <- occurrenceExtension_complete

# Back to the eventCore
#---------------------------------------------                                
# footprintWKT in eventCore
# POLYGON of track
for(i in 1:nrow(eventCore)){
  print(i)
  coords <- data.frame(y=occurrenceExtension[occurrenceExtension$eventID==eventCore[i,]$eventID,]$decimalLatitude,
                       x=occurrenceExtension[occurrenceExtension$eventID==eventCore[i,]$eventID,]$decimalLongitude,
                       year=occurrenceExtension[occurrenceExtension$eventID==eventCore[i,]$eventID,]$year,
                       month=occurrenceExtension[occurrenceExtension$eventID==eventCore[i,]$eventID,]$month,
                       day=occurrenceExtension[occurrenceExtension$eventID==eventCore[i,]$eventID,]$day)
  if(nrow(coords)>0){
    coords$date <- paste(coords$year, coords$month, coords$day, sep="-")
    keep<-c()
    pastdates<-c()
    for(cx in 1:nrow(coords)){
      if(!coords[cx,]$date %in% pastdates){
        keep<-c(keep, cx)
        pastdates <- coords[cx,]$date
      }
    }
    coords <- coords[keep,]
    dates <- coords[order(coords[,"year"], coords[,"month"], coords[,"day"]),]
    mindate <- dates[1,]$date
    maxdate <- dates[nrow(dates),]$date
    coords <- as.matrix(coords[,c("x", "y")])
    polygx <- Orcs::coords2Polygons(coords, hole=FALSE, ID=eventCore[i,]$eventID)
    eventCore[i,]$footprintWKT <- writeWKT(polygx)
    eventCore[i,]$startDayOfYear <- mindate
    eventCore[i,]$endDayOfYear <- maxdate
  }
}

# Quality Check on geographic points
# plot points on a map
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("coordinate QC") +
  geom_point(data = occurrenceExtension, aes(y=decimalLatitude, x=decimalLongitude),
             colour="red",size=0.5)

occurrences$footprintwkt <-  paste(rep("POINT(", nrow(occurrenceExtension)),
                                   occurrenceExtension$decimalLongitude,
                                   rep(",", nrow(occurrenceExtension)),
                                   occurrenceExtension$decimalLatitude,
                                   rep(")",nrow(occurrenceExtension)), sep="")


### additional issues (code added on later)
#---------------------------------------------
events<-eventCore
occurrences<-occurrenceExtension

# issue with eventIDs
# standardize names format
setdiff(unique(occurrences$eventID), unique(events$eventID))
occurrences$eventID<-unname(unlist(sapply(occurrences$eventID, FUN = function(x){
  gsub("14030908 (47)_43922_2004","capture:14030908 (47)_43922_2004",x, fixed=TRUE)})))
occurrences$eventID<-unname(unlist(sapply(occurrences$eventID, FUN = function(x){
  gsub("14030935 (25)_43921_2003","capture:14030935 (25)_43921_2003",x, fixed=TRUE)})))
occurrences$eventID<-unname(unlist(sapply(occurrences$eventID, FUN = function(x){
  gsub("14052370 (152)_55166_2004","capture:14052370 (152)_55166_2004",x, fixed=TRUE)})))
occurrences$eventID<-unname(unlist(sapply(occurrences$eventID, FUN = function(x){
  gsub("14052638 (154)_55167_2004","capture:14052638 (154)_55167_2004",x, fixed=TRUE)})))
occurrences$eventID<-unname(unlist(sapply(occurrences$eventID, FUN = function(x){
  gsub("120425??_20875_2001","capture:120425??_20875_2001",x, fixed=TRUE)})))

# issue eventTime
# all times as UCT
occurrences$eventTime<-unname(unlist(sapply(as.character(occurrences$eventTime), FUN = function(x){
  if(lengths(regmatches(x, gregexpr(":", x)))==1){
    x<-paste(x, ":00", sep="")
  }else{x<-x};return(x)})))
occurrences$eventTime<-unname(unlist(sapply(as.character(occurrences$eventTime), FUN = function(x){
  if(nchar(strsplit(x, ":")[[1]][1])==1){
    x<-paste("0", x, sep="")
  }else{x<-x};return(x)})))
occurrences$eventTime<-unname(unlist(sapply(as.character(occurrences$eventTime), FUN = function(x){
  if(x!=""){
    x<-paste(x, "Z", sep="")
  }else{x<-x};return(x)})))


events$eventTime<-unname(unlist(sapply(as.character(events$eventTime), FUN = function(x){
  if(lengths(regmatches(x, gregexpr(":", x)))==1){
    x<-paste(x, ":00", sep="")
  }else{x<-x};return(x)})))
events$eventTime<-unname(unlist(sapply(as.character(events$eventTime), FUN = function(x){
  if(x!="" && nchar(strsplit(x, ":")[[1]][1])==1){
    x<-paste("0", x, sep="")
  }else{x<-x};return(x)})))
events$eventTime<-unname(unlist(sapply(as.character(events$eventTime), FUN = function(x){
  if(x!=""){
    x<-paste(x, "Z", sep="")
  }else{x<-x};return(x)})))

# issue incorrect interpretation of startDayOfYear and endDayOfYear
events$eventDate<-paste(events$startDayOfYear, events$endDayOfYear, sep="/")
events<-events[,!colnames(events) %in% c("startDayOfYear", "endDayOfYear")]

# issue georeferenceRemarks
events$georeferenceRemarks <- "place of animal capture"

# issue HumanObservation
occ_events<-occurrences[!duplicated(occurrences$eventID),]
occ_events$basisOfRecord<-"humanObservation"
occ_events$dynamicProperties <- ""
occ_events$occurrenceID <- unname(unlist(sapply(occ_events$occurrenceID, FUN = function(x){
  gsub("deployment","capture",x, fixed=TRUE)})))

for(i in 1:nrow(occ_events)){
  xx<-strsplit(occ_events[i,]$occurrenceID, "capture:")[[1]]
  xx<-xx[length(xx)-1]
  occ_events[i,]$occurrenceID<-paste(xx, ":capture:", 
                                     paste(rep(0, 4-nchar(i)), collapse=""), 
                                     i, sep="")
}

#library(taRifx)
#occ_events<-remove.factors(occ_events)

for(i in 1:nrow(occ_events)){
  yy <- as.character(events[events$eventID==occ_events[i,]$eventID, ]$year)
  mm <- as.character(events[events$eventID==occ_events[i,]$eventID, ]$month)
  dd <- as.character(events[events$eventID==occ_events[i,]$eventID, ]$day)
  
  occ_events[i,]$year <- yy
  occ_events[i,]$month <- mm
  occ_events[i,]$day <- dd
  occ_events[i,]$eventDate <- paste(yy, mm, dd, sep="-")
  occ_events[i,]$eventTime <- as.character(events[events$eventID==occ_events[i,]$eventID, ]$eventTime)
  
  occ_events[i,]$decimalLatitude <- as.character(events[events$eventID==occ_events[i,]$eventID, ]$decimalLatitude)
  occ_events[i,]$decimalLongitude <- as.character(events[events$eventID==occ_events[i,]$eventID, ]$decimalLongitude)
}

occurrences<-rbind(occ_events, occurrences)

# issue type
# remove the type col
occurrences<-occurrences[,!colnames(occurrences) %in% "type"]

# issue location_quality
msg <- rep("", nrow(occurrences))
idx <- grepl("location_quality", occurrences$dynamicProperties, fixed=TRUE)

msg[idx] <- sub(".*?location_quality\": (.).*", "\\1", occurrences$dynamicProperties[idx],perl=TRUE)
occurrences$georeferenceVerificationStatus <- msg


#save data to processed folder
#eventCore
#occurrenceExtension
write.csv(events, paste(loc, "03_processed","RAATD_eventCore_completeDataset__5.csv", sep="/"), na="", row.names = FALSE)
write.csv(occurrences, paste(loc, "03_processed","RAATD_occurrenceExtension_completeDataset__5.csv", sep="/"), na="", row.names = FALSE)



