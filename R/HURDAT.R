#' HURDAT
#'
#' Appends a dataframe with the names of any storms in the HURDAT database with a center of circulation within a specified radius of a location.
#'
#' @param Data A data frame where the first column is a date or date_time object.
#' @param lat.loc Numeric vector of length one specifying the latitude (in degrees)of the location of interest.
#' @param lon.loc Numeric vector of length one specifying the longitude (in degrees) of the location of interest.
#' @param rad Vector of length one, specifying radius (in km) about the location of interest to report storm names.
#' @return The dataframe \code{Data} with an additional column containing the named storms.
#' @export
HURDAT <- function(Data,lat.loc,lon.loc,rad){

  #Reading in HURDAT2 database
  #HURDAT2 = read.table('HURDAT2.txt', sep="\t")

  #Identify rows which describe event
  n.char = apply(HURDAT2,1,nchar)

  #Function to extract storm (1) names and (2) duration
  storm.name.cal = function(x) substring(x,19,28)
  storm.duration.cal = function(x) as.numeric(trimws(substring(x,34,36)))

  #Vector of storm (1) names and (2) duations
  storm.name = apply(data.frame(HURDAT2[which(n.char==37),]),1,storm.name.cal)
  storm.duration = apply(data.frame(HURDAT2[which(n.char==37),]),1,storm.duration.cal)

  #Vector to append HURDAT2 dataset with the names of the storms
  storm.names = rep(storm.name, storm.duration)

  #CombinIng name vecor with HURDAT2 data
  HURDAT2 = data.frame(storm.names, HURDAT2[-which(n.char==37),])

  #Functions to extract (1) latitude and (2) longitude of storm
  lat.cal = function(x) as.numeric(trimws(substring(x[[2]],24,27)))
  lon.cal = function(x) as.numeric(trimws(substring(x[[2]],30,35)))

  #Vector of (1) lat and (2) long of storms
  lat = as.numeric(apply(HURDAT2,1,lat.cal))
  lon = as.numeric(apply(HURDAT2,1,lon.cal))

  #Radius of earth at sea level in kms
  R = 6371

  #Converting degrees to radians and finding differences
  d.phi = abs(lat - lat.loc) * pi/180
  d.lambda = abs(lon - lon.loc) * pi/180
  phi = lat * pi/180
  phi.loc = lat.loc * pi/180

  #'haversine' formula to calculate the great-circle distance between two points (distance in km)
  a = sin(d.phi/2) * sin(d.phi/2)  + cos(phi) * cos(phi.loc) * sin(d.lambda/2) * sin(d.lambda/2)
  c = 2 * asin(sqrt(a))
  d = R * c

  #which events pass within 350km of S20_T
  events = unique(storm.names[which(d<rad)])

  #Dates of Tropical cyclones
  storm.dates = trimws(substring(HURDAT2[which(d<rad),][[2]],1,9))

  #Convert to standard date format
  storm.dates = as.Date(storm.dates, tryFormats ="%Y%m%d")

  #TC dates
  storm.dates.df = data.frame(storm.dates,trimws(storm.names[which(d<rad)]))
  colnames(storm.dates.df) = c("Date","Name")

  #Removing duplicated rows
  storm.dates.df = storm.dates.df[!duplicated(storm.dates.df),]

  #Adding storm names to the dataframe input
  #If hourly data remove HH:MM:SS then add it back in
  if(class(Data[,1])[1]=="POSIXct"){
    date_time  = Data[,1]
    Data[,1]<-as.Date(Data[,1])
    names(Data)[1] = "Date"
    Data.HURDAT = left_join(Data,storm.dates.df,by="Date",multiple="first")
    Data.HURDAT[,1] = date_time
  } else{
    Data.HURDAT = left_join(Data,storm.dates.df,by="Date",multiple="first")
  }


  #Return the output
  return(Data.HURDAT)
}
