# This script will hold several utility functions designed to help analyze the output of VAST models

#' ConstructStreamVAST
#' This function constructs objects with class "StreamVAST". Depending on the arguments supplied, it will either
#' initialize the object from countdata and reachdata, or it will add details to an already defined StreamVAST object.
#'
#' @param countdata A dataframe containing count data, such as outputed from AssembleReddData
#' @param reachdata A sf object with data for reaches, such as from AssignReaches
#' @param surveydata A optional data frame or sf object with survey information, such as from MakeSurveyTracks, doesn't do much right now
#' @param countname A column name from countdata to use
#' @param reachname a column name from reachdata to use
#' @param unitconv a value to divide river lengths
#'
#' @return A streamvast object with appropriate formatting
#' @export
#'
#' @examples
ConstructStreamVAST<-function(countdata,reachdata,surveydata,
                              countname,reachname){

  # if the StreamVAST object is missing, we will construct it from countdata and reachdata
  if(missing(countdata) | missing(reachdata)){
    stop("Arguments 'countdata' and 'reachdata' are required to initialize StreamVAST")
  }
  print("Initializing StreamVAST based on countdata and reachdata")
  if(inherits(countdata,"character")){countdata<-utils::read.csv(countdata)}
  if(inherits(reachdata,"character")){reachdata<-sf::st_read(reachdata)}
  if(inherits(reachdata,"sf")==F){stop("Argument reachdata must be a sf object with LINESTRING geometries")}
  if(missing(surveydata)){surveydata<-NA}
  if(inherits(surveydata,"character")){surveydata<-sf::st_transform(sf::st_read(surveydata),"wgs84")}

  # account for different capitalizations of day
  names(countdata)[names(countdata)%in%c("Day","day")]<-"Day"
  if(inherits(surveydata,c("sf","data.frame"))){names(surveydata)[names(surveydata)%in%c("Day","day")]<-"Day"}

  c.names<-names(countdata)
  r.names<-names(reachdata)

  # will try to identify a name for countname and reachname, if absent
  if(missing(countname)){
    namepicker<-which(c("count","redds","redd","fish","abundance","abund")%in%tolower(c.names))
    if(length(namepicker)==0){stop("Could not identify a column for counts. Please specify 'countname' and try again.")}
    countname<-c.names[tolower(c.names)==c("count","redds","redd","fish","abundance","abund")[namepicker[1]]]
  }
  if(missing(reachname)){
    namepicker<-which(c("reach","reachid","reach_id","reach_no","streamid","stream_id","stream_no")%in%tolower(r.names))
    if(length(namepicker)==0){stop("Could not identify a column for reach. Please specify 'reachname' and try again.")}
    reachname<-r.names[tolower(r.names)==c("reach","reachid","reach_id","reach_no","streamid","stream_id","stream_no")[namepicker[1]]]
  }

  # set up some things, notably the midpoints, which requires some transformations
  reachdata<-sf::st_transform(reachdata,"wgs84")
  reach.box<-sf::st_bbox(reachdata)
  EPSG=32700-((sign(mean(reach.box[c(2,4)]))+1)/2)*100+floor((180+mean(reach.box[c(1,3)]))/6)+1
  reachdata.alt<-sf::st_transform(reachdata,EPSG)
  midpoints<-sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=reachdata.alt,sample = .5),
                                                 crs = "wgs84"))
  midpoints.sf<-sf::st_as_sf(as.data.frame(midpoints[,1:2]),coords=1:2,crs="wgs84")

  # Calculate some values that aren't in the data already and add them
  reachdata$Length<-as.numeric(sf::st_length(reachdata))/1000
  reachdata$CenterLon<-midpoints[,1]
  reachdata$CenterLat<-midpoints[,2]

  countdata$Length<-reachdata$Length[match(countdata[,reachname],as.data.frame(reachdata)[,reachname])]
  countdata$Density<-countdata[,countname]/countdata$Effort


  # add stream distance to reachdata
  reachnet<-suppressWarnings(sfnetworks::as_sfnetwork(reachdata))
  root<-LocateRoot(reachdata)

  dist.net<-sfnetworks::st_network_blend(x = reachnet,y = midpoints.sf)
  dist.index<-sf::st_nearest_feature(x = midpoints.sf, y=sfnetworks::activate(dist.net,"nodes"))
  reachdata$streamdist<-as.numeric(sfnetworks::st_network_cost(x=dist.net,from=root)[dist.index])/1000

  count.hab<-reachdata$habitat[match(countdata[,reachname],as.data.frame(reachdata)[,reachname])]
  if(any(count.hab==0)){
    warning(paste("Removed ",sum(count.hab==0)," count records in non-habitat areas"))
    print("Summary of excluded records")
    print(utils::head(countdata[count.hab==0,]))
  }
  countdata<-countdata[count.hab==1,]


  # diagnostics, make sure the data is valid
  overlapplot<-ggplot2::ggplot()+
    ggplot2::geom_sf(data=reachdata$geometry,lwd=5,ggplot2::aes(color="1"))
  if(any(reachdata$habitat==F)){
    overlapplot<-overlapplot+
      ggplot2::geom_sf(data=reachdata[reachdata$habitat==F,],lwd=4,ggplot2::aes(col="2"))
    legend.labels<-c("Spatial Frame","Nonhabitat","Survey Coverage","Data Coverage")
  }else{
    legend.labels<-c("Spatial Frame","Survey Coverage","Data Coverage")
  }
  if(class(surveydata)[1]=="sf"){
    overlapplot<-overlapplot+ggplot2::geom_sf(data=surveydata$geometry,lwd=2,ggplot2::aes(col="3"))
    bluecol<-"royalblue"
  }else{
    bluecol<-"turquoise2"
  }
  overlapplot<-overlapplot+
    ggplot2::geom_sf(data=reachdata$geometry[as.data.frame(reachdata)[,reachname]%in%countdata$Reach],lwd=1,ggplot2::aes(col="4"))+
    ggplot2::theme_bw()+ggplot2::scale_color_manual(values = c("1"="black","2"="gray60","3"="tomato","4"=bluecol),
                                                    labels=legend.labels)+
    ggplot2::theme(legend.title = ggplot2::element_blank())
  print(overlapplot)

  title<-NULL
  vastdata<-NULL
  timetable<-NULL
  timescale<-NULL
  timerange<-NULL
  vastnetwork<-tinyVAST::sfnetwork_mesh(stream = sfnetworks::as_sfnetwork(reachdata))
  vastnetwork$table$dist<-vastnetwork$table$dist/1000
  vastmodel<-NULL
  stats<-list(AIC=NULL,RMSE=NULL,rhoP=NULL,rhoS=NULL)
  preds<-NULL
  eval<-NULL
  spacedata<-NULL
  timedata<-NULL
  aucdata<-NULL
  escapedata<-NULL
  sims<-list(predsims=NULL,evalsims=NULL,aucsims=NULL,auctotalsims=NULL,escapesims=NULL)
  dharmafit<-NULL

  StreamVAST<-list(countdata=countdata,reachdata=reachdata,surveydata=surveydata,
                   title=NULL,countname=countname,reachname=reachname,vastdata=NULL,
                   timetable=NULL,timescale=NULL,vastnetwork=vastnetwork,vastmodel=NULL,
                   stats=list(AIC=NULL,RMSE=NULL,rhoP=NULL,rhoS=NULL),preds=NULL,
                   eval=NULL,spacedata=NULL,timedata=NULL,aucdata=NULL,dharmafit=NULL)

  class(StreamVAST)<-"streamvast"
  return(StreamVAST)
}


#' This function takes a streamvast object defines a temporal sampling frame and formats everything for use with vast
#'
#' @param streamvast a streamvast object, minimally with countdata, reachdata defined
#' @param startdate dates to truncate or expand the data in yyyy-mm-dd format
#' @param enddate dates to truncate or expand the data in yyyy-mm-dd format
#' @param Time character giving the timescale to use
#' @param season a integer vector of length with order start Month, start day, end month, end day for a 'year' or 'season'
#' @param padyear adds zeroes to the start and end of each runyear
#' @param padhabiat adds zeros at regular intervals to areas designated non-habitat
#'
#' @return A streamvast object with formated with the appropriate temporal frame
#' @export
#'
#' @examples
SetTemporalFrame<-function(streamvast,startdate=NA,enddate=NA,Time="Year",
                           season,padyear=F,padhabitat=F){

  if(length(streamvast$timetable)>0){warning("Temporal sampling frame already defined for this object. Overwriting pre-existing temporal frame.")}

  #Get the dates ; fun note, "ifelse" will screw up dates
  data.dates<-as.Date(paste(streamvast$countdata$Year,streamvast$countdata$Day,sep="-"),format="%Y-%j")
  if(is.na(startdate)){
    startdate<-min(data.dates)
  }else{
    startdate<-as.Date(startdate)
  }
  if(is.na(enddate)){
    enddate<-max(data.dates)
  }else{
    enddate<-as.Date(enddate)
  }

  good.data<-streamvast$countdata[data.dates>=startdate & data.dates<=enddate,]
  good.dates<-data.dates[data.dates>=startdate & data.dates<=enddate]

  startyear<-as.integer(format(startdate,"%Y"))
  endyear<-as.integer(format(enddate,"%Y"))

  # a whole bunch of stuff to allow users to specify non-standard seasons or runyears
  if(missing(season)){
    minmonth<-min(as.integer(format(good.dates,"%m")))
    minday<-1
    maxmonth<-max(as.integer(format(good.dates,"%m")))
    maxday<-c(31,28,31,30,31,30,31,31,30,31,30,31)[maxmonth]
    yearbreak<-c(1,1)
  }else{
    if(inherits(season,c("integer","numeric"))==F){stop("Season argument must be an integer vector")}

    minmonth<-season[1]
    minday<-season[2]
    maxmonth<-season[3]
    maxday<-season[4]

    minjul<-as.integer(format(as.Date(paste(1,minmonth,minday,sep="-")),"%j"))
    yearbreak<-as.integer(format(as.Date(paste0(1,"-",minjul),"%Y-%j"),"%m"))
    yearbreak<-c(yearbreak,as.integer(format(as.Date(paste0(1,"-",minjul),"%Y-%j"),"%d")))
  }

  good.years<-as.integer(format(good.dates,"%Y"))
  good.year.breaks<-as.Date(paste(good.years,yearbreak[1],yearbreak[2],sep="-"))
  runyears<-good.years-as.integer(good.dates<good.year.breaks)

  runyear.vec<-min(runyears):max(runyears)
  runyear.start<-as.Date(paste(runyear.vec,minmonth,minday,sep="-"))
  runyear.end<-as.Date(paste(runyear.vec,maxmonth,maxday,sep="-"))

  runyear.end[runyear.end<=runyear.start]<-as.Date(paste(runyear.vec+1,maxmonth,maxday,sep="-"))
  year.table<-data.frame(Runyear=runyear.vec,Start=runyear.start,End=runyear.end)

  # prune data to remove anything that's out of season
  year.table.match<-match(runyears,year.table$Runyear)
  good.vec2<-good.dates>=year.table$Start[year.table.match] & good.dates<=year.table$End[year.table.match]
  good.data2<-good.data[good.vec2,]
  good.dates2<-good.dates[good.vec2]

  year.table$Start[1]<-startdate
  year.table$End[nrow(year.table)]<-enddate

  # Now need to go through and get tables that convert dates into sequential time steps
  if(Time=="Year"){

    year.table$Date<-as.Date("1-1-1")
    for(i in 1:nrow(year.table)){
      year.table$Date[i]<-mean(c(year.table$Start[i],year.table$End[i]))
    }
    time.table<-data.frame(Runyear=year.table$Runyear,
                           Year=as.integer(format(year.table$Date,"%Y")),
                           year.table[,c(2,4,3)])
  }else{

    ############  This whole section is mostly untested  ####################
    for(i in 1:nrow(year.table)){
      tseq<-seq(from=year.table$Start[i],to=year.table$End[i],by=1)

      if(tolower(Time)=="month"){
        time.table0<-unique(data.frame(Runyear=year.table$Runyear[i],
                                       Years=as.integer(format(tseq,"%Y")),
                                       Months=as.integer(format(tseq,"%m"))))
        Refstart<-as.Date(paste(time.table0$Years,time.table0$Months,1,sep="-"))
        Refdate<-as.Date(paste(time.table0$Years,time.table0$Months,15,sep="-"))

        leaps<-as.integer(format(as.Date(paste0(time.table0$Years,"-12-31")),"%j"))==366
        monthdays<-c(31,28,31,30,31,30,31,31,30,31,30,31)[as.integer(time.table0$Months)]
        monthdays[leaps & as.integer(time.table0$Months)==2]<-29
        Refend<-as.Date(paste(time.table0$Years,as.integer(time.table0$Months),monthdays,sep="-"))

        time.table0<-cbind(time.table0,data.frame(Start=Refstart,Date=Refdate,End=Refend))
      }

      if(tolower(Time)%in%c("biweek","week")){

        tdivide<-ifelse(tolower(Time)=="biweek",14,7)

        Refstart<-tseq[which((1:length(tseq)%%tdivide)==1)]
        Refend<-tseq[which((1:length(tseq)%%14)==0)]
        if(length(Refend)!=length(Refstart)){Refend<-c(Refend,tseq[length(tseq)])}

        Refdate<-rep(as.Date("1-1-1"),length(Refstart))
        for(j in 1:length(Refdate)){
          Refdate[j]<-mean(c(Refstart[j],Refend[j]))
        }

        time.table0<-data.frame(Runyear=year.table$Runyear[i],
                                Year=as.integer(format(Refdate,"%Y")),
                                Start=Refstart,Date=Refdate,End=Refend)
      }
      if(i==1){
        time.table<-time.table0
      }else{
        time.table<-rbind(time.table,time.table0)
      }
    }
  }
  time.table$Time<-1:nrow(time.table)

  # we'll have some formatting steps here
  time.table$Month<-as.integer(format(time.table$Date,format="%m"))
  time.table$Day<-as.integer(format(time.table$Date,format="%j"))
  runstart<-stats::aggregate(time.table$Start,by=list(time.table$Runyear),FUN=min)
  time.table$Statday<-as.integer(ceiling(time.table$Date-runstart$x[match(time.table$Runyear,runstart$Group.1)]))+1

  # Now we set up the data
  timematch<-sapply(good.dates2,FUN=function(x){return(which(time.table$Start<=x & time.table$End>=x))})
  reachmatch<-match(good.data2$Reach,streamvast$reachdata$Reach)

  # Get the preliminary output object that has just the original data
  out.data<-data.frame(Time=time.table$Time[timematch],
                       Runyear=time.table$Runyear[timematch],
                       Date=good.dates2,
                       Year=good.data2$Year,
                       Month=as.integer(format(good.dates2,"%m")),
                       Day=as.integer(format(good.dates2,"%d")),
                       Statday=as.integer(good.dates2-time.table$Start[timematch])+1,
                       Reach=good.data2$Reach,
                       Species=streamvast$countname,
                       Dist="obs",
                       Counts=good.data2[,streamvast$countname],
                       Density=good.data2[,streamvast$countname]/good.data2$Effort,
                       Length=good.data2$Length,
                       Effort=good.data2$Effort,
                       Lat=streamvast$reachdata$CenterLat[reachmatch],
                       Lon=streamvast$reachdata$CenterLon[reachmatch],
                       Original=T)

  time.table$Observed<-sapply(time.table$Time,FUN = function(x){return(any(x==out.data$Time))})

  # Now make padding data, if adding zeros to the ends of each season
  if(padyear){
    nreaches<-sum(streamvast$reachdata$habitat)
    reach.vec<-which(streamvast$reachdata$habitat==1)

    front.data<-data.frame(Time=rep(time.table$Time,each=nreaches),
                           Runyear=rep(time.table$Runyear,each=nreaches),
                           Date=rep(time.table$Start,each=nreaches),
                           Year=rep(as.integer(format(time.table$Start,"%Y")),each=nreaches),
                           Month=rep(as.integer(format(time.table$Start,"%m")),each=nreaches),
                           Day=rep(as.integer(format(time.table$Start,"%d")),each=nreaches),
                           Statday=1,
                           Reach=rep(reach.vec,nrow(time.table)),
                           Species=streamvast$countname,
                           Dist="obs",
                           Counts=0,
                           Density=0,
                           Length=rep(streamvast$reachdata$Length[reach.vec],nrow(time.table)),
                           Effort=rep(streamvast$reachdata$Length[reach.vec],nrow(time.table)),
                           Lat=rep(streamvast$reachdata$CenterLat[reach.vec],nrow(time.table)),
                           Lon=rep(streamvast$reachdata$CenterLon[reach.vec],nrow(time.table)),
                           Original=F)

    back.data<-data.frame(Time=rep(time.table$Time,each=nreaches),
                          Runyear=rep(time.table$Runyear,each=nreaches),
                          Date=rep(time.table$End,each=nreaches),
                          Year=rep(as.integer(format(time.table$End,"%Y")),each=nreaches),
                          Month=rep(as.integer(format(time.table$End,"%m")),each=nreaches),
                          Day=rep(as.integer(format(time.table$End,"%d")),each=nreaches),
                          Statday=rep(as.integer(time.table$End-time.table$Start)+1,each=nreaches),
                          Reach=rep(reach.vec,nrow(time.table)),
                          Species=streamvast$countname,
                          Dist="obs",
                          Counts=0,
                          Density=0,
                          Length=rep(streamvast$reachdata$Length[reach.vec],nrow(time.table)),
                          Effort=rep(streamvast$reachdata$Length[reach.vec],nrow(time.table)),
                          Lat=rep(streamvast$reachdata$CenterLat[reach.vec],nrow(time.table)),
                          Lon=rep(streamvast$reachdata$CenterLon[reach.vec],nrow(time.table)),
                          Original=F)

    out.data<-rbind(out.data,front.data,back.data)
  }
  # If desired, zero out nonhabitat areas
  if(padhabitat & any(streamvast$reachdata$habitat==0)){
    badreaches<-which(streamvast$reachdata$habitat==0)
    nbad<-sum(streamvast$reachdata$habitat==0)

    if(tolower(Time)=="year"){
      pad.dates<-NULL
      for(i in 1:nrow(time.table)){
        pad.dates<-as.Date(c(pad.dates,seq(from=time.table$Start[i],to=time.table$End[i],by=15)))
      }
    }else{
      pad.dates<-time.table$Date
    }
    timematch<-sapply(pad.dates,FUN=function(x){return(which(time.table$Start<=x & time.table$End>=x))})

    pad.data<-data.frame(Time=rep(timematch,each=nbad),
                         Runyear=rep(time.table$Runyear[timematch],each=nbad),
                         Date=rep(pad.dates,each=nbad),
                         Year=rep(as.integer(format(pad.dates,"%Y")),each=nbad),
                         Month=rep(as.integer(format(pad.dates,"%m")),each=nbad),
                         Day=rep(as.integer(format(pad.dates,"%d")),each=nbad),
                         Statday=rep(as.integer(pad.dates-time.table$Start[timematch])+1,each=nbad),
                         Reach=rep(badreaches,nrow(time.table)),
                         Species=streamvast$countname,
                         Dist="obs",
                         Counts=0,
                         Density=0,
                         Length=rep(streamvast$reachdata$Length[badreaches],nrow(time.table)),
                         Effort=rep(streamvast$reachdata$Length[badreaches],nrow(time.table)),
                         Lat=rep(streamvast$reachdata$CenterLat[badreaches],nrow(time.table)),
                         Lon=rep(streamvast$reachdata$CenterLon[badreaches],nrow(time.table)),
                         Original=F)

    out.data<-rbind(out.data,pad.data)
  }

  names(out.data)[names(out.data)=="Reach"]<-streamvast$reachname
  names(out.data)[names(out.data)=="Counts"]<-streamvast$countname

  streamvast$vastdata<-out.data
  streamvast$timetable<-time.table
  streamvast$timescale<-Time

  return(streamvast)
}


# Will need a new covariates function, but not immediate priority




#' Get density predictions from a VAST model
#'
#' This function quickly constructs a data frame with VAST density predictions matched to
#' time and space. It also samples and provides upper and lower bounds (90%)
#'
#' @param streamvast A VAST model object
#' @param newdata a data frame with all necessary columns over which to make predictions
#' @param makeauc logical; should auc measures be calculated. Can turn off to save time.
#' @param nsims interger; a number of simulated draws to make from the joint posterior
#' @param bias.correct logical; should the poisson transformation bias be corrected; very slow
#'
#' @return A streamvast object with a dataframe of predictions added
#' @export
#'
#' @examples
VASTpreds<-function(streamvast,newdata,makeauc=T,nsims=100,bias.correct=F){

  if(length(streamvast$preds)>0){warning("Prediction data already found. Overwriting pre-existing data.")}
  if(length(streamvast$eval)>0){warning("Evaluation data already found. Overwriting pre-existing data.")}
  if(length(streamvast$aucdata)>0){warning("AUC or Index data already found. Overwriting pre-existing data.")}

  # first structure the data that will be used for the preddata data frame
  if(missing(newdata)){
    timevec<-streamvast$timetable$Time
    reachvec<-as.data.frame(streamvast$reachdata)[,streamvast$reachname]
    nreach<-length(reachvec)

    # If an auc estimate is desired, need to make sure we have enough data points within each year
    # if the time frame is month, week, etc, that should happen automatrically, so can use the simpler method
    date.dat<-data.frame(Time=vector(),Runyear=vector(),Date=as.Date(NULL))
    if(makeauc & streamvast$timescale=="Year"){
      for(i in 1:nrow(streamvast$timetable)){
        date.dat<-rbind(date.dat,data.frame(Time=streamvast$timetable$Time[i],
                                            Runyear=streamvast$timetable$Runyear[i],
                                            Date=seq(from=streamvast$timetable$Start[i],
                                                     to=streamvast$timetable$End[i],by=15)))
      }
    }else{
      date.dat<-streamvast$timetable[c("Time","Runyear","Date")]
    }

    # Finish filling out the newdata object
    newdata<-data.frame(Time=rep(date.dat$Time,each=nreach),
                        Runyear=rep(date.dat$Runyear,each=nreach),
                        Date=rep(date.dat$Date,each=nreach))
    newdata$Year<-as.integer(format(newdata$Date,"%Y"))
    newdata$Month<-as.integer(format(newdata$Date,"%m"))
    newdata$Day<-as.integer(format(newdata$Date,"%d"))
    newdata$Statday<-as.integer(1+newdata$Date-streamvast$timetable$Start[match(newdata$Time,streamvast$timetable$Time)])
    newdata$Reach<-rep_len(reachvec,length.out=nrow(newdata))

    speciesname<-streamvast$vastmodel$call[["variable_column"]]
    speciesval<-unique(streamvast$vastdata[,speciesname])
    newdata$variable_column<-speciesval
    names(newdata)[ncol(newdata)]<-speciesname

    distname<-streamvast$vastmodel$call[["distribution_column"]]
    distval<-unique(streamvast$vastdata[,distname])
    newdata$distribution_column<-distval
    names(newdata)[ncol(newdata)]<-distname

    reachmatch<-match(newdata$Reach,as.data.frame(streamvast$reachdata)[,streamvast$reachname])
    newdata$Habitat<-streamvast$reachdata$habitat[reachmatch]
    newdata$Length<-streamvast$reachdata$Length[reachmatch]
    newdata$Effort<-newdata$Length
    newdata$Lon<-streamvast$reachdata$CenterLon[reachmatch]
    newdata$Lat<-streamvast$reachdata$CenterLat[reachmatch]
    names(newdata)[names(newdata)=="Reach"]<-streamvast$reachname
  }

  # make predictions; for point predictions, we will use the estimated standard error
  # the preddata object is used for overall predictions, i.e. assuming full survey coverage
  preds<-predict(streamvast$vastmodel,newdata = newdata,what = "p_g",se.fit = T)
  predsims<-Jeremy_sample_variable(obj = streamvast$vastmodel,newdata = newdata,
                                   n_samples = nsims,bias.correct=bias.correct)

  # Especially with the bias correct,sometimes there are a few divide by 0 errors
  good.col<-which(apply(predsims,MARGIN=2,FUN=function(x){return(all(is.na(x)==F))}))
  predsims<-predsims[,good.col]

  predsims2<-predsims/newdata$Length
  evalsims<-tinyVAST::sample_variable(x = streamvast$vastmodel,variable_name = "mu_i",n_samples = nsims)
  evalsims2<-evalsims/streamvast$vastdata$Length

  preddata<-cbind(newdata,data.frame(pred_Count=exp(preds$fit),
                                     pred_Count_SE=exp(preds$se.fit),
                                     pred_Count_lower=exp(preds$fit-1.96*preds$se.fit),
                                     pred_Count_upper=exp(preds$fit+1.96*preds$se.fit),
                                     pred_Density=exp(preds$fit)/newdata$Effort,
                                     pred_Density_SE=exp(preds$se.fit)/newdata$Effort,
                                     pred_Density_lower=exp(preds$fit-1.96*preds$se.fit)/newdata$Effort,
                                     pred_Density_upper=exp(preds$fit+1.96*preds$se.fit)/newdata$Effort))

  #the evaldata object is used for comparisons with observed values, i.e. partial survey coverage
  # consider changing the "pred" prefix to "fit" to avoid confusion
  evals<-predict(streamvast$vastmodel,what="p_g",se.fit=T)
  evaldata<-cbind(streamvast$vastdata,data.frame(fit_Count=exp(evals$fit),
                                                 fit_Count_SE=exp(evals$se.fit),
                                                 fit_Count_lower=exp(evals$fit-1.96*evals$se.fit),
                                                 fit_Count_upper=exp(evals$fit+1.96*evals$se.fit),
                                                 fit_Density=exp(evals$fit)/streamvast$vastdata$Effort,
                                                 fit_Density_SE=exp(evals$se.fit)/streamvast$vastdata$Effort,
                                                 fit_Density_lower=exp(evals$fit-1.96*evals$se.fit)/streamvast$vastdata$Effort,
                                                 fit_Density_upper=exp(evals$fit+1.96*evals$se.fit)/streamvast$vastdata$Effort))

  # now is as good a place as any to grab some summary stats
  streamvast$stats$AIC<-AIC(streamvast$vastmodel)
  streamvast$stats$rhoP<-stats::cor(evaldata$fit_Count,evaldata[,streamvast$countname],method = "pearson")
  streamvast$stats$rhoS<-stats::cor(evaldata$fit_Count,evaldata[,streamvast$countname],method = "spearman")
  streamvast$stats$RMSE<-sqrt(sum((evaldata[,streamvast$countname]-evaldata$fit_Count)^2)/nrow(evaldata))
  streamvast$stats$dev_exp<-streamvast$vastmodel$deviance_explained

  # Make and plot the DHARMa Residuals
  y_ir = replicate(n = 100, expr = streamvast$vastmodel$obj$simulate()$y_i )
  dharmares = DHARMa::createDHARMa( simulatedResponse = y_ir,
                                    observedResponse = streamvast$vastdata[,streamvast$countname],
                                    fittedPredictedResponse = fitted(streamvast$vastmodel) )

  # Add residuals to the eval data frame
  evaldata$Residual_count<-evaldata[,streamvast$countname]-evaldata$fit_Count
  evaldata$Residual_density<-evaldata$Density-evaldata$fit_Density
  evaldata$Residual_pearson<-(evaldata$Density-evaldata$fit_Density)/sqrt(evaldata$fit_Density)
  evaldata$Residual_dharma<-dharmares$scaledResiduals


  # now we will need to compile the aggregate statistics for space
  reachdata.names<-tolower(names(streamvast$reachdata))
  streamname<-names(streamvast$reachdata)[reachdata.names%in%c("strm_na","strm_name","stream","streamname","strmname")]

  spacedata<-data.frame(Reach=as.data.frame(streamvast$reachdata)[,streamvast$reachname],
                        STRM_NAME=as.data.frame(streamvast$reachdata)[,streamname],Habitat=streamvast$reachdata$habitat,
                        Length=streamvast$reachdata$Length,Streamdist=streamvast$reachdata$streamdist,
                        Density=NA,fit_Density=NA,fit_Density_lower=NA,fit_Density_upper=NA,
                        pred_Density=NA,pred_Density_lower=NA,pred_Density_upper=NA)
  names(spacedata)[1]<-streamvast$reachname

  for(i in 1:nrow(spacedata)){
    predspots<-which(preddata[,streamvast$reachname]==spacedata[i,streamvast$reachname])
    evalspots<-which(evaldata[,streamvast$reachname]==spacedata[i,streamvast$reachname])

    if(length(evalspots)>0){
      spacedata$Density[i]<-mean(evaldata[evalspots,"Density"])
      spacedata$fit_Density[i]<-mean(evaldata[evalspots,"fit_Density"])
      sim.eval.means<-apply(matrix(evalsims2[evalspots,],ncol=ncol(evalsims2)),MARGIN = 2,FUN = mean)
      spacedata$fit_Density_lower[i]<-stats::quantile(sim.eval.means,probs=.025,na.rm=T)
      spacedata$fit_Density_upper[i]<-stats::quantile(sim.eval.means,probs=.975,na.rm=T)
    }

    if(length(predspots)>0){
      spacedata$pred_Density[i]<-mean(preddata[predspots,"pred_Density"])
      sim.pred.means<-apply(matrix(predsims2[predspots,],ncol=ncol(predsims2)),MARGIN = 2,FUN = mean)
      spacedata$pred_Density_lower[i]<-stats::quantile(sim.pred.means,probs=.025,na.rm=T)
      spacedata$pred_Density_upper[i]<-stats::quantile(sim.pred.means,probs=.975,na.rm=T)
    }
  }

  # Next we do the same thing for time
  yearlen<-as.integer(floor(mean(streamvast$timetable$End-streamvast$timetable$Start+1)))

  timedata<-data.frame(Time=preddata$Time,Runyear=preddata$Runyear,
                       Date=preddata$Date,Year=preddata$Year,Month=preddata$Month,
                       Day=preddata$Day,Statday=preddata$Statday,
                       Dayseq=preddata$Statday+((preddata$Runyear-min(preddata$Runyear))*yearlen),
                       Density=NA,fit_Density=NA,fit_Density_lower=NA,fit_Density_upper=NA,
                       pred_Density=NA,pred_Density_lower=NA,pred_Density_upper=NA)
  timedata<-unique(timedata)
  timematch<-sapply(evaldata$Date,FUN=function(x){return(which.min(abs(as.integer(x-timedata$Date))))})

  for(i in 1:nrow(timedata)){
    predspots<-which(preddata$Date==timedata$Date[i])
    evalspots<-which(timematch==i)

    if(length(evalspots)>0){
      timedata$Density[i]<-mean(evaldata[evalspots,"Density"])
      timedata$fit_Density[i]<-mean(evaldata[evalspots,"fit_Density"])
      sim.eval.means<-apply(matrix(evalsims2[evalspots,],ncol=ncol(evalsims2)),MARGIN = 2,FUN = mean)
      timedata$fit_Density_lower[i]<-stats::quantile(sim.eval.means,probs=.025,na.rm=T)
      timedata$fit_Density_upper[i]<-stats::quantile(sim.eval.means,probs=.975,na.rm=T)
    }

    if(length(predspots)>0){
      timedata$pred_Density[i]<-mean(preddata[predspots,"pred_Density"])
      sim.pred.means<-apply(matrix(predsims2[predspots,],ncol=ncol(predsims2)),MARGIN = 2,FUN = mean)
      timedata$pred_Density_lower[i]<-stats::quantile(sim.pred.means,probs=.025,na.rm=T)
      timedata$pred_Density_upper[i]<-stats::quantile(sim.pred.means,probs=.975,na.rm=T)
    }
  }

  if(makeauc){
    # Now that time is set, we can also set up the AUC data objects
    aucdata<-data.frame(Runyear=rep(sort(unique(streamvast$timetable$Runyear)),each=nrow(streamvast$reachdata)),
                        Reach=rep(sort(unique(as.data.frame(streamvast$reachdata)[,streamvast$reachname])),
                                  times=length(unique(streamvast$timetable$Runyear))),
                        AUC=NA,pred_AUC=NA,pred_AUC_lower=NA,pred_AUC_upper=NA,pred_AUC_SD=NA)

    auctotals<-data.frame(Runyear=sort(unique(streamvast$timetable$Runyear)),
                          AUC=NA,pred_AUC=NA,pred_AUC_lower=NA,pred_AUC_upper=NA,
                          pred_AUC_SD=NA)

    aucsims<-cbind(data.frame(Runyear=aucdata$Runyear,Reach=aucdata$Reach),
                   matrix(data=NA,nrow=nrow(aucdata),ncol=ncol(predsims)))
    auctotalsims<-cbind(data.frame(Runyear=auctotals$Runyear),
                        matrix(data=NA,nrow=nrow(auctotals),ncol=ncol(predsims)))

    for(y in 1:nrow(auctotals)){

      predspots<-which(preddata$Runyear==auctotals$Runyear[y])
      yearaucdata<-aucdata[aucdata$Runyear==auctotals$Runyear[y],]

      yearpreddata<-preddata[predspots,]
      yearpredsims<-matrix(predsims[predspots,],ncol=ncol(predsims))
      yeartotalpredsim<-matrix(nrow=nrow(yearaucdata),ncol=ncol(predsims))

      for(r in 1:nrow(yearaucdata)){

        predspots2<-which(yearpreddata[,streamvast$reachname]==yearaucdata$Reach[r])
        yearreachpreds<-yearpreddata[predspots2,]

        index<-which(aucdata$Runyear==auctotals$Runyear[y] & aucdata$Reach==yearaucdata$Reach[r])

        if(length(predspots2)>0){
          aucdata$pred_AUC[index]<-MESS::auc(x=yearreachpreds$Statday,y=yearreachpreds$pred_Count)

          yearreachpredaucs<-apply(yearpredsims[predspots2,],MARGIN=2,
                                   FUN=function(x){return(MESS::auc(x=yearreachpreds$Statday,y=x))})
          yeartotalpredsim[r,]<-yearreachpredaucs

          aucdata$pred_AUC_lower[index]<-stats::quantile(yearreachpredaucs,probs = .025)
          aucdata$pred_AUC_upper[index]<-stats::quantile(yearreachpredaucs,probs = .975)
          aucdata$pred_AUC_SD[index]<-sqrt(stats::var(yearreachpredaucs))

          aucsims[index,3:(ncol(predsims)+2)]<-yearreachpredaucs
        }
      }
      auctotals$AUC[y]<-sum(aucdata$AUC[aucdata$Runyear==auctotals$Runyear[y]],na.rm=T)
      auctotals$pred_AUC[y]<-sum(aucdata$pred_AUC[aucdata$Runyear==auctotals$Runyear[y]],na.rm=T)

      auctotalsims[y,2:(ncol(yeartotalpredsim)+1)]<-apply(yeartotalpredsim,MARGIN = 2,FUN = sum)

      auctotals$pred_AUC_lower[y]<-stats::quantile(auctotalsims[y,2:ncol(auctotalsims)],probs=.025,na.rm=T)
      auctotals$pred_AUC_upper[y]<-stats::quantile(auctotalsims[y,2:ncol(auctotalsims)],probs=.975,na.rm=T)
      auctotals$pred_AUC_SD[y]<-sqrt(stats::var(unlist(auctotalsims[y,2:ncol(auctotalsims)])))
    }
    names(aucdata)[names(aucdata)=="Reach"]<-streamvast$reachname
    names(auctotals)[names(auctotals)=="Reach"]<-streamvast$reachname
    names(aucsims)[names(aucsims)=="Reach"]<-streamvast$reachname
    names(auctotalsims)[names(auctotalsims)=="Reach"]<-streamvast$reachname

  }else{
    aucdata<-NULL
    auctotals<-NULL
    aucsims<-NULL
    auctotalsims<-NULL
  }

  # Update all the names so that everything matches
  names(preddata)[names(preddata)=="pred_Count"]<-paste0("pred_",streamvast$countname)
  names(preddata)[names(preddata)=="pred_Count_SE"]<-paste0("pred_",streamvast$countname,"_SE")
  names(preddata)[names(preddata)=="pred_Count_upper"]<-paste0("pred_",streamvast$countname,"_upper")
  names(preddata)[names(preddata)=="pred_Count_lower"]<-paste0("pred_",streamvast$countname,"_lower")

  names(evaldata)[names(evaldata)=="fit_Count"]<-paste0("fit_",streamvast$countname)
  names(evaldata)[names(evaldata)=="fit_Count_SE"]<-paste0("fit_",streamvast$countname,"_SE")
  names(evaldata)[names(evaldata)=="fit_Count_upper"]<-paste0("fit_",streamvast$countname,"_upper")
  names(evaldata)[names(evaldata)=="fit_Count_lower"]<-paste0("fit_",streamvast$countname,"_lower")
  names(evaldata)[names(evaldata)=="Residual_count"]<-paste0("Residual_",streamvast$countname)

  names(spacedata)[names(spacedata)=="Reach"]<-streamvast$reachname


  streamvast$preds<-preddata
  streamvast$eval<-evaldata
  streamvast$spacedata<-spacedata
  streamvast$timedata<-timedata
  streamvast$dharmafit<-dharmares
  streamvast$aucdata<-list(aucdata=aucdata,auctotals=auctotals)
  streamvast$sims<-list(predsims=predsims,evalsims=evalsims,aucsims=aucsims,
                        auctotalsims=auctotalsims,escapesims=NULL)
  return(streamvast)
}


#' Makes a line plot that shows tributaries by stream distance
#'
#' @param streamvast  a streamvast object with predictions or evaluations
#' @param plotvar     the column name in data to be plotted
#' @param streamname  the column name in the reaches with the stream name
#' @param usepreds    logical, should the prediction or evaluation data be used, if absent, the function will search
#' @param ylab        character, a label for the y axis
#' @param title       a title to display above the graph
#' @param show.names  character; a vector of names to highlight, useful for large watersheds
#'
#' @return A plot of the variable with respect ot stream distance
#' @export
#'
#' @examples
plotStream<-function(streamvast,plotvar="pred_Density",ylab,streamname,usepreds=T,title,show.names="all"){

  reachdata<-as.data.frame(streamvast$reachdata)
  if(streamname%in%names(reachdata)==F){stop("Streamname not found. Check names and try again")}
  if(missing(plotvar)){plotvar<-"pred_Density"}

  if(missing(ylab)){ylab<-plotvar}

  if(usepreds){
    vardata<-streamvast$preds
  }else{
    vardata<-streamvast$eval
  }
  if(plotvar%in%names(vardata)==F){stop("Plotvar not found in indicated dataset. Check settings and try again.")}

  if(show.names[1]=="all"){
    show.names<-unique(reachdata[,streamname])
    show.names[-which(show.names=="NA")]
  }

  # Note, can't use the already established parent in this case
  segdata<-data.frame(Reach=reachdata[,streamvast$reachname],name=reachdata[,streamname],var=0,
                      dist=reachdata$streamdist,parent=reachdata$parent,parentvar=0,
                      parentdist=reachdata$streamdist[match(reachdata$parent,reachdata[,streamvast$reachname])])
  segdata$name[segdata$name%in%show.names==F]<-NA

  vardata2<-stats::aggregate(vardata[,plotvar],by=list(Reach=vardata[,streamvast$reachname]),FUN=mean)
  segdata$var[match(vardata2[,streamvast$reachname],segdata$Reach)]<-vardata2$x

  segdata$parentvar<-segdata$var[segdata$parent]
  segdata$parentdist<-segdata$dist[segdata$parent]
  segdata$parentvar[is.na(segdata$parentvar)]<-0
  segdata$parentdist[is.na(segdata$parentdist)]<-0

  outplot<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=segdata,ggplot2::aes(x=dist,xend=parentdist,y=var,yend=parentvar),col="gray50")+
    ggplot2::geom_point(data=segdata,ggplot2::aes(x=dist,y=var),shape=16,col="gray50")+
    ggplot2::geom_segment(data=subset(segdata,is.na(name)==F),
                          ggplot2::aes(x=dist,xend=parentdist,y=var,yend=parentvar,col=name),
                          linewidth=1)+
    ggplot2::geom_point(data=subset(segdata,is.na(name)==F),
                        ggplot2::aes(x=dist,y=var,col=name),shape=16)+
    ggplot2::theme_bw()+ggplot2::theme(legend.title = ggplot2::element_blank())+
    viridis::scale_color_viridis(discrete = T)+
    ggplot2::xlab("Stream Distance (km)")+ggplot2::ylab(ylab)

  if(missing(title)==F){outplot<-outplot+ggplot2::ggtitle(title)}

  return(outplot)
}



plotStream2<-function(streamvast,streamname,usepreds=T,title,show.names="all"){

  reachdata<-as.data.frame(streamvast$reachdata)
  if(streamname%in%names(reachdata)==F){stop("Streamname not found. Check names and try again")}
  if(show.names[1]=="all"){
    show.names<-unique(reachdata[,streamname])
    show.names[-which(show.names=="NA")]
  }

  # Note, can't use the already established parent in this case
  segdata<-data.frame(Reach=reachdata[,streamvast$reachname],name=reachdata[,streamname],
                      dist=reachdata$streamdist,parent=reachdata$parent,parentvar=NA,parentvar_lower=NA,parentvar_upper=NA,
                      parentdist=reachdata$streamdist[match(reachdata$parent,reachdata[,streamvast$reachname])])
  segdata$name[segdata$name%in%show.names==F]<-NA

  segdata<-cbind(segdata,streamvast$spacedata[match(streamvast$spacedata[,streamvast$reachname],segdata$Reach),
                                              c("pred_Density","pred_Density_lower","pred_Density_upper")])


  segdata$parentvar<-segdata$pred_Density[segdata$parent]
  segdata$parentvar_lower<-segdata$pred_Density_lower[segdata$parent]
  segdata$parentvar_upper<-segdata$pred_Density_upper[segdata$parent]

  segdata$parentvar[is.na(segdata$parentvar)]<-0
  segdata$parentvar_lower[is.na(segdata$parentvar_lower)]<-0
  segdata$parentvar_upper[is.na(segdata$parentvar_upper)]<-0
  segdata$parentdist[is.na(segdata$parentdist)]<-0

  outplot<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=segdata,ggplot2::aes(x=dist,xend=parentdist,y=pred_Density,yend=parentvar),col="gray50")+
    ggplot2::geom_point(data=segdata,ggplot2::aes(x=dist,y=pred_Density),shape=16,col="gray50")+
    ggplot2::geom_segment(data=subset(segdata,is.na(name)==F),
                          ggplot2::aes(x=dist,xend=parentdist,y=pred_Density,yend=parentvar,col=name),
                          linewidth=1)+
    ggplot2::geom_segment(data=subset(segdata,is.na(name)==F),
                          ggplot2::aes(x=dist,xend=parentdist,y=pred_Density_lower,yend=parentvar_lower,col=name),
                          linewidth=1,linetype=3)+
    ggplot2::geom_segment(data=subset(segdata,is.na(name)==F),
                          ggplot2::aes(x=dist,xend=parentdist,y=pred_Density_upper,yend=parentvar_upper,col=name),
                          linewidth=1,linetype=3)+
    ggplot2::geom_point(data=subset(segdata,is.na(name)==F),
                        ggplot2::aes(x=dist,y=pred_Density,col=name),shape=16)+
    ggplot2::theme_bw()+ggplot2::theme(legend.title = ggplot2::element_blank())+
    viridis::scale_color_viridis(discrete = T)+
    ggplot2::xlab("Stream Distance (km)")+ggplot2::ylab("Predicted Density")

  if(missing(title)==F){outplot<-outplot+ggplot2::ggtitle(title)}

  return(outplot)
}

#' Maps, the desired data, with options for facetting
#'
#' @param streamvast   a streamvast object with a prediction data frame
#' @param mapvar       the column in data to map
#' @param facet        a column to determine the facets, or a vector of values
#' @param FUN          a function like mean or median to summarize data
#' @param legendname   character; a name for the legend
#' @param background   a sf object to put in the background
#' @param subset       a logical vector indicating which pred values to use
#' @param make.labels  should the reaches be labelled, turn off for multiple years
#' @param xaxis.breaks Vector of x values for axis lines; useful to reduce crowding
#' @param yaxis.breaks Vector of y values for axis lines; useful to reduce crowding
#' @param palette      a viridis color palette for plotting
#' @param max          a maximum value to cap the scale; values above the max are set equal to the max
#'
#' @return A colored heatmap of the desired variable mapped over the stream network
#' @export
#'
#' @examples
plotPredictionMap<-function(streamvast,mapvar="pred_Density",facet=NA,FUN="mean",background,subset,
                            make.labels=F,xaxis.breaks=NA, yaxis.breaks=NA,palette="turbo",legendname,
                            max=Inf){

  if(missing(subset)){subset<-rep(T,nrow(streamvast$preds))}

  # Identify which column from preds to plot
  pred.col<-which(tolower(names(streamvast$preds))==tolower(mapvar))
  if(length(pred.col)==0){
    pred.cols<-which(unlist(lapply(strsplit(names(streamvast$preds),split="_"),FUN=function(x){return(x[1]=="pred")})))
    stop("Please choose mapvar from: ",paste(names(streamvast$preds)[pred.cols],collapse = ", "))
  }

  if(is.na(facet)==F){
    facet.name<-"group"

    if(length(facet)==1){
      if(is.character(facet)){
        facet.name<-facet
      }
      datareach<-stats::aggregate(streamvast$preds[subset,pred.col],by=list(Reach=streamvast$preds[subset,streamvast$reachname],
                                                                            streamvast$preds[subset,facet]),FUN=FUN)
    }
    names(datareach)[2]<-facet.name
    facet.vec<-rep(sort(unique(datareach[,facet])),each=nrow(streamvast$reachdata))
  }else{
    datareach<-stats::aggregate(streamvast$preds[subset,pred.col],by=list(Reach=streamvast$preds[subset,streamvast$reachname]),FUN=FUN)
    facet.vec<-NA
    facet.name<-"group"
  }

  # add labels, gets rather complicated
  if(make.labels){
    midpoints<-as.data.frame(streamvast$reachdata)[,c("CenterLon","CenterLat")]
    midpoints$Reach<-as.data.frame(streamvast$reachdata)[,streamvast$reachname]
    midpoints<-st_as_sf(midpoints,coords=1:2,crs="wgs84")

    labeldat<-data.frame(X=rep(NA,nrow(streamvast$reachdata)),Y=NA,labels=as.data.frame(streamvast$reachdata)[,streamvast$reachname])
    meanlength<-mean(sf::st_length(streamvast$reachdata))
    distmat<-sf::st_distance(streamvast$reachdata)
    for(i in 1:nrow(streamvast$reachdata)){

      #grab all nearby features, combine to a multi-line, compute the buffer, and pick a point on the buffer
      nearby<-which(distmat[i,]<meanlength*1.5)
      nearby.line<-sf::st_as_sf(sf::st_combine(streamvast$reachdata[nearby,]))

      buf<-sf::st_cast(sf::st_buffer(nearby.line,dist=meanlength/2),"LINESTRING")
      nearest.point<-sf::st_as_sf(sf::st_nearest_points(y = midpoints.sf[i,],x=buf))

      labeldat$X[i]<-sf::st_coordinates(nearest.point)[1,1]
      labeldat$Y[i]<-sf::st_coordinates(nearest.point)[1,2]
    }
  }

  # Set up the ggplot data
  reachplot<-sf::st_as_sf(rep(sf::st_geometry(streamvast$reachdata),times=length(facet.vec)))
  reachplot$reach<-rep(as.data.frame(streamvast$reachdata)[,streamvast$reachname],times=length(facet.vec))
  reachplot$vastid<-streamvast$reachdata$vastid
  reachplot$habitat<-streamvast$reachdata$habitat
  reachplot$group<-facet.vec
  reachplot$preds<-NA
  for(i in 1:nrow(reachplot)){
    if(reachplot$habitat[i]==T){
      if(length(facet.vec)>1){
        reachplot$preds[i]<-datareach$x[datareach$Reach==reachplot$reach[i] &
                                          datareach[,facet.name]==reachplot$group[i]]
      }else{
        reachplot$preds[i]<-datareach$x[datareach$Reach==reachplot$reach[i]]
      }
    }
  }
  bounds<-sf::st_bbox(sf::st_transform(reachplot,crs = "wgs84"))

  reachplot$preds[reachplot$preds>max]<-max

  if(mapvar=="pred_Density"){plotname0<-"Density"}
  if(mapvar==paste0("pred_",streamvast$countname)){plotname0<-streamvast$countname}

  if(FUN=="mean"){plot.name<-paste0("Mean ",mapvar)}
  if(FUN=="var"){plot.name<-paste0(mapvar, " Variance")}

  if(missing(legendname)){legendname<-paste(FUN,mapvar)}

  if(missing(background)){
    map<-ggplot2::ggplot()
  }else{
    map<-ggplot2::ggplot()+
      ggplot2::geom_sf(data=sf::st_transform(background,crs=sf::st_crs(streamvast$reachdata)),col="gray50",alpha=.5)
  }

  map<-map+
    ggplot2::geom_sf(data=reachplot,lwd=2)+
    ggplot2::geom_sf(data=reachplot,ggplot2::aes(col=preds),lwd=1.25)+ggplot2::theme_bw()+
    viridis::scale_color_viridis(option = palette,name=legendname,na.value="grey75")
  if(length(facet.vec)>1){
    map<-map+ggplot2::facet_wrap(~group)
  }
  if(is.na(xaxis.breaks[1])==F){
    map<-map+ggplot2::scale_x_continuous(breaks = xaxis.breaks,limits = sf::st_bbox(reachplot)[c(1,3)])
  }else{
    map<-map+ggplot2::scale_x_continuous(limits = sf::st_bbox(reachplot)[c(1,3)])+
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  if(is.na(yaxis.breaks[1])==F){
    map<-map+ggplot2::scale_y_continuous(breaks = yaxis.breaks,limits = sf::st_bbox(reachplot)[c(2,4)])
  }else{
    map<-map+ggplot2::scale_y_continuous(limits = sf::st_bbox(reachplot)[c(2,4)])+
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }
  if(make.labels){
    map<-map+ggplot2::geom_text(data=labeldat,ggplot2::aes(x=X,y=Y,label=labels),check_overlap = T)
  }
  return(map)
}


#' Plot the dHARMa scaled residuals, accounts for dummy data
#'
#' @param streamvast A streamvast object with fitted model and predictions
#' @param span A smoothing parameter for the line
#'
#' @return a plot mimicing the style of dHARMa, but accounts for dummy data
#' @export
#'
#' @examples
Dharmaplot<-function(streamvast,span=.1){

  scale.ranks<-rank(streamvast$dharmafit$fittedPredictedResponse)/
    length(streamvast$dharmafit$fittedPredictedResponse)

  bins<-0:100/100
  median.vec<-rep(NA,100)
  lower.vec<-rep(NA,100)
  upper.vec<-rep(NA,100)
  for(i in 1:(length(bins)-1)){
    start=bins[i]
    stop=bins[i+1]
    rank.picks<-which(scale.ranks>start & scale.ranks<=stop)
    median.vec[i]<-stats::median(streamvast$dharmafit$scaledResiduals[rank.picks])
    lower.vec[i]<-stats::quantile(streamvast$dharmafit$scaledResiduals[rank.picks],probs=.25)
    upper.vec[i]<-stats::quantile(streamvast$dharmafit$scaledResiduals[rank.picks],probs=.75)
  }

  pt.dat<-data.frame(x=scale.ranks,y=streamvast$dharmafit$scaledResiduals)
  line.dat<-data.frame(x=1:100/100, y=stats::predict(stats::loess(y~x,data=data.frame(x=1:100/100,y=median.vec),
                                                                  span = span),newdata = data.frame(x=1:100/100)))
  line.dat2<-data.frame(x=1:100/100, y=stats::predict(stats::loess(y~x,data=data.frame(x=1:100/100,y=lower.vec),
                                                                   span = span),newdata = data.frame(x=1:100/100)))
  line.dat3<-data.frame(x=1:100/100, y=stats::predict(stats::loess(y~x,data=data.frame(x=1:100/100,y=upper.vec),
                                                                   span = span),newdata = data.frame(x=1:100/100)))

  outplot<-ggplot2::ggplot()+ggplot2::geom_point(data=pt.dat,ggplot2::aes(x=x,y=y),alpha=.25)+
    ggplot2::geom_line(data=line.dat,ggplot2::aes(x=x,y=y),col=2,linewidth=1.5)+
    ggplot2::geom_line(data=line.dat2,ggplot2::aes(x=x,y=y),col=2,linewidth=1)+
    ggplot2::geom_line(data=line.dat3,ggplot2::aes(x=x,y=y),col=2,linewidth=1)+
    ggplot2::geom_hline(yintercept=c(.05,.25,.4,.5,.6,.75,.95),col=2,linetype=2)+
    ggplot2::theme_bw()+ggplot2::xlab("Rank-transformed Predictions")+ggplot2::ylab("Scaled Residuals")

  return(outplot)
}


# For internal use by the preds functions. This is my hacked together variant that I use to sample
# the joint posterior distribution and simulataneously generate predictions at user specified points
# as opposed to the points in the original data
#' Title
#'
#' @param obj a tinyVAST model object
#' @param newdata a dataset for making new predictions
#' @param what what predictions are desired, almost always "mu_g"
#' @param n_samples integer of samples
#' @param sample_fixed logical; almost always T
#' @param seed integer; a seed for the random draws
#' @param bias.correct logical; should the transformation bias correction be applied; very slow
#'
#' @return an object with simulated predictions on newdata
#' @export
#'
#' @examples
Jeremy_sample_variable<-function(obj,newdata,what="mu_g",n_samples=100,sample_fixed=TRUE,seed=123456,bias.correct=F){

  require(Matrix)
  if (!("jointPrecision" %in% names(obj$sdrep))) {
    stop("jointPrecision not present in x$sdrep; please re-run with `getJointPrecision=TRUE`")
  }
  ParHat = obj$obj$env$parList()
  Intersect = intersect(names(obj$rep), names(ParHat))
  if (isFALSE(all.equal(obj$rep[Intersect], ParHat[Intersect]))) {
    stop("Duplicate entries in `Obj$report()` and `Obj$env$parList()` are not identical when calling `sample_variable`")
  }
  Output = c(obj$rep, ParHat)
  # if (isFALSE(variable_name %in% names(Output))) {
  #   stop(variable_name, " not found in `Obj$report()` or `Obj$env$parList()`; please choose check your requested variable name from available list: ",
  #        paste(names(Output), collapse = ", "))
  # }
  rmvnorm_prec <- function(mu, prec, n.sims, seed) {
    set.seed(seed)
    z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
    L <- Matrix::Cholesky(prec, super = TRUE)
    z <- solve(L, z, system = "Lt")
    z <- solve(L, z, system = "Pt")
    z <- as.matrix(z)
    return(mu + z)
  }
  if (sample_fixed == TRUE) {
    u_zr = rmvnorm_prec(mu = obj$obj$env$last.par.best, prec = obj$sdrep$jointPrecision,
                        n.sims = n_samples, seed = seed)
  } else {
    u_zr = obj$obj$env$last.par.best %o% rep(1, n_samples)
    MC = obj$obj$env$MC(keep = TRUE, n = n_samples, antithetic = FALSE)
    u_zr[obj$obj$env$random, ] = attr(MC, "samples")
  }
  message("# Obtaining samples from predictive distribution for variable ",
          what)
  for (rI in 1:n_samples) {
    if (rI%%max(1, floor(n_samples/10)) == 0) {
      message("  Running sample ", rI, " of ", n_samples)
    }
    Report = obj$obj$report(par = u_zr[, rI])
    ParHat = obj$obj$env$parList(x = u_zr[, rI][obj$obj$env$lfixed()],
                                 par = u_zr[, rI])
    if (isFALSE(all.equal(obj$rep[Intersect], ParHat[Intersect]))) {
      stop("Duplicate entries in `obj$obj$report()` and `obj$obj$env$parList()` are not identical when calling `sample_variable`")
    }

    # perhaps add an if statement dependent on newdata
    # this part taken from tinyVAST.predict
    tmb_data2 = add_predictions( object = obj,
                                 newdata = newdata,
                                 remove_origdata = FALSE )

    # Rebuild object
    newobj = TMB::MakeADFun( data = tmb_data2,    # insert newdata
                             parameters = ParHat, # insert new parameters
                             map = obj$tmb_inputs$tmb_map,
                             random = obj$tmb_inputs$tmb_random,
                             profile = obj$internal$control$profile,
                             hessian = T,
                             DLL = "tinyVAST" )
    newobj$env$beSilent()
    if(bias.correct){
      newobj.sd<-TMB::sdreport(obj = newobj,
                               getReportCovariance = F,
                               bias.correct = F)
      out0 <- exp(newobj.sd$value-((newobj.sd$sd)^2)/2)
    }else{
      out0<-newobj$report()[[what]]
    }
    if(rI==1){
      out<-as.matrix(out0)
    }else{
      out<-cbind(out,out0) # cbind is too slow, ought to change
    }
  }
  return(out)
}
