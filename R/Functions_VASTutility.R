# This script will hold several utility functions designed to help analyze the output of VAST models

#' ConstructStreamVAST
#' This function constructs objects with class "StreamVAST". Depending on the arguments supplied, it will either
#' initialize the object from countdata and reachdata, or it will add details to an already defined StreamVAST object.
#'
#' @param StreamVAST An StreamVAST object to append details
#' @param countdata A dataframe containing count data, such as outputed from AssembleReddData
#' @param reachdata A sf object with data for reaches, such as from AssignReaches
#' @param surveydata A optional data frame or sf object with survey information, such as from MakeSurveyTracks, doesn't do much right now
#' @param countname A column name from countdata to use
#' @param reachname a column name from reachdata to use
#'
#' @return A streamvast object with appropriate formatting
#' @export
#'
#' @examples
ConstructStreamVAST<-function(countdata,reachdata,surveydata,
                              countname,reachname,unitconv){

  # if the StreamVAST object is missing, we will construct it from countdata and reachdata
  if(missing(countdata) | missing(reachdata)){
    stop("Arguments 'countdata' and 'reachdata' are required to initialize StreamVAST")
  }
  print("Initializing StreamVAST based on countdata and reachdata")
  if(class(countdata)[1]=="character"){countdata<-utils::read.csv(countdata)}
  if(class(reachdata)[1]=="character"){reachdata<-sf::st_read(reachdata)}
  if(class(reachdata)[1]!="sf"){stop("Argument reachdata must be a sf object with LINESTRING geometries")}
  if(missing(surveydata)){surveydata<-NA}

  if(missing(unitconv)){
    x1<-sf::st_length(reachdata)[1]
    x2<-x1
    units(x2)<-units::make_units(km)
    unitconv<-as.numeric(x2)/as.numeric(x1)
    if(abs(unitconv-1)>.1){print("Converting length units to kilometers (Default)")}
  }

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

  #Fix an occasion naming issue
  if("parent.distance"%in%names(countdata)){
    names(countdata)[names(countdata)=="parent.distance"]<-"prnt_ds"
    }

  # Calculate some values that aren't in the data already
  reachdata$Length<-as.numeric(sf::st_length(reachdata))*unitconv  # convert feet to km
  countdata$Area<-reachdata$Length[countdata[,reachname]]
  countdata$Effort<-countdata$Effort*unitconv
  countdata$prnt_ds<-as.numeric(countdata$prnt_ds*unitconv)
  reachdata$prnt_ds<-as.numeric(reachdata$prnt_ds*unitconv)
  countdata$Density<-countdata[,countname]/countdata$Effort

  # account for different capitalizations of day
  names(countdata)[names(countdata)%in%c("Day","day")]<-"Day"

  data.dates<-as.Date(paste(countdata$Year,countdata$Day,sep="_"),format="%Y_%j")
  countdata$Month<-format(data.dates,format="%b")
  countdata$Week<-ceiling(countdata$Day/7)

  # add stream distance to reachdata and add vastids that will account for non-habitat areas
  reachdata$streamdist<-NA
  reachdata$vastid<-NA
  reachdata$vastparent<-NA
  reachdata$streamdist[is.na(reachdata$parent)]<-reachdata$Length[is.na(reachdata$parent)]/2
  reachdata$vastid[is.na(reachdata$parent)]<-1


  needs.doing<-as.data.frame(reachdata)[is.na(reachdata$streamdist),reachname]
  already.done<-as.data.frame(reachdata)[is.na(reachdata$streamdist)==F,reachname]
  while(length(needs.doing)>0){
    up.next<-min(as.data.frame(reachdata)[as.data.frame(reachdata)[,reachname]%in%needs.doing & reachdata$parent%in%already.done,reachname])
    index<-which(as.data.frame(reachdata)[,reachname]==up.next)
    parent.index<-which(as.data.frame(reachdata)[,reachname]==reachdata$parent[index])
    reachdata$streamdist[index]<-reachdata$streamdist[parent.index]+reachdata$prnt_ds[index]

    if(reachdata$habitat[index]==T){
      reachdata$vastid[index]<-max(reachdata$vastid,na.rm=T)+1
      reachdata$vastparent[index]<-reachdata$vastid[parent.index]
    }

    needs.doing<-needs.doing[needs.doing!=up.next]
    already.done<-c(already.done,up.next)
  }

  habitat.reaches<-reachdata[reachdata$habitat==T,]
  countdata$vastid<-reachdata$vastid[match(countdata$Reach,reachdata$Reach)]
  countdata$Vastparent<-reachdata$vastparent[match(countdata$Reach,reachdata$Reach)]

  # Next set up the spatial inputs. These are consistent and do not change unless you load a new data set
  # we designate the prediction and sampling knots as being at the midpoint of each reach
  midpoints<-as.data.frame(sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=habitat.reaches,sample = .5),
                                                               crs = "wgs84")))[,1:2]
  midpoints$ID<-habitat.reaches$vastid

  # important format stuff
  # child is the same as reach id
  # parent is the reach downstream of the child
  # The root has parent = 0 and dist = Inf
  vastinput<-midpoints[order(midpoints$ID),]
  names(vastinput)<-c("Lon","Lat","child_i")
  vastinput$Area_km2<-habitat.reaches$Length[match(vastinput$child_i,habitat.reaches$vastid)]
  vastnetworkLL<-vastinput[,1:2]
  vastnetworkLL$child_s<-vastinput$child_i
  vastnetworkLL$parent_s<-habitat.reaches$vastparent[match(vastnetworkLL$child_s,habitat.reaches$vastid)]
  vastnetworkLL$dist_s<-habitat.reaches$prnt_ds[match(vastinput$child_i,habitat.reaches$vastid)]
  vastnetworkLL$parent_s[is.na(vastnetworkLL$parent_s)]<-0
  vastnetworkLL$dist_s[vastnetworkLL$parent_s==0]<-Inf
  vastnetwork<-vastnetworkLL[,3:5]

  if(missing(surveydata)){surveydata<-NA}
  if(class(surveydata)[1]=="character"){surveydata<-sf::st_read(surveydata)}
  if(any(class(surveydata)%in%c("sf","data.frame"))){names(surveydata)[names(surveydata)%in%c("Day","day")]<-"Day"}
  if(any(countdata$habitat==F)){
    warning(paste("Removed ",sum(countdata$habitat==F)," count records in non-habitat areas"))
    print("Summary of excluded records")
    print(utils::head(subset(countdata,habitat==F)))
  }
  countdata<-subset(countdata,habitat==T)

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
  vastmodel<-NULL
  stats<-list(AIC=NULL,RMSE=NULL,rhoP=NULL,rhoS=NULL)
  preds<-NULL
  eval<-NULL
  spacedata<-NULL
  timedata<-NULL
  aucdata<-NULL
  dharmafit<-NULL
  covariates<-list(covariatedata=NULL,pformula=stats::formula(~0),pconfig=NULL,
                                           dformula=stats::formula(~0),dconfig=NULL)

  StreamVAST<-list(countdata=countdata,reachdata=reachdata,surveydata=surveydata,
                   vastinput=vastinput,vastnetwork=vastnetwork,vastnetworkLL=vastnetworkLL,
                   title=NULL,countname=countname,reachname=reachname,vastdata=NULL,
                   timetable=NULL,timescale=NULL,vastmodel=NULL,unitconv=unitconv,
                   stats=list(AIC=NULL,RMSE=NULL,rhoP=NULL,rhoS=NULL),preds=NULL,
                   eval=NULL,spacedata=NULL,timedata=NULL,aucdata=NULL,dharmafit=NULL,
                   covariates=list(covariatedata=NULL,pformula=stats::formula(~0),pconfig=NULL,
                                   dformula=stats::formula(~0),dconfig=NULL))

  class(StreamVAST)<-"streamvast"
  return(StreamVAST)
}


#' This function takes a streamvast object defines a temporal sampling frame and formats everything for use with vast
#'
#' @param streamvast a streamvast object, minimally with countdata, reachdata defined
#' @param startdate dates to truncate or expand the data in yyyy-mm-dd format
#' @param enddate dates to truncate or expand the data in yyyy-mm-dd format
#' @param padzero should survey with zero observations be added at start/end dates
#' @param Time character giving the timescale to use
#'
#' @return A streamvast object with formated with the appropriate temporal frame
#' @export
#'
#' @examples
SetTemporalFrame<-function(streamvast,startdate=NA,enddate=NA,padzero=T,Time="Year"){

  if(length(streamvast$timetable)>0){warning("Temporal sampling frame already defined for this object. Overwriting pre-existing temporal frame.")}

  if(tolower(Time)=="year"){timeformat<-"%Y"}
  if(tolower(Time)=="month"){timeformat<-"%Y-%m"}
  if(tolower(Time)%in%c("biweek","week","day")){timeformat<-"%Y-%j"}
  tdivide<-ifelse(tolower(Time)=="biweek",14,ifelse(tolower(Time)=="week",7,1))

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

  data.dates2<-c(startdate,data.dates,enddate)

  # prune data, if necessary
  good.data<-streamvast$countdata[data.dates>=startdate & data.dates<=enddate,]
  good.dates<-data.dates[data.dates>=startdate & data.dates<=enddate]
  habitat.reaches<-streamvast$reachdata[streamvast$reachdata$habitat==T,]
  habitat.ids<-habitat.reaches$vastid

  # Now need to go through and get tables that convert dates into sequential time steps
  minyear<-as.integer(format(startdate,"%Y"))
  maxyear<-as.integer(format(enddate,"%Y"))
  nyears<-length(minyear:maxyear)
  if(Time=="Year"){
    time.table<-data.frame(Time=1:nyears,Year=minyear:maxyear,Key=minyear:maxyear,
                           Refdate=as.Date(paste(minyear:maxyear,round(mean(good.data$Day)),sep="-"),format="%Y-%j"))
    good.data$Key<-good.data$Year
  }else{
    firstt<-ceiling(as.integer(format(startdate,strsplit(timeformat,split="-")[[1]][2]))/tdivide)
    lastt<-ceiling(as.integer(format(enddate,strsplit(timeformat,split="-")[[1]][2]))/tdivide)
    mint<-min(ceiling(as.integer(format(data.dates2,strsplit(timeformat,split="-")[[1]][2]))/tdivide))
    maxt<-max(ceiling(as.integer(format(data.dates2,strsplit(timeformat,split="-")[[1]][2]))/tdivide))
    firstyear<-firstt:maxt
    lastyear<-mint:lastt
    ntimesteps<-length(mint:maxt)*(nyears-2) +length(firstyear)+length(lastyear)

    time.table<-data.frame(Time=1:ntimesteps,
                           X1=c(firstyear,rep(mint:maxt,times=nyears-2),lastyear),
                           Year=c(rep(minyear,length(firstyear)),
                                  rep((minyear+1):(maxyear-1),each=length(mint:maxt)),
                                  rep(maxyear,length(lastyear))))
    names(time.table)[2]<-Time

    time.table$Key<-paste(time.table$Year,time.table[,Time],sep="-")

    if(tolower(Time)=="month"){time.table$Refdate<-as.Date(paste(time.table$Year,time.table[,Time],15,sep="-"),format="%Y-%m-%d")}
    if(tolower(Time)=="biweek"){
      days<-ifelse((time.table[,Time]-1)*14+7>365,365,(time.table[,Time]-1)*14+7)
      time.table$Refdate<-as.Date(paste(time.table$Year,days,sep="-"),format="%Y-%j")
    }
    if(tolower(Time)=="week"){
      days<-ifelse((time.table[,Time]-1)*7+4>365,365,(time.table[,Time]-1)*7+4)
      time.table$Refdate<-as.Date(paste(time.table$Year,days,sep="-"),format="%Y-%j")
    }
    if(tolower(Time)=="day"){time.table$Refdate<-as.Date(paste(time.table$Year,time.table[,Time],sep="-"),format="%Y-%j")}

    good.data$Key<-paste(good.data$Year,ceiling(as.integer(format(good.dates,format=strsplit(timeformat,split="-")[[1]][2]))/tdivide),sep="-")
  }
  time.table$Day<-as.integer(format(time.table$Refdate,format="%j"))

  # Set up some reach info
  midpoints<-as.data.frame(sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=habitat.reaches,sample = .5),crs = "wgs84")))[,1:2]
  midpoints$vastid<-habitat.ids

  # Get the preliminary output object that has just the original data
  out.data<-data.frame(Time=time.table$Time[match(good.data$Key,time.table$Key)],
                       Year=good.data$Year,
                       Day=good.data$Day,
                       Month=good.data$Month,
                       Reach=good.data$Reach,
                       vastid=good.data$vastid,
                       Counts=good.data[,streamvast$countname],
                       Density=good.data[,streamvast$countname]/good.data$Effort,
                       Area=good.data$Area,
                       Effort=good.data$Effort,
                       Lat=midpoints[match(x=good.data$vastid,midpoints$vastid),2],
                       Lon=midpoints[match(x=good.data$vastid,midpoints$vastid),1],
                       original=T,
                       dummy=F)

  # Now make padding data, if adding zeros to the ends of each season
  if(padzero){
    if(tolower(Time)=="year"){stop("Padding with zeros doesn't make sense with an annual time scale")}

    front.time<-stats::aggregate(time.table$Time,by=list(time.table$Year),FUN=min)
    back.time<-stats::aggregate(time.table$Time,by=list(time.table$Year),FUN=max)

    # for each year, we check which reaches were sampled in the earliest and latest time frame for that year
    for(y in 1:nyears){
      yeardata<-subset(out.data,Year==((minyear:maxyear)[y]))

      if(nrow(yeardata)>0){
        front.check<-stats::aggregate(yeardata$Time,by=list(yeardata$vastid),FUN=function(x){return(min(x)==front.time$x[y])})
        back.check<-stats::aggregate(yeardata$Time,by=list(yeardata$vastid),FUN=function(x){return(max(x)==back.time$x[y])})

        not.sampled<-habitat.ids[which(habitat.ids%in%yeardata$vastid==F)]
        needs.front<-sort(c(not.sampled,front.check$Group.1[front.check$x==F]))
        needs.back<-sort(c(not.sampled,back.check$Group.1[back.check$x==F]))
      }else{
        needs.front<-habitat.ids
        needs.back<-habitat.ids
      }
      if(length(needs.front)>0){
        front.data<-data.frame(Time=front.time$x[y],
                               Year=(minyear:maxyear)[y],
                               Day=time.table$Day[time.table$Time==front.time$x[y]],
                               Month=format(time.table$Refdate[match(front.time$x[y],time.table$Time)],"%b"),
                               Reach=as.data.frame(habitat.reaches)[match(needs.front,habitat.reaches$vastid),streamvast$reachname],
                               vastid=needs.front,
                               Counts=0,
                               Density=0,
                               Area=as.numeric(streamvast$reachdata$Length[needs.front]),
                               Effort=as.numeric(streamvast$reachdata$Length[needs.front]),
                               Lat=midpoints[match(needs.front,habitat.ids),2],
                               Lon=midpoints[match(needs.front,habitat.ids),1],
                               original=F,
                               dummy=F)
        out.data<-rbind(out.data,front.data)
      }
      if(length(needs.back)>0){
        back.data<-data.frame(Time=back.time$x[y],
                              Year=(minyear:maxyear)[y],
                              Day=time.table$Day[time.table$Time==back.time$x[y]],
                              Month=format(time.table$Refdate[match(back.time$x[y],time.table$Time)],"%b"),
                              Reach=as.data.frame(habitat.reaches)[match(needs.back,habitat.reaches$vastid),streamvast$reachname],
                              vastid=needs.back,
                              Counts=0,
                              Density=0,
                              Area=as.numeric(streamvast$reachdata$Length[needs.back]),
                              Effort=as.numeric(streamvast$reachdata$Length[needs.back]),
                              Lat=midpoints[match(needs.back,habitat.ids),2],
                              Lon=midpoints[match(needs.back,habitat.ids),1],
                              original=F,
                              dummy=F)
        out.data<-rbind(out.data,back.data)
      }
    }
  }

  # Now to make the dummy data that signals VAST to interpolate a prediction
  checklist<-data.frame(Time=rep(time.table$Time,nrow(habitat.reaches)),
                        vastid=rep(habitat.reaches$vastid,each=nrow(time.table)),
                        Reach=rep(as.data.frame(habitat.reaches)[,streamvast$reachname],each=nrow(time.table)))
  checklist$Year<-time.table$Year[match(checklist$Time,time.table$Time)]
  checklist$Key<-time.table$Key[match(checklist$Time,time.table$Time)]
  checklist$Refdate<-time.table$Refdate[match(checklist$Time,time.table$Time)]
  checklist<-cbind(checklist,time.table[match(checklist$Time,time.table$Time),Time])
  names(checklist)[ncol(checklist)]<-Time

  missing.vec<-!apply(checklist[,1:2],MARGIN = 1,FUN=function(x){return(any(out.data$Time==x[1] &
                                                                              out.data$vastid==x[2]))})

  missing.data<-data.frame(Time=checklist$Time[missing.vec],
                           Year=checklist$Year[missing.vec],
                           Day=as.integer(format(checklist$Refdate[missing.vec],"%j")),
                           Month=format(checklist$Refdate[missing.vec],"%b"),
                           Reach=checklist$Reach[missing.vec],
                           vastid=checklist$vastid[missing.vec],
                           Counts=0,
                           Density=0,
                           Area=habitat.reaches$Length[match(checklist$vastid[missing.vec],habitat.ids)],
                           Effort=habitat.reaches$Length[match(checklist$vastid[missing.vec],habitat.ids)],
                           Lat=midpoints$Y[match(checklist$vastid[missing.vec],midpoints$vastid)],
                           Lon=midpoints$X[match(checklist$vastid[missing.vec],midpoints$vastid)],
                           original=F,
                           dummy=T)

  out.data<-rbind(out.data,missing.data)
  # quick fix to months
  min.month<-min(which(month.abb%in%out.data$Month))
  max.month<-max(which(month.abb%in%out.data$Month))
  out.data$Month<-factor(out.data$Month,levels = month.abb[min.month:max.month])

  # add a bit to determine which time frames are original
  x<-aggregate(out.data$original[out.data$dummy==F],
               by=list(Time=out.data$Time[out.data$dummy==F]),FUN=max)
  time.table$Original<-F
  time.table$Original[time.table$Time%in%x$Time[which(x$x==1)]]<-T

  names(out.data)[names(out.data)=="Reach"]<-streamvast$reachname
  names(out.data)[names(out.data)=="Counts"]<-streamvast$countname

  streamvast$vastinput<-streamvast$vastinput[match(out.data$vastid,streamvast$vastinput$child_i),]
  streamvast$vastdata<-out.data
  streamvast$timetable<-time.table
  streamvast$timescale<-Time

  return(streamvast)
}


# This function adds covariates to a streamvast object
# There are two methods: If you wish to manually supply covariates, use the covariatedata argument (must match to vastdata)
# alternately, supply tables to match spatially, temporally, or both and the function will construct covariatedata for you
# formulas default toward being off, so the user must supply at least one
# the function will provide a basic default if the configurations are left blank

#' Add habitat covariates to a streamvast object
#'
#' @param streamvast A streamvast object with temporal and spatial frames defined
#' @param pform A one-sided formula for covariates affecting the probability of occurance
#' @param dform A one-sided formula for covariates affecting density given presence
#' @param pconfig A data frame of control values for pform: See ?VAST::make_data for details
#' @param dconfig A data frame of control values for dform: See ?VAST::make_data for details
#' @param spcovars A dataframe or sf object with spatially referenced covariate data
#' @param tempcovars A dataframe with temporally referenced covariate data
#' @param covariatedata A dataframe with covariates matched to vastdata
#' @param append Logical; Should the covariates be appended to any pre-existing covariate data in the object
#' @param center Logical; Should the covariates be centered to a mean value of zero
#' @param scale Logical; Should the covariates be scaled to have standard deviation of one
#'
#' @return A streamvast object with covariate features appropriately formatted
#' @export
#'
#' @examples
SetVastCovariates<-function(streamvast,pform,dform,pconfig=NULL,dconfig=NULL,spcovars,tempcovars,covariatedata,
                            append=T,center=T,scale=T){

  # checks; may wish to add more later on
  if(missing(pform) & missing(dform)){stop("Must specify at least one formula pform or dform, or else no point to covariates!")}
  if(missing(pform)){pform=stats::formula(~0)}
  if(missing(dform)){dform=stats::formula(~0)}
  if(inherits(pform,"formula")==F){pform<-stats::as.formula(pform)}
  if(inherits(dform,"formula")==F){pform<-stats::as.formula(dform)}


  if(is.null(streamvast$vastdata)){stop("Please format the time frame for the model (streamvast::SetTemporalFrame) before adding covariates!")}

  out.data<-streamvast$vastdata[,c("Time",streamvast$reachname,"vastid","Lat","Lon","original","dummy")]
  exclude.names<-c(names(out.data),"Year") # Note, we will change time to year at the end of the function

  #check for preexisting covariate data
  old.covs<-NULL
  if(is.null(streamvast$covariates$covariatedata)==F){
    if(append){
      print("Appending new covariates to pre-existing data")
      old.covs<-names(streamvast$covariates$covariatedata)
      out.data<-cbind(out.data,subset(streamvast$covariates$covariatedata,select=which(old.covs%in%exclude.names==F)))
      exclude.names<-unique(c(exclude.names,old.covs))
    }
    if(append==F){
      print("Overwriting pre-existing covariate data")
    }
  }

  # if covariatedata has been specified, transfer it over first
  if(missing(covariatedata)==F){
    if(nrow(covariatedata)!=nrow(streamvast$vastdata)){stop("Covariatedata must have the same number of rows as vastdata!")}
    if(inherits(covariatedata,"sf")){covariatedata<-sf::st_drop_geometry(covariatedata)}
    new.covs<-names(covariatedata)

    good.new.names<-which(new.covs%in%exclude.names==F)
    new.data<-subset(covariatedata,select=good.new.names)

    if(center){
      number.covs<-which(unlist(lapply(new.data,inherits,what=c("integer","numeric"))))
      if(length(number.covs)>0){
        num.data<-subset(new.data,select=number.covs)
        new.data[,number.covs]<-num.data-rep(apply(num.data,2,mean),each=nrow(new.data))
      }
    }
    if(scale){
      number.covs<-which(unlist(lapply(new.data,inherits,what=c("integer","numeric"))))
      if(length(number.covs)>0){
        num.data<-subset(new.data,select=number.covs)
        new.data[,number.covs]<-num.data/rep(sqrt(apply(num.data,2,var)),each=nrow(new.data))
      }
    }
    out.data<-cbind(out.data,new.data)
    exclude.names<-unique(c(exclude.names,new.covs))
  }

  # next we'll do spatial covariates
  if(missing(spcovars)==F){
    if(inherits(spcovars,"sf")){spcovars<-sf::st_drop_geometry(spcovars)}
    space.names<-names(spcovars)
    has.reachname<-streamvast$reachname%in%space.names
    has.vastid<-"vastid"%in%space.names
    match.name<-ifelse(has.reachname,streamvast$reachname,ifelse(has.vastid,"vastid",NA))
    if(is.na(match.name)){
      stop("Could not match spcovars to vastdata. The best way to do this is to include a column with a reach identifier is spcovars.")
    }
    good.space.names<-which(space.names%in%exclude.names==F)
    match.vec<-match(out.data[,match.name],spcovars[,match.name])
    space.data<-as.data.frame(subset(spcovars,select=good.space.names)[match.vec,])
    names(space.data)<-space.names[good.space.names]

    if(center){
      number.covs<-which(unlist(lapply(space.data,inherits,what=c("integer","numeric"))))
      if(length(number.covs)>0){
        num.data<-subset(space.data,select=number.covs)
        space.data[,number.covs]<-num.data-rep(apply(num.data,2,mean),each=nrow(space.data))
      }
    }
    if(scale){
      number.covs<-which(unlist(lapply(space.data,inherits,what=c("integer","numeric"))))
      if(length(number.covs)>0){
        num.data<-subset(space.data,select=number.covs)
        space.data[,number.covs]<-num.data/rep(sqrt(apply(num.data,2,var)),each=nrow(space.data))
      }
    }
    out.data<-cbind(out.data,space.data)
    exclude.names<-unique(c(exclude.names,space.names))
  }

  #Next we'll do temporal covariates
  if(missing(tempcovars)==F){
    if(inherits(tempcovars,"sf")){tempcovars<-sf::st_drop_geometry(tempcovars)}
    temp.names<-names(tempcovars)
    has.time<-"Time"%in%temp.names
    has.timescale<-streamvast$timescale%in%temp.names
    match.name<-ifelse(has.time,"Time",ifelse(has.timescale,streamvast$timescale,NA))
    if(is.na(match.name)){
      stop("Could not match tempcovars to vastdata. The best way to do this is to include a column with a 'Time' identifier to tempcovars")
    }
    good.time.names<-which(temp.names%in%exclude.names==F)
    match.vec<-match(out.data[,match.name],tempcovars[,match.name])
    time.data<-subset(tempcovars,select=good.time.names)[match.vec,]

    if(center){
      number.covs<-which(unlist(lapply(time.data,inherits,what=c("integer","numeric"))))
      if(length(number.covs)>0){
        num.data<-subset(time.data,select=number.covs)
        time.data[,number.covs]<-num.data-rep(apply(num.data,2,mean),each=nrow(time.data))
      }
    }
    if(scale){
      number.covs<-which(unlist(lapply(time.data,inherits,what=c("integer","numeric"))))
      if(length(number.covs)>0){
        num.data<-subset(time.data,select=number.covs)
        time.data[,number.covs]<-num.data/rep(sqrt(apply(num.data,2,var)),each=nrow(time.data))
      }
    }
    out.data<-cbind(out.data,time.data)
    exclude.names<-unique(c(exclude.names,temp.names))
  }

  names(out.data)[names(out.data)=="Time"]<-"Year"
  streamvast$covariates<-list(pformula=pform,pconfig=pconfig,dformula=dform,dconfig=dconfig,covariatedata=out.data)
  return(streamvast)
}


# This runs a number of variations of a basic vast model and returns the one with the best AIC
# laster versions may offer some choice as to which stat to use to evaluate
# also covariates are not implemented at the moment
#' Title
#'
#' @param streamvast       a streamvast object with a defined temporal frame
#' @param vastsettings     a list of settings for a vast model
#' @param optimize         logical, should the algorthim run multiple models with varying settings
#' @param maxiter          a maximum number of iterations to try before giving up
#'
#' @return A streamvast object with a fitted model attached
#' @export
#'
#' @examples
RunVAST<-function(streamvast,vastsettings,optimize=T,maxiter=3,startpar){

  if(length(streamvast$preds)>0){warning("VAST model already present. Overwriting pre-existing model.")}
  if(missing(startpar)){startpar<-NULL}

  if(is.null(vastsettings)){
    if(missing(vastsettings)){
      vastsettings <- FishStatsUtils::make_settings(Region="stream_network",
                                                    zone=10,                      # we'll need to generalize this eventually, not sure how
                                                    purpose = "index2",
                                                    fine_scale=F,
                                                    n_x = sum(streamvast$reachdata$habitat==T),
                                                    FieldConfig = c("Omega1"= 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                                                    RhoConfig = c("Beta1" = 2, "Epsilon1" = 2, "Beta2" = 2, "Epsilon2" = 2),
                                                    OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
                                                    ObsModel = c(5, 0),
                                                    bias.correct = F,
                                                    use_anisotropy = F)
      vastsettings$Method<-"Stream_network"
      vastsettings$grid_size_km<-mean(streamvast$reachdata$Length[streamvast$reachdata$habitat==T])
    }
  }else{
    optimize<-F
  }

  nmodels<-ifelse(optimize,6,1)
  if(optimize){fit.list<-list(NA,NA,NA,NA,NA,NA)}
  if(optimize==F){fit.list<-list(NA)}
  aic.vec<-rep(NA,nmodels)
  converge.vec<-rep(F,nmodels)
  settings.list<-list()

  for(i in 1:nmodels){
    # will update settings to test alternate versions
    if(i==2){
      vastsettings$FieldConfig = c("Omega1"= 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 0)
      vastsettings$RhoConfig = c("Beta1" = 2, "Epsilon1" = 2, "Beta2" = 2, "Epsilon2" = 0)
    }
    if(i==3){
      vastsettings$FieldConfig = c("Omega1"= 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
    }
    if(i==4){
      vastsettings$FieldConfig = c("Omega1"= 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
      vastsettings$RhoConfig = c("Beta1" = 2, "Epsilon1" = 2, "Beta2" = 2, "Epsilon2" = 2)
      vastsettings$ObsModel = c(7,0)
    }
    if(i==5){
      vastsettings$FieldConfig = c("Omega1"= 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 0)
      vastsettings$RhoConfig = c("Beta1" = 2, "Epsilon1" = 2, "Beta2" = 2, "Epsilon2" = 0)
    }
    if(i==6){
      vastsettings$FieldConfig = c("Omega1"= 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
    }
    settings.list[[i]]<-vastsettings

    not.converged<-T
    niter=1

    while(not.converged & niter<=maxiter){
      try(fit.list[[i]]<-FishStatsUtils::fit_model("settings" = vastsettings,
                                                   Lat_i = streamvast$vastdata$Lat,
                                                   Lon_i = streamvast$vastdata$Lon,
                                                   t_i = as.integer(streamvast$vastdata$Time),
                                                   b_i = units::as_units(streamvast$vastdata[,streamvast$countname],"count"),
                                                   a_i = units::as_units(streamvast$vastdata$Effort,"km"),
                                                   test_fit = FALSE,
                                                   PredTF_i=streamvast$vastdata$dummy,
                                                   input_grid=streamvast$vastinput,
                                                   Network_sz=streamvast$vastnetwork,
                                                   Network_sz_LL=streamvast$vastnetworkLL,
                                                   covariate_data = streamvast$covariates$covariatedata,
                                                   X1config_cp = streamvast$covariates$pconfig,
                                                   X1_formula = streamvast$covariates$pformula,
                                                   X2config_cp = streamvast$covariates$dconfig,
                                                   X2_formula = streamvast$covariates$dformula,
                                                   newtonsteps = 3,
                                                   startpar=startpar,
                                                   Aniso=F,
                                                   getsd=T,
                                                   getHessian=T,
                                                   max_cells=nrow(streamvast$vastdata),
                                                   Options = c('treat_nonencounter_as_zero' = F),
                                                   control=list(eval.max=500,iter.max=500,trace=1)))
      if(is.na(fit.list[[i]])[1]==F){
        if(length(fit.list[[i]]$Report)>1){
          not.converged<-F
          converge.vec[i]<-T
          aic.vec[i]<-fit.list[[i]]$parameter_estimates$AIC
        }else{
          fit.list[[i]]<-NA
        }
      }
      niter<-niter+1
    }
  }
  if(all(converge.vec==F)){stop("Failed to produce a successful model. Review settings and data structure.")}

  streamvast$vastmodel<-fit.list[[which.min(aic.vec)]]
  streamvast$vastsettings<-settings.list[[which.min(aic.vec)]]
  streamvast$stats$AIC<-min(aic.vec,na.rm=T)

  return(streamvast)
}


#' Get density predictions from a VAST model
#'
#' This function quickly constructs a data frame with VAST density predictions matched to
#' time and space. It also samples and provides upper and lower bounds (90%)
#'
#' @param streamvast A VAST model object
#'
#' @return A streamvast object with a dataframe of predictions added
#' @export
#'
#' @examples
VASTpreds<-function(streamvast){

  if(length(streamvast$preds)>0){warning("Prediction data already found. Overwriting pre-existing data.")}
  if(length(streamvast$eval)>0){warning("Evaluation data already found. Overwriting pre-existing data.")}
  if(length(streamvast$aucdata)>0){warning("AUC or Index data already found. Overwriting pre-existing data.")}

  #Setup
  print("Simulating predicted density distribution")
  report.dims<-dim(streamvast$vastmodel$Report$D_gct)
  habitat.reaches<-subset(streamvast$reachdata,habitat==T)

  preddata<-data.frame(Time=rep(1:report.dims[3],each=report.dims[1]),
                       vastid=rep(1:report.dims[1],times=report.dims[3]),
                       Reach=NA,
                       Length=rep(habitat.reaches$Length,times=report.dims[3]),
                       Date=rep(streamvast$timetable$Refdate,each=report.dims[1]),
                       Year=rep(streamvast$timetable$Year,each=report.dims[1]),
                       Month=rep(format(streamvast$timetable$Refdate,format="%b"),each=report.dims[1]),
                       Day=rep(streamvast$timetable$Day,each=report.dims[1]),
                       Density=NA,
                       Density_upper=NA,
                       Density_lower=NA,
                       Density_stdev=NA,
                       Count=NA,
                       Count_upper=NA,
                       Count_lower=NA,
                       Count_stdev=NA)

  preddata$Reach<-as.data.frame(habitat.reaches)[match(preddata$vastid,habitat.reaches$vastid),streamvast$reachname]

  # simulate the data for CIs
  pred.sim<-FishStatsUtils::sample_variable(Sdreport=streamvast$vastmodel$parameter_estimates$SD,
                                            Obj=streamvast$vastmodel$tmb_list$Obj,variable_name="D_gct",n_samples = 100)
  # quick loop to grab the values
  for(i in 1:report.dims[3]){
    spots<-(1+(i-1)*report.dims[1]):(i*report.dims[1])
    preddata$Density[spots]<-streamvast$vastmodel$Report$D_gct[,1,i]
    preddata$Density_lower[spots]<-apply(pred.sim[,1,i,],MARGIN = 1,FUN = stats::quantile,probs=.05,na.rm=T)
    preddata$Density_upper[spots]<-apply(pred.sim[,1,i,],MARGIN = 1,FUN = stats::quantile,probs=.95,na.rm=T)
    preddata$Density_stdev[spots]<-sqrt(apply(pred.sim[,1,i,],MARGIN = 1,FUN = stats::var))
  }

  preddata$Count<-preddata$Density*preddata$Length
  preddata$Count_upper<-preddata$Density_upper*preddata$Length
  preddata$Count_lower<-preddata$Density_lower*preddata$Length
  preddata$Count_stdev<-preddata$Density_stdev*preddata$Length

  names(preddata)[names(preddata)=="Reach"]<-streamvast$reachname
  names(preddata)[names(preddata)=="Count"]<-streamvast$countname
  names(preddata)[names(preddata)=="Count_upper"]<-paste0(streamvast$countname,"_upper")
  names(preddata)[names(preddata)=="Count_lower"]<-paste0(streamvast$countname,"_lower")
  names(preddata)[names(preddata)=="Count_stdev"]<-paste0(streamvast$countname,"_stdev")

  # now do evaldata
  print("Simulating predicted count distribution and residuals")
  eval.sim<-FishStatsUtils::sample_variable(Sdreport=streamvast$vastmodel$parameter_estimates$SD,
                                            Obj=streamvast$vastmodel$tmb_list$Obj,variable_name="D_i",n_samples = 100)
  obs.sim.data<-eval.sim[streamvast$vastdata$dummy==F,]
  obs.data<-subset(streamvast$vastdata,dummy==F)
  obs.D_i<-streamvast$vastmodel$Report$D_i[streamvast$vastdata$dummy==F]

  evaldata<-data.frame(Time=obs.data$Time,
                       Reach=obs.data[,streamvast$reachname],
                       Year=obs.data$Year,
                       Month=obs.data$Month,
                       Day=obs.data$Day,
                       Count=obs.data[,streamvast$countname],
                       Density=obs.data$Density,
                       Effort=obs.data$Effort,
                       pCount=obs.D_i,
                       pCount_upper=apply(obs.sim.data,MARGIN = 1,FUN=stats::quantile,probs=.95,na.rm=T),
                       pCount_lower=apply(obs.sim.data,MARGIN = 1,FUN=stats::quantile,probs=.05,na.rm=T),
                       pDensity=obs.D_i/obs.data$Effort,
                       pDensity_upper=apply(obs.sim.data/obs.data$Effort,MARGIN = 1,FUN=stats::quantile,probs=.95,na.rm=T),
                       pDensity_lower=apply(obs.sim.data/obs.data$Effort,MARGIN = 1,FUN=stats::quantile,probs=.05,na.rm=T),
                       Residual_count=NA,
                       Residual_density=NA,
                       Residual_pearson=NA,
                       Residual_dharma=NA)

  streamvast$stats$rhoP<-stats::cor(evaldata[6],evaldata[9],method = "pearson")
  streamvast$stats$rhoS<-stats::cor(evaldata[6],evaldata[9],method = "spearman")
  streamvast$stats$RMSE<-sqrt(sum((evaldata[6]-evaldata[9])^2)/nrow(evaldata))

  #Seems like a good place to do the residuals
  # note that due to capability issues,  dHARMa requires the vast model to be stored in an object called "fit" in
  # .GlobalEnv. Therefore, we do some back end work to make sure we don't overwrite the user's data
  # first we save a backup using a string that is unlikely to be there, then return the original from the backup
  had.fit<-F
  if(exists("fit")){
    name.nobody_would.ever_use<-fit
    had.fit<-T
  }
  fit<<-streamvast$vastmodel
  dharmaRes <-summary(streamvast$vastmodel, what="residuals", type=1)
  graphics::par(mfrow=c(1,1))
  rm(fit,envir = .GlobalEnv)
  if(had.fit){fit<<-name.nobody_would.ever_use}

  evaldata$Residual_count<-evaldata$Count-evaldata$pCount
  evaldata$Residual_density<-evaldata$Density-evaldata$pDensity
  evaldata$Residual_pearson<-(evaldata$Density-evaldata$pDensity)/sqrt(evaldata$pDensity)
  evaldata$Residual_dharma<-dharmaRes$scaledResiduals[streamvast$vastdata$dummy==F]

  names(evaldata)[names(evaldata)=="Count"]<-streamvast$countname
  names(evaldata)[names(evaldata)=="pCount"]<-paste0("p",streamvast$countname)
  names(evaldata)[names(evaldata)=="pCount_upper"]<-paste0("p",streamvast$countname,"_upper")
  names(evaldata)[names(evaldata)=="pCount_lower"]<-paste0("p",streamvast$countname,"_lower")
  names(evaldata)[names(evaldata)=="Residual_count"]<-paste0("Residual",streamvast$countname)

  # now we will need to compile the aggregate statistics, like reach averages and AUCs
  spacedata<-data.frame(Reach=as.data.frame(streamvast$reachdata)[,streamvast$reachname],Habitat=streamvast$reachdata$habitat,
                        Length=streamvast$reachdata$Length,
                        Density=NA,Count=NA,
                        pDensity=NA,pDensity_lower=NA,pDensity_upper=NA,pCount=NA,pCount_lower=NA,pCount_upper=NA,
                        tDensity=NA,tDensity_lower=NA,tDensity_upper=NA,tCount=NA,tCount_lower=NA,tCount_upper=NA)
  names(spacedata)[1]<-streamvast$reachname


  timedata<-data.frame(Time=streamvast$timetable$Time,tframe=streamvast$timetable[,streamvast$timescale],
                       Year=streamvast$timetable$Year,Month=format(streamvast$timetable$Refdate,format="%b"),
                       Day=streamvast$timetable$Day,Original=streamvast$timetable$Original,
                       Density=NA,Count=NA,
                       pDensity=NA,pDensity_lower=NA,pDensity_upper=NA,pCount=NA,pCount_lower=NA,pCount_upper=NA,
                       tDensity=NA,tDensity_lower=NA,tDensity_upper=NA,tCount=NA,tCount_lower=NA,tCount_upper=NA)

  if(tolower(streamvast$timescale)=="year"){timedata<-timedata[-3]}
  if(tolower(streamvast$timescale)=="month"){timedata<-timedata[-4]}
  if(tolower(streamvast$timescale)=="day"){timedata<-timedata[-5]}
  names(timedata)[2]<-streamvast$timescale

  # Now that time is set, we can also set up the AUC data objects
  aucdata<-data.frame(Year=rep(sort(unique(streamvast$timetable$Year)),each=nrow(habitat.reaches)),
                      Reach=rep(as.data.frame(habitat.reaches)[,streamvast$reachname],times=length(unique(streamvast$timetable$Year))),
                      AUC=NA,pAUC=NA,pAUC_lower=NA,pAUC_upper=NA,
                      tAUC=NA,tAUC_lower=NA,tAUC_upper=NA)

  auctotals<-data.frame(Year=sort(unique(streamvast$timetable$Year)),AUC=NA,pAUC=NA,pAUC_lower=NA,pAUC_upper=NA,
                        tAUC=NA,tAUC_lower=NA,tAUC_upper=NA)

  # we could also put an index in here, but it's low priority for the moment


  # Now we will start filling in spacedata
  for(i in 1:nrow(habitat.reaches)){

    # first do space
    active.reach<-as.data.frame(habitat.reaches)[i,streamvast$reachname]
    spot<-which(spacedata[,streamvast$reachname]==active.reach)

    r.simdata<-pred.sim[i,1,,]
    obs.rows<-which(obs.data[,streamvast$reachname]==active.reach)
    r.obs.data<-obs.data[obs.rows,]

    # figure out the observed densities and counts
    if(nrow(r.obs.data)>0){
      spacedata$Density[spot]<-mean(r.obs.data$Density)
      spacedata$Count[spot]<-mean(r.obs.data[,streamvast$countname])
    }else{
      spacedata$Density[spot]<-NA
      spacedata$Count[spot]<-NA
    }

    # then scaled densities and counts
    if(length(obs.rows)>0){
      r.data<-matrix(obs.sim.data[obs.rows,],nrow=length(obs.rows),ncol=100)
      r.data.density<-r.data/obs.data$Effort[obs.rows]

      spacedata$pCount[spot]<-mean(obs.D_i[obs.rows])
      spacedata$pCount_lower[spot]<-stats::quantile(apply(r.simdata,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
      spacedata$pCount_upper[spot]<-stats::quantile(apply(r.simdata,MARGIN=2,FUN=mean),prob=.95,na.rm=T)

      spacedata$pDensity[spot]<-mean(obs.D_i[obs.rows]/obs.data$Effort[obs.rows])
      spacedata$pDensity_lower[spot]<-stats::quantile(apply(r.data.density,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
      spacedata$pDensity_upper[spot]<-stats::quantile(apply(r.data.density,MARGIN=2,FUN=mean),prob=.95,na.rm=T)
    }

    # then total densities and counts
    spacedata$tDensity[spot]<-as.numeric(mean(streamvast$vastmodel$Report$D_gct[i,1,]))
    spacedata$tDensity_lower[spot]<-quantile(apply(r.data,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
    spacedata$tDensity_upper[spot]<-quantile(apply(r.data,MARGIN=2,FUN=mean),prob=.95,na.rm=T)

    spacedata$tCount[spot]<-spacedata$tDensity[spot]*spacedata$Length[spot]
    spacedata$tCount_upper[spot]<-spacedata$tDensity_upper[spot]*spacedata$Length[spot]
    spacedata$tCount_lower[spot]<-spacedata$tDensity_lower[spot]*spacedata$Length[spot]
  }
  # update the names
  names(spacedata)[which(names(spacedata)=="Count")]<-streamvast$countname
  names(spacedata)[which(names(spacedata)=="pCount")]<-paste0("p",streamvast$countname)
  names(spacedata)[which(names(spacedata)=="pCount_lower")]<-paste0("p",streamvast$countname,"_lower")
  names(spacedata)[which(names(spacedata)=="tpCount_upper")]<-paste0("p",streamvast$countname,"_upper")
  names(spacedata)[which(names(spacedata)=="tCount")]<-paste0("t",streamvast$countname)
  names(spacedata)[which(names(spacedata)=="tCount_lower")]<-paste0("t",streamvast$countname,"_lower")
  names(spacedata)[which(names(spacedata)=="tCount_upper")]<-paste0("t",streamvast$countname,"_upper")

  # Next, we'll do time
  for(i in 1:nrow(timedata)){

    active.time<-timedata$Time[i]
    active.rows<-which(obs.data$Time==active.time)

    t.data<-pred.sim[,1,i,]
    t.data.count<-t.data*habitat.reaches$Length
    t.obs.data<-obs.data[active.rows,]

    if(nrow(t.obs.data)>0){
      timedata$Density[i]<-mean(t.obs.data$Density)
      timedata$Count[i]<-mean(t.obs.data[,streamvast$countname])
    }else{
      timedata$Density[i]<-NA
      timedata$Count[i]<-NA
    }

    if(length(active.rows)>0){
      t.obs.sim.data<-matrix(obs.sim.data[active.rows,],nrow=length(active.rows),ncol=100)
      t.obs.data.density<-t.obs.sim.data/obs.data$Effort[active.rows]

      timedata$pCount[i]<-mean(obs.D_i[active.rows])
      timedata$pCount_lower[i]<-stats::quantile(apply(t.obs.sim.data,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
      timedata$pCount_upper[i]<-stats::quantile(apply(t.obs.sim.data,MARGIN=2,FUN=mean),prob=.95,na.rm=T)

      timedata$pDensity[i]<-mean(obs.D_i[active.rows]/obs.data$Effort[active.rows])
      timedata$pDensity_lower[i]<-stats::quantile(apply(t.obs.data.density,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
      timedata$pDensity_upper[i]<-stats::quantile(apply(t.obs.data.density,MARGIN=2,FUN=mean),prob=.95,na.rm=T)
    }

    timedata$tDensity[i]<-mean(streamvast$vastmodel$Report$D_gct[,1,i])
    timedata$tDensity_lower[i]<-quantile(apply(t.data,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
    timedata$tDensity_upper[i]<-quantile(apply(t.data,MARGIN=2,FUN=mean),prob=.95,na.rm=T)

    timedata$tCount[i]<-as.numeric(mean(streamvast$vastmodel$Report$D_gct[,1,i]*habitat.reaches$Length))
    timedata$tCount_lower[i]<-quantile(apply(t.data.count,MARGIN=2,FUN=mean),prob=.05,na.rm=T)
    timedata$tCount_upper[i]<-quantile(apply(t.data.count,MARGIN=2,FUN=mean),prob=.95,na.rm=T)
  }

  names(timedata)[names(timedata)=="Count"]<-streamvast$countname
  names(timedata)[names(timedata)=="pCount"]<-paste0("p",streamvast$countname)
  names(timedata)[names(timedata)=="pCount_lower"]<-paste0("p",streamvast$countname,"_lower")
  names(timedata)[names(timedata)=="pCount_upper"]<-paste0("p",streamvast$countname,"_upper")
  names(timedata)[names(timedata)=="tCount"]<-paste0("t",streamvast$countname)
  names(timedata)[names(timedata)=="tCount_lower"]<-paste0("t",streamvast$countname,"_lower")
  names(timedata)[names(timedata)=="tCount_upper"]<-paste0("t",streamvast$countname,"_upper")

  #Start filling the aucdata in here, which will require another loop
  for(y in 1:length(unique(streamvast$timetable$Year))){

    active.year<-sort(unique(streamvast$timetable$Year))[y]
    active.times<-streamvast$timetable$Time[streamvast$timetable$Year==active.year]
    day.vec<-streamvast$timetable$Day[streamvast$timetable$Time%in%active.times]
    active.rows<-which(obs.data$Time%in%active.times)

    year.total.sim<-matrix(nrow=nrow(habitat.reaches),ncol=100)
    year.total.sim.obs<-matrix(nrow=nrow(habitat.reaches),ncol=100)

    for(r in 1:nrow(habitat.reaches)){

      index<-which(aucdata$Year==active.year & aucdata$Reach==as.data.frame(habitat.reaches)[r,streamvast$reachname])
      active.id<-habitat.reaches$vastid[r]
      active.area<-habitat.reaches$Length[r]

      # Remember that pred.sim comes from D_gct and is a density; don't forget to multiple by length
      # Whereas obs.sim.data comes from D_i and is a count; no multiplication needed
      reach.preds<-streamvast$vastmodel$Report$D_gct[r,1,active.times]*active.area
      reach.obs<-obs.data[obs.data$Time%in%active.times & obs.data$vastid==active.id,]
      reach.obs.preds<-obs.D_i[obs.data$Time%in%active.times & obs.data$vastid==active.id]

      year.sim<-pred.sim[r,1,active.times,]*active.area
      year.obs.sim<-obs.sim.data[obs.data$Time%in%active.times & obs.data$vastid==active.id,]

      aucdata$AUC[index]<-MESS::auc(x = reach.obs$Day,y = reach.obs[,streamvast$countname])
      aucdata$pAUC[index]<-MESS::auc(x= reach.obs$Day,y=reach.obs.preds)
      aucdata$tAUC[index]<-MESS::auc(x = day.vec,y = reach.preds)

      year.aucs<-apply(year.sim,MARGIN = 2,FUN = function(y){return(MESS::auc(x=day.vec,y=y))})
      year.obs.aucs<-apply(year.obs.sim,MARGIN = 2,FUN = function(y){return(MESS::auc(x=reach.obs$Day,y=y))})

      aucdata$pAUC_lower[index]<-stats::quantile(year.obs.aucs,probs = .05,na.rm=T)
      aucdata$pAUC_upper[index]<-stats::quantile(year.obs.aucs,probs = .95,na.rm=T)
      aucdata$tAUC_lower[index]<-stats::quantile(year.aucs,probs = .05,na.rm=T)
      aucdata$tAUC_upper[index]<-stats::quantile(year.aucs,probs = .95,na.rm=T)

      year.total.sim[r,]<-year.aucs
      year.total.sim.obs[r,]<-year.obs.aucs
    }
    auctotals$AUC[y]<-sum(aucdata$AUC[aucdata$Year==active.year])
    auctotals$pAUC[y]<-sum(aucdata$pAUC[aucdata$Year==active.year])
    auctotals$tAUC[y]<-sum(aucdata$tAUC[aucdata$Year==active.year])

    auctotals$pAUC_lower[y]<-stats::quantile(apply(year.total.sim.obs,MARGIN = 2,FUN = sum),probs=.05,na.rm=T)
    auctotals$pAUC_upper[y]<-stats::quantile(apply(year.total.sim.obs,MARGIN = 2,FUN = sum),probs=.95,na.rm=T)
    auctotals$tAUC_lower[y]<-stats::quantile(apply(year.total.sim,MARGIN = 2,FUN = sum),probs=.05,na.rm=T)
    auctotals$tAUC_upper[y]<-stats::quantile(apply(year.total.sim,MARGIN = 2,FUN = sum),probs=.95,na.rm=T)
  }

  streamvast$preds<-preddata
  streamvast$eval<-evaldata
  streamvast$spacedata<-spacedata
  streamvast$timedata<-timedata
  streamvast$dharmafit<-dharmaRes
  streamvast$aucdata<-list(aucdata=aucdata,auctotals=auctotals)
  return(streamvast)
}


#' Makes a line plot that shows tributaries by stream distance
#'
#' @param streamvast a streamvast object with predictions or evaluations
#' @param plotvar the column name in data to be plotted
#' @param streamname the column name in the reaches with the stream name
#' @param usepreds logical, should the prediction or evaluation data be used, if absent, the function will search
#' @param title a title to display above the graph
#' @param show.names
#'
#' @return A plot of the variable with respect ot stream distance
#' @export
#'
#' @examples
plotStream<-function(streamvast,plotvar,streamname,usepreds,title,show.names="all"){

  if(missing(plotvar)){plotvar<-streamvast$countname}

  if(missing(usepreds)){
    in.preds<-plotvar%in%names(streamvast$preds)
    in.eval<-plotvar%in%names(streamvast$eval)

    if(in.preds & in.eval==F){usepreds<-T}
    if(in.preds & in.eval){
      usepreds<-T
      warnings(paste0("Column ",plotvar," found in both prediction and evaluation data. Defaulting to prediction data.
                      Specify 'usepreds=F' to use evaluation data instead."))
    }
    if(in.preds==F & in.eval){usepreds<-F}
    if(in.preds==F & in.eval==F){stop(paste0("Could not find ",plotvar," in either prediction or evaluation data"))}
  }

  if(usepreds){
    vardata<-streamvast$preds
  }else{
    vardata<-streamvast$eval
  }

  reachdata<-as.data.frame(streamvast$reachdata)

  if(show.names[1]=="all"){
    show.names<-unique(reachdata[,streamname])
    show.names[-which(show.names=="NA")]
  }

  # Note, can't use the already established parent in this case
  segdata<-data.frame(Reach=reachdata$Reach,name=reachdata[,streamname],var=0,dist=reachdata$streamdist,
                      parent=match(reachdata$from,reachdata$to),parentvar=0,parentdist=NA)
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
    ggplot2::scale_color_viridis(discrete = T)+
    ggplot2::xlab("Stream Distance (km)")+ggplot2::ylab(plotvar)

  if(missing(title)==F){outplot<-outplot+ggplot2::ggtitle(title)}

  return(outplot)
}


#' Maps, the desired data, with options for facetting
#'
#' @param streamvast a streamvast object with a prediction data frame
#' @param mapvar the column in data to map
#' @param facet a column to determine the facets, or a vector of values
#' @param background a sf object to put in the background
#' @param subset a logical vector indicating which pred values to use
#' @param make.labels should the reaches be labelled, turn off for multiple years
#' @param xaxis.breaks Vector of x values for axis lines; useful to reduce crowding
#' @param yaxis.breaks Vector of y values for axis lines; useful to reduce crowding
#' @param palette a viridis color palette for plotting
#'
#' @return A colored heatmap of the desired variable mapped over the stream network
#' @export
#'
#' @examples
plotPredictionMap<-function(streamvast,mapvar="Density",facet=NA,FUN="mean",background,subset,
                            make.labels=F,xaxis.breaks=NA, yaxis.breaks=NA,palette="turbo",
                            max=Inf){

  if(missing(subset)){subset<-rep(T,nrow(streamvast$preds))}

  if(is.na(facet)==F){
    facet.name<-"group"

    if(length(facet)==1){
      if(is.character(facet)){
        facet.name<-facet
      }
      datareach<-stats::aggregate(streamvast$preds[subset,mapvar],by=list(Reach=streamvast$preds[subset,streamvast$reachname],
                                                                    streamvast$preds[subset,facet]),FUN=FUN)
      }
    names(datareach)[2]<-facet.name
    facet.vec<-rep(sort(unique(datareach[,facet])),each=nrow(streamvast$reachdata))
  }else{
    datareach<-stats::aggregate(streamvast$preds[subset,mapvar],by=list(Reach=streamvast$preds[subset,streamvast$reachname]),FUN=FUN)
    facet.vec<-NA
    facet.name<-"group"
  }

  # add labels, gets rather complicated
  if(make.labels){
    midpoints<-as.data.frame(sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=streamvast$reachdata,sample = .5),
                                                                 crs = sf::st_crs(streamvast$reachdata))))[,1:2]
    midpoints$Reach<-as.data.frame(streamvast$reachdata)[,streamvast$reachname]
    midpoints.sf<-sf::st_as_sf(midpoints,coords=1:2,crs=sf::st_crs(streamvast$reachdata))

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

  FUN.name<-mapvar
  if(FUN=="mean"){FUN.name<-paste0("Mean ",mapvar)}
  if(FUN=="var"){FUN.name<-paste0(mapvar, " Variance")}
  if(mapvar=="Density"){FUN.name<-paste0(FUN.name,"\n (",streamvast$countname,"/km)")}
  if(mapvar=="Density_stdev"){FUN.name<-"Mean Standard\nDeviation in Density"}

  if(missing(background)){
    map<-ggplot2::ggplot()
  }else{
    map<-ggplot2::ggplot()+
      ggplot2::geom_sf(data=background,col="gray50",alpha=.5)
  }

  map<-map+
    ggplot2::geom_sf(data=reachplot,lwd=2)+
    ggplot2::geom_sf(data=reachplot,ggplot2::aes(col=preds),lwd=1.25)+ggplot2::theme_bw()+
    viridis::scale_color_viridis(option = palette,name=FUN.name,na.value="grey75")
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

  good.picks<-which(streamvast$vastdata$dummy==F)
  good.ranks<-scale.ranks[good.picks]
  good.resids<-streamvast$dharmafit$scaledResiduals[good.picks]

  bins<-0:100/100
  median.vec<-rep(NA,100)
  for(i in 1:(length(bins)-1)){
    start=bins[i]
    stop=bins[i+1]
    rank.picks<-which(good.ranks>start & good.ranks<=stop)
    median.vec[i]<-stats::median(good.resids[rank.picks])
  }

  pt.dat<-data.frame(x=good.ranks,y=good.resids)
  line.dat<-data.frame(x=1:100/100, y=stats::predict(stats::loess(y~x,data=data.frame(x=1:100/100,y=median.vec),
                                                    span = span),newdata = data.frame(x=1:100/100)))

  outplot<-ggplot2::ggplot()+ggplot2::geom_point(data=pt.dat,ggplot2::aes(x=x,y=y),alpha=.25)+
    ggplot2::geom_line(data=line.dat,ggplot2::aes(x=x,y=y),col=2,linewidth=1.5)+
    ggplot2::geom_hline(yintercept=c(.05,.25,.4,.5,.6,.75,.95),col=2,linetype=2)+
    ggplot2::theme_bw()+ggplot2::xlab("Rank-transformed Predictions")+ggplot2::ylab("Scaled Residuals")

  return(outplot)
}


