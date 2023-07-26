# This script will hold several utility functions designed to help analyze the output of VAST models

#let's make a function wrapper for all the preliminary setup
# outputs a list with the necessay components, and a diagnostic plots
#' Final formatting of Stream data for VAST
#'
#' Provides some final formatting and convenience steps before using VAST
#'
#' @param counts A data frame with count data, presumably from AssembleReddData, or a filename
#' @param reaches A sf lines object, presumably from other StreamVAST functions, or a filename
#' @param surveys An optional set of surveys, presumably from MakeSurveyTracks, or a filename
#' @param unit.conv A conversion used for distance units, can set to 1
#'
#' @return A list with several elements necessary for VAST
#' @export
#'
#' @examples
FormatStreamData<-function(counts,reaches,surveys,unit.conv=.0003048){

  # read in the data, either directly or from a file
  countdata<-counts
  if(class(counts)=="character"){countdata<-utils::read.csv(counts)}
  reachdata<-reaches
  if(class(reaches)=="character"){reachdata<-sf::st_read(reaches)}
  surveydata<-surveys
  if(class(surveys)=="character"){surveydata<-sf::st_read(surveys)}

  # Calculate some values that aren't in the data already
  reachdata$Length<-as.numeric(sf::st_length(reachdata))*unit.conv  # convert feet to km
  countdata$Area<-reachdata$Length[match(countdata$Reach,reachdata$reachid)]
  countdata$Effort<-countdata$Effort*unit.conv

  # This needs to be made more general
  if(all(c("Redds_NR","Redds_SV")%in%names(countdata))){
    countdata$Redds_Total<-countdata$Redds_NR+countdata$Redds_SV
    countdata$Density_Total<-countdata$Redds_Total/countdata$Effort
  }else{
    countdata$Density_Total<-countdata$Redds/countdata$Effort
  }

  # account for different capitalizations of day
  names(countdata)[names(countdata)%in%c("Day","day")]<-"Day"
  names(surveydata)[names(surveydata)%in%c("Day","day")]<-"Day"

  data.dates<-as.Date(paste(countdata$Year,countdata$Day,sep="_"),format="%Y_%j")
  countdata$Month<-format(data.dates,format="%b")
  countdata$Week<-ceiling(countdata$Day/7)

  # Next set up the spatial inputs. These are consistent and do not change unless you load a new data set
  # we designate the prediction and sampling knots as being at the midpoint of each reach
  midpoints<-as.data.frame(sf::st_coordinates(sf::st_line_sample(x=reachdata,sample = .5)))
  midpoints<-as.data.frame(sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=reachdata,sample = .5),
                                                               crs = "wgs84")))[,1:2]
  midpoints$Reach<-reachdata$reachid

  # add stream distance to reachdata
  reachdata$streamdist<-NA
  reachdata$streamdist[is.na(reachdata$parent)]<-reachdata$Length[is.na(reachdata$parent)]/2
  needs.doing<-reachdata$reachid[is.na(reachdata$streamdist)]
  already.done<-reachdata$reachid[is.na(reachdata$streamdist)==F]
  while(length(needs.doing>0)){
    up.next<-min(reachdata$reachid[reachdata$reachid%in%needs.doing & reachdata$parent%in%already.done])
    index<-which(reachdata$reachid==up.next)
    parent.index<-which(reachdata$reachid==reachdata$parent[index])
    reachdata$streamdist[index]<-reachdata$streamdist[parent.index]+(reachdata$Length[parent.index]/2)+
      (reachdata$Length[index]/2)

    needs.doing<-needs.doing[needs.doing!=up.next]
    already.done<-c(already.done,up.next)
  }

  # important format stuff
  # child is the same as reach id
  # parent is the reach downstream of the child
  # The root has parent = 0 and dist = Inf

  vastinput<-midpoints[order(midpoints$Reach),]
  names(vastinput)<-c("Lon","Lat","child_i")
  vastinput$Area_km2=reachdata$Length[match(vastinput$child_i,reachdata$reachid)]
  vastnetworkLL<-vastinput[,1:2]
  vastnetworkLL$child_s<-vastinput$child_i
  vastnetworkLL$parent_s<-reachdata$parent[match(vastnetworkLL$child_s,reachdata$reachid)]
  vastnetworkLL$dist_s<-reachdata$prnt_ds[match(vastinput$child_i,reachdata$reachid)]*unit.conv
  vastnetworkLL$parent_s[is.na(vastnetworkLL$parent_s)]<-0
  vastnetworkLL$dist_s[vastnetworkLL$parent_s==0]<-Inf
  vastnetwork<-vastnetworkLL[,3:5]

  # diagnostics, make sure the data is valid
  plot1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=countdata,ggplot2::aes(x=Day,y=Year),size=3,col=1)+
    ggplot2::geom_point(data=surveydata,ggplot2::aes(x=Day,y=Year),size=1,col=2)+ggplot2::theme_bw()
  print(plot1)

  plot3<-ggplot2::ggplot()+
    ggplot2::geom_sf(data=reachdata$geometry,lwd=4)+
    ggplot2::geom_sf(data=surveydata$geometry,lwd=2,col=2)+
    ggplot2::geom_sf(data=reachdata$geometry[reachdata$reachid%in%countdata$Reach],lwd=1,col=4)+ggplot2::theme_bw()
  print(plot3)

  return(list(countdata,reachdata,surveydata,vastinput,vastnetwork,vastnetworkLL))
}


#' Get density predictions from a VAST model
#'
#' This function quickly constructs a data frame with VAST density predictions matched to
#' time and space. It also samples and provides upper and lower bounds (90%)
#'
#' @param model A VAST model object
#' @param time A name for the time variable for the output data frame
#' @param space A name for the space variable for the output data frame
#' @param time.table An optional table matching time to a date, used to provide year, month, and day
#' @param sim.data If you already have a simulated data, save time by supplying it here
#'
#' @return A data frame with density predictions with uncertainty match to space and time
#' @export
#'
#' @examples
VASTpreds<-function(model,time="Time",space="Reach",time.table=NA,sim.data){

  #Setup
  report.dims<-dim(model$Report$D_gct)
  year.names<-names(model$Report$D_gct[1,1,])
  if(is.null(year.names)){year.names<-1:report.dims[3]}

  out.data<-data.frame(Time=rep(as.integer(year.names),each=report.dims[1]),
                       space=rep(1:(report.dims[1]),times=length(year.names)))
  names(out.data)<-c(time,space)

  # if supplied, expand the date information
  if(is.data.frame(time.table)){
    match.vec<-match(out.data$Time,time.table$Time)
    if(time!="Year"){out.data$Year<-time.table$Year[match.vec]}
    if(time!="Month"){out.data$Month<-as.integer(format(time.table$Refdate[match.vec],format="%m"))}
    if(time!="Day"){out.data$Day<-as.integer(format(time.table$Refdate[match.vec],format="%j"))}
  }

  out.data<-cbind(out.data,data.frame(preds=rep(NA,nrow(out.data)),lower=rep(NA,nrow(out.data)),
                                      upper=rep(NA,nrow(out.data))))
  # simulate the data for CIs
  if(missing(sim.data)){
    sim.data<-FishStatsUtils::sample_variable(Sdreport=model$parameter_estimates$SD,Obj=model$tmb_list$Obj,
                              variable_name="D_gct",n_samples = 100)
  }
  # quick loop to grab the values
  for(i in 1:report.dims[3]){
    spots<-(1+(i-1)*report.dims[1]):(i*report.dims[1])
    out.data$preds[spots]<-model$Report$D_gct[,1,i]
    out.data$lower[spots]<-apply(sim.data[,1,i,],MARGIN = 1,FUN = stats::quantile,probs=.05)
    out.data$upper[spots]<-apply(sim.data[,1,i,],MARGIN = 1,FUN = stats::quantile,probs=.95)
  }

  # covert density into redd predictions
  input.grid<-unique(model$input_args$extra_args$input_grid)

  out.data$Length<-input.grid$Area_km2[match(out.data$Reach,input.grid$child_i)]
  out.data$Redds<-out.data$preds*out.data$Length
  out.data$lowerRedds<-out.data$lower*out.data$Length
  out.data$upperRedds<-out.data$upper*out.data$Length

  return(out.data)
}


#' Get predictions matched to observations from a vast model
#'
#' @param model A VAST model object
#' @param data The model data used to fit the VAST model
#' @param effort A name for the column with effort information
#' @param time.table An optional table matching time to a date, used to provide year, month, and day
#'
#' @return A data frame of observed densities and counts matched to predictions
#' @export
#'
#' @examples
VASTeval<-function(model,data,effort="Effort",time.table=NA){

  sim.data<-FishStatsUtils::sample_variable(Sdreport=model$parameter_estimates$SD,Obj=model$tmb_list$Obj,
                            variable_name="D_i",n_samples = 100)

  out.data<-data.frame(data,
                       pDensity=model$Report$D_i/data[,effort],
                       pDensity_lower=apply(sim.data/data[,effort],MARGIN = 1,FUN=stats::quantile,probs=.05),
                       pDensity_upper=apply(sim.data/data[,effort],MARGIN = 1,FUN=stats::quantile,probs=.95),
                       pRedds=model$Report$D_i)
  out.data$pRedds_lower<-apply(sim.data,MARGIN = 1,FUN=stats::quantile,probs=.05)
  out.data$pRedds_upper<-apply(sim.data,MARGIN = 1,FUN=stats::quantile,probs=.95)

  return(out.data)
}



#' This function takes data from the Load Stream data function and formats it for use with VAST
#'
#' @param data a data frame with time, space, effort, and observations
#' @param reachdata A sf object, presumably from LoadStreamData
#' @param unit.conv a unit conversion for distances
#' @param reddname the name in data that corresponds to the Redd count
#' @param startdate dates to truncate or expand the data in yyyy-mm-dd format
#' @param enddate dates to truncate or expand the data in yyyy-mm-dd format
#' @param padzero should survey with zero observations be added at start/end dates
#' @param Time character giving the timescale to use
#'
#' @return
#' @export
#'
#' @examples
MakeVASTData<-function(data, reachdata,unit.conv=1, reddname="Redds",
                       startdate=NA,enddate=NA,padzero=T, Time="Year"){

  if(Time=="Year"){timeformat<-"%Y"}
  if(Time=="Month"){timeformat<-"%Y-%m"}
  if(Time%in%c("Biweek","Week","Day")){timeformat<-"%Y-%j"}

  tdivide<-ifelse(Time=="Biweek",14,ifelse(Time=="Week",7,1))

  #Get the dates ; fun note, "ifelse" will screw up dates
  data.dates<-as.Date(paste(data$Year,data$Day,sep="-"),format="%Y-%j")

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
  data<-data[data.dates>=startdate & data.dates<=enddate,]

  # Now need to go through and get tables that convert dates into sequential time steps
  minyear<-as.integer(format(startdate,"%Y"))
  maxyear<-as.integer(format(enddate,"%Y"))
  nyears<-length(minyear:maxyear)
  if(Time=="Year"){
    time.table<-data.frame(Time=1:nyears,Year=minyear:maxyear,Key=minyear:maxyear,
                           Refdate=as.Date(paste(minyear:maxyear,round(mean(data$Day)),sep="-"),format="%Y-%j"))
    data$Key<-data$Year
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

    if(Time=="Month"){time.table$Refdate<-as.Date(paste(time.table$Year,time.table[,Time],15,sep="-"),format="%Y-%m-%d")}
    if(Time=="Biweek"){
      days<-ifelse((time.table[,Time]-1)*14+7>365,365,(time.table[,Time]-1)*14+7)
      time.table$Refdate<-as.Date(paste(time.table$Year,days,sep="-"),format="%Y-%j")
    }
    if(Time=="Week"){
      days<-ifelse((time.table[,Time]-1)*7+4>365,365,(time.table[,Time]-1)*7+4)
      time.table$Refdate<-as.Date(paste(time.table$Year,days,sep="-"),format="%Y-%j")
    }
    if(Time=="Day"){time.table$Refdate<-as.Date(paste(time.table$Year,time.table[,Time],sep="-"),format="%Y-%j")}

    data$Key<-paste(data$Year,ceiling(as.integer(format(data.dates,format=strsplit(timeformat,split="-")[[1]][2]))/tdivide),sep="-")
  }


  # Set up some reach info
  midpoints<-as.data.frame(sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=reachdata,sample = .5),crs = "wgs84")))[,1:2]
  midpoints$Reach<-reachdata$reachid

  # Get the preliminary output object that has just the original data
  out.data<-data.frame(Time=time.table$Time[match(data$Key,time.table$Key)],
                       Year=data$Year,
                       Day=data$Day,
                       Month=data$Month,
                       Reach=data$Reach,
                       Redds=data[,reddname],
                       Density=data[,reddname]/data$Effort,
                       Area=data$Area,
                       Effort=data$Effort,
                       Lat=midpoints[match(x=data$Reach,midpoints$Reach),2],
                       Lon=midpoints[match(x=data$Reach,midpoints$Reach),1],
                       original=T,
                       dummy=F)

  # Now make padding data, if adding zeros to the ends of each season
  if(padzero){
    if(Time=="Year"){stop("Padding with zeros doesn't make sense with an annual time scale")}

    # lot's of tricksy formatting with dates
    frontpad<-as.Date(paste(minyear:maxyear,format(startdate,"%j"),sep="-"),format="%Y-%j")
    frontpad.time0<-ceiling(as.integer(strsplit(format(frontpad,timeformat),split="-")[[1]][2])/tdivide)
    frontpad.time<-time.table$Time[match(paste((minyear:maxyear),frontpad.time0,sep="-"),time.table$Key)]

    backpad<-as.Date(paste(minyear:maxyear,format(enddate,"%j"),sep="-"),format="%Y-%j")
    backpad.time0<-ceiling(as.integer(strsplit(format(backpad,timeformat),split="-")[[1]][2])/tdivide)
    backpad.time<-time.table$Time[match(paste((minyear:maxyear),backpad.time0,sep="-"),time.table$Key)]

    for(y in 1:nyears){
      yeardata<-subset(out.data,Year==((minyear:maxyear)[y]))
      yeardata$Date<-as.Date(paste(yeardata$Year,yeardata$Day,sep="-"),format="%Y-%j")

      if(nrow(yeardata)>0){
        front.check<-stats::aggregate(yeardata$Date,by=list(yeardata$Reach),FUN=function(x){return(min(x)<=frontpad[y])})
        back.check<-stats::aggregate(yeardata$Date,by=list(yeardata$Reach),FUN=function(x){return(min(x)>=backpad[y])})

        not.sampled<-which(1:nrow(reachdata)%in%yeardata$Reach==F)
        needs.front<-c(not.sampled,front.check$Group.1[front.check$x==F])
        needs.back<-c(not.sampled,back.check$Group.1[front.check$x==F])
      }else{
        needs.front<-reachdata$reachid
        needs.back<-reachdata$reachid
      }
      if(length(needs.front)>0){
        front.data<-data.frame(Time=frontpad.time[y],
                               Year=(minyear:maxyear)[y],
                               Day=as.integer(format(time.table$Refdate[match(frontpad.time[y],time.table$Time)],"%j")),
                               Month=format(time.table$Refdate[match(frontpad.time[y],time.table$Time)],"%b"),
                               Reach=needs.front,
                               Redds=0,
                               Density=0,
                               Area=as.numeric(reachdata$Length[match(needs.front,reachdata$reachid)]),
                               Effort=as.numeric(reachdata$Length[match(needs.front,reachdata$reachid)]),
                               Lat=midpoints[match(needs.front,reachdata$reachid),2],
                               Lon=midpoints[match(needs.front,reachdata$reachid),1],
                               original=F,
                               dummy=F)
        out.data<-rbind(out.data,front.data)
      }
      if(length(needs.back)>0){
        back.data<-data.frame(Time=backpad.time[y],
                              Year=(minyear:maxyear)[y],
                              Day=as.integer(format(time.table$Refdate[match(backpad.time[y],time.table$Time)],"%j")),
                              Month=format(time.table$Refdate[match(backpad.time[y],time.table$Time)],"%b"),
                              Reach=needs.back,
                              Redds=0,
                              Density=0,
                              Area=as.numeric(reachdata$Length[match(needs.back,reachdata$reachid)]),
                              Effort=as.numeric(reachdata$Length[match(needs.back,reachdata$reachid)]),
                              Lat=midpoints[match(needs.back,reachdata$reachid),2],
                              Lon=midpoints[match(needs.back,reachdata$reachid),1],
                              original=F,
                              dummy=F)
        out.data<-rbind(out.data,back.data)
      }
    }
  }

  # Now to make the dummy data that signals VAST to interpolate a prediction
  checklist<-data.frame(Time=rep(time.table$Time,nrow(reachdata)),Reach=rep(1:nrow(reachdata),each=nrow(time.table)))
  checklist$Year<-time.table$Year[match(checklist$Time,time.table$Time)]
  checklist$Key<-time.table$Key[match(checklist$Time,time.table$Time)]
  checklist$Refdate<-time.table$Refdate[match(checklist$Time,time.table$Time)]
  checklist<-cbind(checklist,time.table[match(checklist$Time,time.table$Time),Time])
  names(checklist)[ncol(checklist)]<-Time

  missing.vec<-!apply(checklist[,1:2],MARGIN = 1,FUN=function(x){return(any(out.data$Time==x[1] &
                                                                              out.data$Reach==x[2]))})

  missing.data<-data.frame(Time=checklist$Time[missing.vec],
                           Year=checklist$Year[missing.vec],
                           Day=as.integer(format(checklist$Refdate[missing.vec],"%j")),
                           Month=format(checklist$Refdate[missing.vec],"%b"),
                           Reach=checklist$Reach[missing.vec],
                           Redds=0,
                           Density=0,
                           Area=reachdata$Length[match(checklist$Reach[missing.vec],reachdata$reachid)],
                           Effort=reachdata$Length[match(checklist$Reach[missing.vec],reachdata$reachid)],
                           Lat=midpoints$Y[match(checklist$Reach[missing.vec],midpoints$Reach)],
                           Lon=midpoints$X[match(checklist$Reach[missing.vec],midpoints$Reach)],
                           original=F,
                           dummy=T)


  out.data<-rbind(out.data,missing.data)
  # quick fix to months
  min.month<-min(which(month.abb%in%out.data$Month))
  max.month<-max(which(month.abb%in%out.data$Month))
  out.data$Month<-factor(out.data$Month,levels = month.abb[min.month:max.month])

  return(list(out.data,time.table))
}


#' Maps, the desired data, with options for facetting
#'
#' @param data a data frame with predictions and any grouping data
#' @param mapvar the column in data to map
#' @param reaches an sf object with reach data
#' @param facet a column to determine the facets, or a vector of values
#' @param reach.names the column in reaches that holds the id numbers
#' @param make.labels should the reaches be labelled, turn off for multiple years
#' @param xaxis.breaks Vector of x values for axis lines; useful to reduce crowding
#' @param yaxis.breaks Vector of y values for axis lines; useful to reduce crowding
#' @param palette a viridis color palette for plotting
#'
#' @return A colored heatmap of the desired variable mapped over the stream network
#' @export
#'
#' @examples
plotPredictionMap<-function(data, mapvar="preds", reaches,facet=NA,reach.names="reachid",
                            make.labels=T,xaxis.breaks=NA, yaxis.breaks=NA,palette="turbo"){

  if(is.na(facet)==F){
    facet.name<-"group"

    if(length(facet)==1){
      if(is.character(facet)){
        facet.name<-facet
      }
      datareach<-stats::aggregate(data[,mapvar],by=list(Reach=data$Reach,data[,facet]),FUN=mean)
    }
    if(length(facet)==nrow(data)){
      datareach<-stats::aggregate(data[,mapvar],by=list(Reach=data$Reach,facet),FUN=mean)
    }
    names(datareach)[2]<-facet.name
    facet.vec<-rep(sort(unique(datareach[,facet])),each=nrow(reaches))
  }else{
    datareach<-stats::aggregate(data[,mapvar],by=list(Reach=data$Reach),FUN=mean)
    facet.vec<-NA
    facet.name<-"group"
  }

  # add labels, gets rather complicated
  if(make.labels){

    midpoints<-as.data.frame(sf::st_coordinates(sf::st_transform(sf::st_line_sample(x=reaches,sample = .5),crs = sf::st_crs(reaches))))[,1:2]
    midpoints$Reach<-reaches[reach.names]
    midpoints.sf<-sf::st_as_sf(midpoints,coords=1:2,crs=sf::st_crs(reaches))

    labeldat<-data.frame(X=rep(NA,nrow(reaches)),Y=NA,labels=reaches$reachid)
    meanlength<-mean(sf::st_length(reaches))
    distmat<-sf::st_distance(reaches)
    for(i in 1:nrow(reaches)){

      #grab all nearby features, combine to a multi-line, compute the buffer, and pick a point on the buffer
      nearby<-which(distmat[i,]<meanlength*1.5)
      nearby.line<-sf::st_as_sf(sf::st_combine(reaches[nearby,]))

      buf<-sf::st_cast(sf::st_buffer(nearby.line,dist=meanlength/2),"LINESTRING")
      nearest.point<-sf::st_as_sf(sf::st_nearest_points(y = midpoints.sf[i,],x=buf))

      labeldat$X[i]<-sf::st_coordinates(nearest.point)[1,1]
      labeldat$Y[i]<-sf::st_coordinates(nearest.point)[1,2]
    }
  }

  # Set up the ggplot data
  reachplot<-sf::st_as_sf(rep(reaches$geometry,times=length(facet.vec)))
  reachplot$reach<-rep(reaches$reachid,times=length(facet.vec))
  reachplot$group<-facet.vec
  reachplot$preds<-NA
  for(i in 1:nrow(reachplot)){
    if(length(facet.vec)>1){
      reachplot$preds[i]<-datareach$x[datareach$Reach==reachplot$reach[i] &
                                       datareach[,facet.name]==reachplot$group[i]]
    }else{
      reachplot$preds[i]<-datareach$x[datareach$Reach==reachplot$reach[i]]
    }
  }

  bounds<-sf::st_bbox(sf::st_transform(reachplot,crs = "wgs84"))

  map<-ggplot2::ggplot()+
    ggplot2::geom_sf(data=reachplot,lwd=3)+
    ggplot2::geom_sf(data=reachplot,ggplot2::aes(col=preds),lwd=2)+ggplot2::theme_bw()+
    viridis::scale_color_viridis(option = palette,name=mapvar)
  if(length(facet.vec)>1){
    map<-map+ggplot2::facet_wrap(~group)
  }
  if(is.na(xaxis.breaks[1])==F){
    map<-map+ggplot2::scale_x_continuous(breaks = xaxis.breaks)
  }
  if(is.na(yaxis.breaks[1])==F){
    map<-map+ggplot2::scale_y_continuous(breaks = yaxis.breaks)
  }
  if(make.labels){
    map<-map+ggplot2::geom_text(data=labeldat,ggplot2::aes(x=X,y=Y,label=labels),check_overlap = T)
  }
  return(map)
}


#' Makes a line plot that shows tributaries by stream distance
#'
#' @param data a data frame with a variable to plot matched to reach
#' @param plotvar the column name in data to be plotted
#' @param reaches a sf object containing the reach data
#' @param streamdist name of a column in reaches with stream distance
#' @param rnames name or id of reaches
#' @param main the name of the mainstem or most important stream, which will be highlighted
#'
#' @return A plot of the variable with respect ot stream distance
#' @export
#'
#' @examples
plotStream<-function(data,plotvar,streamdist="streamdist",rnames="STRM_NA",reaches,main=NA){

  pointdata<-stats::aggregate(data[,plotvar],by=list(Reach=preds$Reach),FUN=mean)
  pointdata$streamdist<-as.data.frame(reaches)$streamdist[match(pointdata$Reach,reaches$reachid)]
  pointdata$name<-as.data.frame(reaches)[match(pointdata$Reach,reaches$reachid),rnames]
  names(pointdata)[2]<-"Density"

  parent.table<-data.frame(stream=unique(reaches$STRM_NA),
                           parent=as.integer(sapply(unique(reaches$STRM_NA),
                                                    FUN=function(x){min(reaches$parent[reaches$STRM_NA==x])})))

  joiners<-data.frame(Reach=parent.table$parent,
                      Density=pointdata$Density[match(parent.table$parent,pointdata$Reach)],
                      streamdist=pointdata$streamdist[match(parent.table$parent,pointdata$Reach)],
                      name=parent.table$stream)

  linedata<-rbind(pointdata,stats::na.omit(joiners))

  if(is.na(main)==F){
    pointdata$main<-pointdata$name==main
    linedata$main<-linedata$name==main

    outplot<-ggplot2::ggplot()+
      ggplot2::geom_point(data=pointdata,ggplot2::aes(x=streamdist,y=Density,col=name,size=main),shape=16)+
      ggplot2::geom_line(data=linedata,ggplot2::aes(x=streamdist,y=Density,col=name,group=name,lwd=main))+
      ggplot2::scale_size_manual(values=c(1.5,3),guide="none")+
      ggplot2::scale_linewidth_manual(values=c(1,1.6),guide="none")+
      ggplot2::theme_bw()+ggplot2::theme(legend.title = ggplot2::element_blank())+
      ggplot2::xlab("Stream Distance (km)")+ggplot2::ggtitle(main)
    outplot

  }else{
    outplot<-ggplot2::ggplot()+
      ggplot2::geom_point(data=pointdata,ggplot2::aes(x=streamdist,y=Density,col=name),shape=16)+
      ggplot2::geom_line(data=linedata,ggplot2::aes(x=streamdist,y=Density,col=name,group=name))+
      ggplot2::theme_bw()+ggplot2::theme(legend.title = ggplot2::element_blank())+
      ggplot2::xlab("Stream Distance (km)")
  }

  return(outplot)
}


#' Assemble estimates of total population index by time
#'
#' @param model A fitted VAST model
#' @param eval A dataframe produced by VASTeval
#'
#' @return A data frame with the population index and CIs
#' @export
#'
#' @examples
VASTtotals<-function(model,eval){

  out.data<-stats::aggregate(eval[,"Redds"],by=list(Time=eval$Time),FUN=sum)
  names(out.data)[2]<-"Redds"

  sim.data<-FishStatsUtils::sample_variable(Sdreport=model$parameter_estimates$SD,Obj=model$tmb_list$Obj,
                            variable_name="Index_ctl",n_samples = 100)

  out.data$pRedds<-as.numeric(model$Report$Index_ctl[1,out.data$Time,1])
  out.data$pRedds_lower<-apply(sim.data[1,out.data$Time,1,],MARGIN = 1,FUN=stats::quantile,probs = .05)
  out.data$pRedds_upper<-apply(sim.data[1,out.data$Time,1,],MARGIN = 1,FUN=stats::quantile,probs = .95)

  return(out.data)
}
