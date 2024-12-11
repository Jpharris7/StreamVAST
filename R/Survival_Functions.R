#' Convert redd data into survival data set
#'
#' This function takes a dataframe of repeated redd observations over time and converts them into
#' a format suitable for modeling redd survival.
#'
#' @param streamvast    A streamvast object with reachdata, countdata, and surveydata
#' @param redds         a dataframe containing observations of redds in each row
#' @param redd.ids      character; the name of the column with consistent id for each read over time
#' @param redd.crs      a crs for the redds
#' @param redd.coords   A vector of length 2 with the names of the columns for longitude and latitude
#' @param survey.tol    numeric; a max tolerance for matching redds to survey lines
#' @param redd.status   character; the name of the column with redd status (should use NR, SV, NV coding)
#' @param min.surveys   A minimum number of surveys for observations of a reach to be included
#' @param buffer        A number of days to be added to redds found on the first or last day of a season
#'
#' @return a dataframe formated for use with INLA survival models
#' @export
#'
#' @examples
MakeReddSurvival<-function(streamvast,redds,redd.ids="redd_name_txt",redd.crs="wgs84",
                           redd.coords=c("lon","lat"),survey.tol=250,
                           redd.status="redd_status_code",min.surveys=1,buffer=14){

  reach.box<-sf::st_bbox(sf::st_transform(streamvast$reachdata,crs=redd.crs))
  good.redds<-redds[redds[,redd.coords[1]]>=reach.box[1] & redds[,redd.coords[1]]<=reach.box[3] &
                      redds[,redd.coords[2]]>=reach.box[2] & redds[,redd.coords[2]]<=reach.box[4],]
  redd.names<-unique(good.redds[,redd.ids])

  survival.data<-data.frame(Redd=redd.names,Year=NA,min.start=NA,max.start=NA,min.end=NA,max.end=NA,
                            min.duration=NA,max.duration=NA,Reach=NA,lon=NA,lat=NA,complete=NA)
  surveys<-streamvast$surveydata
  names(surveys)[names(surveys)=="Day"]<-"day"

  # We will also screen out any reaches that receieve a single survey that year
  coverage.table<-stats::aggregate(streamvast$countdata$Year,by=list(Year=streamvast$countdata$Year,
                                                                     Reach=streamvast$countdata[,streamvast$reachname]),
                                   FUN=length)
  good.coverage<-coverage.table[coverage.table$x>min.surveys,]

  pb<-utils::txtProgressBar(min = 0, max = length(redd.names), initial = 0, style=3)
  for(i in 1:length(redd.names)){

    r.data<-redds[which(redds[,redd.ids]==redd.names[i]),]

    survival.data$Year[i]<-r.data$Year[1]
    survival.data$lon[i]<-r.data[1,redd.coords[1]]
    survival.data$lat[i]<-r.data[1,redd.coords[2]]

    # check that the redd contains a valid history
    if(sum(r.data[,redd.status]=="NR")==1 & sum(r.data[,redd.status]=="NV")<=1 &
       is.na(r.data[1,redd.coords[1]])==F){

      redd.sf<-sf::st_as_sf(as.data.frame(survival.data[i,]),coords=c("lon","lat"),crs=redd.crs)
      redd.sf2<-sf::st_transform(redd.sf,crs=sf::st_crs(streamvast$reachdata))
      closest<-sf::st_nearest_feature(redd.sf2,streamvast$reachdata)
      survival.data$Reach[i]<-as.data.frame(streamvast$reachdata)[closest,streamvast$reachname]

      #check that the reach was covered enough to be useful
      is.good<-any(good.coverage$Year==r.data$Year[1] & good.coverage$Reach==closest)

      if(is.good){
        # let's fix it to use the actual survey tracks
        y.surveys<-surveys[which(surveys$Year==r.data$Year[1]),]
        dists<-as.numeric(sf::st_length(sf::st_nearest_points(redd.sf2,y.surveys)))
        best.dist<-min(dists)

        # Check that the dist is close enough to be valid
        if(best.dist<=survey.tol){
          r.surveys<-y.surveys[dists==best.dist,]
          day.vec<-sort(r.surveys$day)

          # the max start is the day the redd is observed
          survival.data$max.start[i]<-r.data$day[r.data[,redd.status]=="NR"]

          # the min start is the day after the previous survey, or buffer a few days if there is none
          found.index<-which(day.vec==survival.data$max.start[i])
          if(length(found.index)==1 & nrow(r.surveys)>1){

            survival.data$complete[i]<-TRUE

            if(found.index==1){
              survival.data$min.start[i]<-survival.data$max.start[i]-buffer
              survival.data$complete[i]<-FALSE
            }else{
              survival.data$min.start[i]<-day.vec[found.index-1]+1
            }

            #the min end is the day after the last record where it was visible
            last.visible.day<-max(r.data$day[r.data[,redd.status]%in%c("NR","SV")])
            survival.data$min.end[i]<-last.visible.day+1

            # the max end is the day before it was recorded as NV or NO, after the last visible day,
            # or 14 after the last still-visible record
            not.visible.days<-r.data$day[r.data[,redd.status]%in%c("NV","NO")]
            potential.last.days<-not.visible.days[not.visible.days>last.visible.day]

            if(length(potential.last.days)>0){
              survival.data$max.end[i]<-min(potential.last.days)-1
            }else{
              survival.data$max.end[i]<-last.visible.day+buffer
              survival.data$complete[i]<-FALSE
            }
          }
        }
      }
    }
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)
  survival.data$min.duration<-survival.data$min.end-survival.data$max.start
  survival.data$max.duration<-survival.data$max.end-survival.data$min.start

  return(stats::na.omit(survival.data))
}




#' Make an adjacency matrix from an sf object
#'
#' @param reaches a sf object with LINESTRING geometry representing a stream network
#'
#' @return an adjacency matrix used by INLA for spatial relationships
#' @export
#'
#' @examples
AdjacencyMatrix<-function(reaches){

  out.mat<-matrix(data=0,nrow=nrow(reaches),ncol=nrow(reaches))

  for(i in 1:nrow(reaches)){
    parent<-reaches$parent[i]
    children<-reaches$Reach[which(reaches$parent==reaches$Reach[i])]
    if(is.na(parent)==F){out.mat[i,parent]<-1}
    if(length(children)>0){out.mat[i,children]<-1}
  }
  return(out.mat)
}

#' Summarize the results of a survival model in tabular form
#'
#' @param model A fitted INLA survival model
#' @param extrayears any years not present in the data for which an estimate is desired
#'
#' @return an adjacency matrix used by INLA for spatial relationships
#' @export
#'
#' @examples
SurvivalTable<-function(model,extrayears){
  if(inherits(model,"inla")==F){stop("Requires a fitted inla model")}
  fam<-model$.args$family
  if(fam%in%c("gammasurv","weibullsurv")==F){
    stop("Sorry, the current version only supports gammasurv and weibull surv distributions")
  }

  shape<-model$summary.hyperpar[1,1]

  surv.table<-data.frame(Year=rep(model$summary.random$Year$ID,each=nrow(model$summary.random$Reach)),
                         Reach=rep(model$summary.random$Reach$ID,times=nrow(model$summary.random$Year)),
                         Shape=shape,Scale=NA,MedLife=NA)

  for(i in 1:nrow(surv.table)){

    year.eff<-model$summary.random$Year$mean[model$summary.random$Year$ID==surv.table$Year[i]]
    reach.eff<-model$summary.random$Reach$mean[model$summary.random$Reach$ID==surv.table$Reach[i]]

    if(fam=="gammasurv"){
      surv.table$Scale[i]<-exp(model$summary.fixed[1,1]+year.eff+reach.eff)/shape
      surv.table$MedLife[i]<-stats::qgamma(.5,shape = shape,scale = surv.table$Scale[i])
    }
    if(fam=="weibullsurv"){
      surv.table$Scale[i]<-(1/(model$summary.fixed[1,1]+year.eff+reach.eff))^(1/shape)
      surv.table$MedLife[i]<-stats::qweibull(.5,shape = shape,scale = surv.table$Scale[i])
    }
  }

  # If necessary, extrapolate to additional years, as an average of adjacent years
  if(missing(extrayears)==F){

    #screen out any years already in the data
    extrayears2<-extrayears[extrayears%in%surv.table$Year==F]
    if(length(extrayears2>0)){
      extra.table<-data.frame(Year=rep(extrayears2,each=nrow(model$summary.random$Reach)),
                              Reach=rep(model$summary.random$Reach$ID,times=length(extrayears2)),
                              Shape=shape,Scale=NA,MedLife=NA)

      datayears<-sort(unique(surv.table$Year))
      for(i in 1:nrow(extra.table)){

        closest.left<-max(datayears[datayears<extra.table$Year[i]])
        closest.right<-min(datayears[datayears>extra.table$Year[i]])
        extra.table$Scale[i]<-mean(surv.table$Scale[surv.table$Year%in%c(closest.left,closest.right) &
                                                      surv.table$Reach==extra.table$Reach[i]])

        if(fam=="gammasurv"){extra.table$MedLife[i]<-stats::qgamma(.5,shape = shape,scale = extra.table$Scale[i])}
        if(fam=="weibullsurv"){extra.table$MedLife[i]<-stats::qweibull(.5,shape = shape,scale = extra.table$Scale[i])}
      }
      surv.table<-rbind(surv.table,extra.table)
    }else{
      warning("All values for extrayears already present in model.")
    }
  }
  surv.table<-surv.table[order(surv.table$Year,surv.table$Reach),]

  return(surv.table)
}


#' Plot a histogram of survival estimates across years and reaches
#'
#' @param table A table of survival estimates, usually from'SurvivalTable'
#' @param year numeric; a subset of years to include in the histogram
#' @param reach numeric; a subset of reach numbers to include in the histogram
#' @param title character; a title to pass to ggtitle()
#'
#' @return a ggplot object for the histogram
#' @export
#'
#' @examples
plotSurvivalHistogram<-function(table,year="all",reach="all",title){

  if(year[1]!="all"){table<-subset(table,Year%in%year)}
  if(reach[1]!="all"){table<-subset(table,Reach%in%reach)}

  outplot<-ggplot2::ggplot()+
    ggplot2::geom_histogram(data = table,ggplot2::aes(MedLife),bins=20)+
    ggplot2::geom_vline(xintercept=stats::quantile(table$MedLife,probs=c(.05,.5,.95)),
                        col=2,linetype=c(3,2,3),linewidth=c(1,1.25,1))+ggplot2::xlab("Days")+
    ggplot2::theme_bw()
  if(missing(title)==F){
    outplot<-outplot+ggplot2::ggtitle(title)
  }
  return(outplot)
}


#' Plot survival curves from a INLA survival model
#'
#' @param model a fitted INLA survival model
#' @param data a dataframe of redd observations, usually from MakeReddSurvival
#' @param year numeric; a subset of years to include in the plot
#' @param mult a numeric multiplier, used to backtransform values, if necessary
#' @param reach numeric; a subset of reach numbers to include in the plot
#' @param title character; a title to pass to ggtitle()
#'
#' @return a ggplot showing the
#' @export
#'
#' @examples
plotSurvivalCurves<-function(model,data,year="all",reach="all",title,mult=1){

  surv.table<-SurvivalTable(model)
  surv.table$MedLife<-surv.table$MedLife*mult

  gg.seg.table<-data.frame(Year=NA,Reach=NA, group=rep(1:nrow(surv.table),each=100),
                           Day=rep(1:100,times=nrow(surv.table)),Surv=NA)
  fam<-model$.args$family
  xseq<-1:100

  for(i in 1:nrow(surv.table)){
    start<-1+(100*(i-1))
    end<-100*i
    gg.seg.table$Year[start:end]<-surv.table$Year[i]
    gg.seg.table$Reach[start:end]<-surv.table$Reach[i]
    if(fam=="gammasurv"){
      gg.seg.table$Surv[start:end]<-1-stats::pgamma(xseq[1:100]/mult,shape=model$summary.hyperpar[1,1],
                                                    scale=surv.table$Scale[i])
    }
    if(fam=="weibullsurv"){
      gg.seg.table$Surv[start:end]<-1-stats::pweibull(xseq[1:100]/mult,shape=model$summary.hyperpar[1,1],
                                                      scale=surv.table$Scale[i])
    }
  }
  if(year[1]!="all"){gg.seg.table<-subset(gg.seg.table,Year%in%year)}
  if(reach[1]!="all"){gg.seg.table<-subset(gg.seg.table,Reach%in%reach)}

  scale<-exp(model$summary.fixed[1,1])/model$summary.hyperpar[1,1]

  out.plot<-ggplot2::ggplot()+ggplot2::theme_bw()+
    ggplot2::geom_line(data=gg.seg.table,ggplot2::aes(x=Day,y=Surv,group=group),col=2,alpha=.05)+
    ggplot2::geom_point(data=data.frame(Day=stats::quantile(data$mean.duration,probs=0:100/100),Surv=1-0:100/100),
                        ggplot2::aes(x=Day,y=Surv))+
    ggplot2::geom_errorbarh(data=data.frame(Surv=1-0:100/100,xmin=stats::quantile(data$min.duration,probs=0:100/100),
                                            xmax=stats::quantile(data$max.duration,probs=0:100/100)),
                            ggplot2::aes(y=Surv,xmin=xmin,xmax=xmax),alpha=.5,col="gray50")
  if(fam=="gammasurv"){
    out.plot<-out.plot+
      ggplot2::geom_line(data=data.frame(Day=stats::qgamma((0:100)/100,shape = model$summary.hyperpar[1,1],
                                                           scale = scale),Surv=1-0:100/100),
                         ggplot2::aes(x=Day*mult,y=Surv),linewidth=1)+ggplot2::ylab("Survival")
  }
  if(fam=="weibullsurv"){
    out.plot<-out.plot+
      ggplot2::geom_line(data=data.frame(Day=stats::qweibull((0:100)/100,shape = model$summary.hyperpar[1,1],
                                                             scale = scale),Surv=1-0:100/100),
                         ggplot2::aes(x=Day*mult,y=Surv),linewidth=1)+ggplot2::ylab("Survival")
  }
  if(missing(title)==F){
    out.plot<-out.plot+ggplot2::ggtitle(title)
  }
  return(out.plot)
}


#' Make a table of escapement values based on the AUC method
#'
#' @param streamvast a fitted streamvast object with auc predictions
#' @param survival a table of survival values, usually from SurvivalTable
#' @param fixed.survival numeric; sets survival to a fixed value for all years/reaches
#' @param mult numeric: a value to multiply the number of unique redds by, typically a spawner/redds ratio
#'
#' @return a table of escapement values by year and reach
#' @export
#'
#' @examples
MakeEscapement<-function(streamvast,survival,fixed.survival=F,mult=1.62,years="all",reaches="all"){

  if(is.null(streamvast$escapedata)==F){warning("Overwriting previous escapement data")}

  escape<-streamvast$aucdata$aucdata
  escape$Redds<-NA
  escape$Redds_lower<-NA
  escape$Redds_upper<-NA

  reachname<-streamvast$reachname

  if(years[1]=="all"){
    good.years<-sort(unique(escape$Runyear))
  }else{
    good.years<-sort(unique(years[years%in%escape$Runyear]))
  }
  if(reaches[1]=="all"){
    good.reaches<-sort(unique(escape[,reachname]))
  }else{
    good.reaches<-sort(unique(reaches[reaches%in%escape[,reachname]]))
  }

  good.escape<-escape[escape$Runyear%in%good.years & escape[,reachname]%in%good.reaches,]

  if(fixed.survival){
    if(is.numeric(survival)==F | length(survival)>1){stop(" If using 'fixed.survival = TRUE',
                 then 'survival' should be a single number representing the average redd life in days")}
    good.escape$pred_Redds<-good.escape$pred_AUC/survival
    good.escape$Redds_SD<-good.escape$pred_AUC/survival
    good.escape$pred_Redds_lower<-good.escape$pred_AUC_lower/survival
    good.escape$pred_Redds_upper<-good.escape$pred_AUC_upper/survival

    yvec<-sort(unique(streamvast$timetable$Runyear))
    rvec<-sort(unique(as.data.frame(streamvast$reachdata)[,streamvast$reachname]))
    survival<-data.frame(Year=rep(yvec,each=length(rvec)),Reach=rep(rvec,length(yvec)),MedLife=survival)

  }else{
    for(i in 1:nrow(good.escape)){
      survival.match<-which(survival$Year==good.escape$Runyear[i] &
                              survival[,streamvast$reachname]==good.escape[i,streamvast$reachname])
      good.escape$Redds[i]<-good.escape$pred_AUC[i]/survival$MedLife[survival.match]
      good.escape$Redds_lower[i]<-good.escape$pred_AUC_lower[i]/survival$MedLife[survival.match]
      good.escape$Redds_upper[i]<-good.escape$pred_AUC_upper[i]/survival$MedLife[survival.match]
      good.escape$Redds_SD[i]<-good.escape$pred_AUC_SD[i]/survival$MedLife[survival.match]
    }
  }
  good.escape$Escape<-good.escape$Redds*mult
  good.escape$Escape_lower<-good.escape$Redds_lower*mult
  good.escape$Escape_upper<-good.escape$Redds_upper*mult
  good.escape$Escape_SD<-good.escape$Redds_SD*mult

  # now make the year totals, which will require dipping back into the simulations
  yearvec<-sort(unique(good.escape$Runyear))
  reachvec<-sort(unique(good.escape[,streamvast$reachname]))

  escapetotals<-data.frame(Runyear=yearvec,Escape=NA,Escape_lower=NA,Escape_median=NA,
                           Escape_upper=NA,Escape_SD=NA)

  escapesims<-cbind(streamvast$sims$aucsims[,c("Runyear",reachname)],
                    matrix(data=NA,ncol=ncol(streamvast$sims$aucsims)-2,
                           nrow=nrow(streamvast$sims$aucsims)))
  escapetotalsims<-cbind(Runyear=yearvec,matrix(data=NA,ncol=ncol(streamvast$sims$auctotalsims)-1,
                                                nrow=length(yearvec)))

  for(y in 1:nrow(escapetotals)){
    picks<-which(escapesims$Runyear==escapetotals$Runyear[y] &
                   escapesims[,reachname]%in%good.reaches)

    yearsim<-streamvast$sims$aucsims[picks,]
    yearsurvival<-subset(survival,Year==escapetotals$Runyear[y])

    yearreachsurv<-yearsurvival$MedLife[match(yearsim[,reachname],yearsurvival$Reach)]
    yearescapesim<-(yearsim[,3:ncol(yearsim)]/yearreachsurv)*mult
    yeartotals<-apply(yearescapesim,MARGIN=2,FUN=sum)

    escapesims[picks,3:ncol(escapesims)]<-yearescapesim
    escapetotalsims[y,2:ncol(escapetotalsims)]<-yeartotals

    escapetotals$Escape[y]<-sum(good.escape$Escape[good.escape$Runyear==escapetotals$Runyear[y]])
    escapetotals$Escape_lower[y]<-stats::quantile(yeartotals,probs=.025)
    escapetotals$Escape_median[y]<-stats::quantile(yeartotals,probs=.5)
    escapetotals$Escape_upper[y]<-stats::quantile(yeartotals,probs=.975)
    escapetotals$Escape_SD[y]<-sqrt(stats::var(yeartotals))
  }

  streamvast$escapedata$escapedata<-good.escape[,c("Runyear",reachname,"Escape",
                                                   "Escape_lower","Escape_upper","Escape_SD")]
  streamvast$escapedata$escapetotals<-escapetotals
  streamvast$sims$escapesims<-escapesims
  streamvast$sims$escapetotalsims<-as.data.frame(escapetotalsims)

  return(streamvast)
}


#' Make a plot of yearly escapement values
#'
#' @param escape A table of escapement values, usually from MakeEscapement
#' @param obs.escape A vector of observed or reported escapement values; must be equal to the # of years in escape
#' @param median should the plot show the predicted escapement based on the mle or the median of the posterior
#' @param title character; a title to be passed to ggtitle
#'
#' @return either a ggplot object or a dataframe depending on the 'plot' option
#' @export
#'
#' @examples
plotEscapement<-function(streamvast,obs.escape,title,ribbons=NA,median=F){

  escapetotals<-streamvast$escapedata$escapetotals
  escapetotals$type<-"Predicted"

  if(median){escapetotals$Escape<-escapetotals$Escape_median}

  # should we plot a confidence ribbon
  if(is.na(ribbons[1])==F){
    # are we setting a number of ribbons
    if(length(ribbons)==1 & ribbons[1]>=1){
      lows<-seq(from=.025,to=.5,length.out=ribbons+1)
      highs<-seq(from=.5,to=.975,length.out=ribbons+1)
    }
    # or specifying a number of ribbons
    if(all(ribbons<1) & all(ribbons>0)){
      lows<-sort(ribbons[ribbons<=.5])
      lows<-unique(c(lows,.5))
      highs<-sort(ribbons[ribbons>=.5])
      highs<-unique(c(.5,highs))
    }
    highs<-rev(highs)
    if(length(lows)!=length(highs)){stop("Ribbon argument must be symmetric")}
    nribs<-length(lows)-1

    ribbon.data<-data.frame(Runyear=rep(escapetotals$Runyear,nribs),
                            ribbon=as.factor(rep(1:nribs,each=nrow(escapetotals))),
                            ymin=NA,ymax=NA,ymin2=NA,ymax2=NA)

    for(i in 1:nrow(ribbon.data)){
      yearsim<-streamvast$sims$escapetotalsims[streamvast$sims$escapetotalsims$Runyear==ribbon.data$Runyear[i],2:ncol(streamvast$sims$escapetotalsims)]
      rib<-as.integer(as.character(ribbon.data$ribbon[i]))
      ribbon.data$ymin[i]<-stats::quantile(unlist(yearsim),probs=lows[rib])
      ribbon.data$ymax[i]<-stats::quantile(unlist(yearsim),probs=lows[rib+1])
      ribbon.data$ymax2[i]<-stats::quantile(unlist(yearsim),probs=highs[rib])
      ribbon.data$ymin2[i]<-stats::quantile(unlist(yearsim),probs=highs[rib+1])
    }
  }

  if(missing(obs.escape)){
    obs.totals<-data.frame(Year=vector(),Escapement=vector(),lower=vector(),
                           upper=vector(),type=vector())
  }else{
    if(any(c("Year","Escapement")%in%names(obs.escape)==F)){
      stop("Argument obs.escape must include columns named 'Year' and 'Escapement'!")
    }
    obs.totals<-data.frame(Runyear=obs.escape$Year,Escape=obs.escape$Escapement,Escape_lower=NA,
                           Escape_median=NA,Escape_upper=NA,Escape_SD=NA,type="Reported")
  }

  escape.plot<-rbind(escapetotals,obs.totals)
  escape.plot$type<-factor(escape.plot$type,levels = c("Reported","Predicted"))

  outplot<-ggplot2::ggplot()+
    ggplot2::geom_point(data=escape.plot,ggplot2::aes(x=Runyear,y=Escape,col=type))+
    ggplot2::geom_line(data=escape.plot,ggplot2::aes(x=Runyear,y=Escape,col=type))+
    ggplot2::geom_line(data=escape.plot,ggplot2::aes(x=Runyear,y=Escape_lower,group=type),col=2,linetype=3)+
    ggplot2::geom_line(data=escape.plot,ggplot2::aes(x=Runyear,y=Escape_upper,group=type),col=2,linetype=3)
  if(is.na(ribbons[1])==F){
    outplot<-outplot+
      ggplot2::geom_ribbon(data=ribbon.data,ggplot2::aes(x=Runyear,ymin=ymin,ymax=ymax,alpha=ribbon),fill=2,show.legend = F)+
      ggplot2::geom_ribbon(data=ribbon.data,ggplot2::aes(x=Runyear,ymin=ymin2,ymax=ymax2,alpha=ribbon),fill=2,show.legend = F)+
      scale_alpha_discrete(range=c(.1,.75))
  }
  outplot<-outplot+
    ggplot2::scale_x_continuous(n.breaks=length(unique(escapetotals$Runyear)))+
    ggplot2::scale_color_manual(values=1:2,drop=F)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.title = ggplot2::element_blank(),legend.position = "inside",
                   legend.position.inside = c(.02,.98),legend.justification = c(0,1))+
    ggplot2::ylab("Escapement")+ggplot2::ylim(c(0,NA))
  if(missing(title)==F){outplot<-outplot+ggplot2::ggtitle(title)}

  return(outplot)
}


#' Make a dataframe with outputs from the VAST model for reporting
#'
#' @param streamvast a streamvast object with completed predictions
#' @param watershed  character; a name for the watershed
#'
#' @return
#' @export
#'
#' @examples
EscapeReport<-function(streamvast,watershed){

  report.data<-data.frame(Time=streamvast$preds$Time,Watershed=rep(watershed,nrow(streamvast$preds)),
                          Stream=NA,Year=streamvast$preds$Year,
                          Reach=streamvast$preds$Reach,Month=streamvast$preds$Month,Day=streamvast$preds$Day,
                          StartLon=NA,StartLat=NA,EndLon=NA,EndLat=NA,
                          ObsDensity=NA,ObsRedds=NA,ObsEffort=NA,PredDensity=round(streamvast$preds$pred_Density,3),
                          Pred_LowerCI=round(streamvast$preds$pred_Density_lower,3),
                          Pred_UpperCI=round(streamvast$preds$pred_Density_upper,3))

  reaches.gps<-sf::st_transform(streamvast$reachdata,"wgs84")

  for(i in 1:nrow(report.data)){
    reach.match<-which(reaches.gps$Reach==report.data$Reach[i])

    reach.coords<-sf::st_coordinates(reaches.gps[reach.match,])
    report.data$Stream[i]<-reaches.gps$STRM_NA[reach.match]
    report.data$StartLon[i]<-round(reach.coords[1,1],2)
    report.data$StartLat[i]<-round(reach.coords[1,2],2)
    report.data$EndLon[i]<-round(reach.coords[nrow(reach.coords),2],1)
    report.data$EndLat[i]<-round(reach.coords[nrow(reach.coords),2],2)

    vast.match<-which(streamvast$vastdata$Time==report.data$Time[i] &
                        streamvast$vastdata$Reach==report.data$Reach[i])
    is.dummy<-any(streamvast$vastdata$dummy[vast.match])
    if(is.dummy){
      report.data$ObsDensity[i]<-NA
      report.data$ObsRedds[i]<-NA
      report.data$ObsEffort[i]<-0
    }else{
      report.data$ObsDensity[i]<-round(mean(streamvast$vastdata$Density[vast.match]),3)
      report.data$ObsRedds[i]<-sum(streamvast$vastdata$Redds[vast.match])
      report.data$ObsEffort[i]<-round(sum(streamvast$vastdata$Effort[vast.match]),3)
    }
  }
  return(report.data[,-1])
}


#' Make a table with relevant outputs from the VAST and survival models
#'
#' @param streamvast      A streamvast
#' @param survival        A table of survival values; usually from SurvivalTable
#' @param escape          A table of escapement values, usually from MakeEscapement
#' @param watershed       character; a name for the watershed
#'
#' @return
#' @export
#'
#' @examples
EscapeReport2<-function(streamvast,survival,escape,watershed){

  report.data<-data.frame(Watershed=rep(watershed,nrow(streamvast$aucdata$aucdata)),
                          Stream=NA,Year=streamvast$aucdata$aucdata$Year,
                          Reach=streamvast$aucdata$aucdata$Reach,StartLon=NA,
                          StartLat=NA,EndLon=NA,EndLat=NA,ObsDensity=NA,
                          ObsRedds=NA,ObsEffort=NA,PredDensity=NA,Survival=NA,
                          PredRedds=NA,Escapement=NA,Escape_UpperCI=NA,Escape_LowerCI=NA)

  reaches.gps<-sf::st_transform(streamvast$reachdata,"wgs84")

  for(i in 1:nrow(report.data)){
    reach.match<-which(reaches.gps$Reach==report.data$Reach[i])

    reach.coords<-sf::st_coordinates(reaches.gps[reach.match,])
    report.data$Stream[i]<-reaches.gps$STRM_NA[reach.match]
    report.data$StartLon[i]<-round(reach.coords[1,1],2)
    report.data$StartLat[i]<-round(reach.coords[1,2],2)
    report.data$EndLon[i]<-round(reach.coords[nrow(reach.coords),2],1)
    report.data$EndLat[i]<-round(reach.coords[nrow(reach.coords),2],2)

    vast.match<-which(streamvast$vastdata$Year==report.data$Year[i] &
                        streamvast$vastdata$Reach==report.data$Reach[i] &
                        streamvast$vastdata$dummy==F)
    if(length(vast.match)==0){
      report.data$ObsDensity[i]<-NA
      report.data$ObsRedds[i]<-NA
      report.data$ObsEffort[i]<-0
    }else{
      report.data$ObsDensity[i]<-round(mean(streamvast$vastdata$Density[vast.match]),3)
      report.data$ObsRedds[i]<-sum(streamvast$vastdata$Redds[vast.match])
      report.data$ObsEffort[i]<-round(sum(streamvast$vastdata$Effort[vast.match]),3)
    }

    pred.match<-which(streamvast$preds$Year==report.data$Year[i] &
                        streamvast$preds$Reach==report.data$Reach[i])
    report.data$PredDensity[i]<-round(mean(vast$preds$pred_Density[pred.match]),2)

    if(is.data.frame(survival)){
      surv.match<-which(survival$Year==report.data$Year[i] &
                          survival$Reach==report.data$Reach[i])
      report.data$Survival[i]<-round(survival$MedLife[surv.match],2)
    }else{
      report.data$Survival[i]<-survival
    }

    escape.match<-which(escape$Year==report.data$Year[i] &
                          escape$Reach==report.data$Reach[i])
    report.data$PredRedds[i]<-round(escape$Redds[escape.match],2)
    report.data$Escapement[i]<-round(escape$Escape[escape.match],2)
    report.data$Escape_LowerCI[i]<-round(escape$Escape_lower[escape.match],2)
    report.data$Escape_UpperCI[i]<-round(escape$Escape_upper[escape.match],2)
  }
  return(report.data)
}
