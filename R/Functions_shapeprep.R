#' FindRoots
#'
#' This function locates all vertices connected to the shape specified in mainstem
#'
#' @param shape An sf object with linestring geometry and Identifier data
#' @param mainstem a name or unique id for the mainstem that basins branch from
#'
#' @return An sf object with points for each potential root
#' @export
#'
#' @examples
#'
#'
FindRoots<-function(shape,mainstem=NA){

  # First make sure the mainstem is identified properly
  name.col<-which(apply(shape,MARGIN = 2,FUN = function(x){return(any(x==mainstem))}))
  identifier<-mainstem

  if(length(name.col)!=1){
    stop("Invalid identifier for mainstem or id; Either no match, or more than one match was found")
  }

  #we could check for an id key or something, but it's easier to just make one ourselves
  shape$IDKEY<-1:nrow(shape)

  # then separate the mainstem from the rest
  shape.main<-shape[which(as.data.frame(shape)[,name.col]==identifier),]
  shape.main<-sf::st_union(shape.main)
  shape.sub<-shape[-which(as.data.frame(shape)[,name.col]==identifier),]

  # record any points of contact
  Root.points<-sf::st_intersection(y=shape.main,x=shape.sub)

  # now we need to worry about side channels, sloughs, islands, and other minor diversions
  # if it intersects the mainstem multiple times, then its probably a side channel
  side.ids<-Root.points$IDKEY[which(sf::st_geometry_type(Root.points)!="POINT")]
  side.channels<-shape.sub[shape.sub$IDKEY%in%side.ids,]

  # Let's also set up the non-side roots so that they are done and out of the way
  main.points<-Root.points[sf::st_geometry_type(Root.points)=="POINT",]
  main.points$parent<-mainstem
  main.points<-data.frame(as.data.frame(sf::st_coordinates(main.points)),
                          as.data.frame(main.points)[,name.col],
                          as.data.frame(main.points)[,"parent"])
  names(main.points)[3:4]<-c(names(shape)[name.col],"parent")

  # now, we check each side channel, splitting the side channel from the rest of the tree
  # in theory, this could get endlessly recursive, with side channels to side channels, etc
  # but we will stop after one layer and spit out a warning instead

  # Loop through all the side channels
  side.points.out<-main.points[0,]
  for(i in 1:nrow(side.channels)){
    # now, we check each side channel, splitting the side channel from the rest of the tree
    shape.side<-side.channels[i,]
    shape.side.sub<-shape.sub[which(shape.sub$IDKEY!=shape.side$IDKEY[i]),]

    # and find the intersections
    side.points<-sf::st_intersection(y=shape.side,x=shape.side.sub)

    if(nrow(side.points)>0){
      side.points$parent<-as.data.frame(shape.side)[,name.col]

      # in theory, this could get endlessly recursive, with side channels to side channels, etc
      # but we will stop after one layer and spit out a warning instead
      side.points<-side.points[sf::st_geometry_type(side.points)=="POINT",]

      if(any(sf::st_geometry_type(side.points)!="POINT")){
        warning(paste0("Potential additional roots along side channel ",as.data.frame(shape.side)[1,name.col],
                       " ; feature # ",as.data.frame(shape.side)[1,"IDKEY"]))
      }

      # clean it up and get it ready for output
      side.points<-data.frame(as.data.frame(sf::st_coordinates(side.points)),
                              as.data.frame(side.points)[,name.col],
                              as.data.frame(side.points)[,"parent"])
      names(side.points)[3:4]<-c(names(shape)[name.col],"parent")
      side.points$parent[is.na(side.points$parent)]<-paste0(mainstem," side channel")

      side.points.out<-rbind(side.points.out,side.points)
    }
  }
  Root.coords<-rbind(main.points,side.points.out)

  out.shape<-sf::st_as_sf(Root.coords,coords=c("X","Y"),crs=sf::st_crs(shape))

  plot(sf::st_geometry(shape.sub),col="gray50")
  plot(sf::st_sf(shape.main),col="blue",add=T,lwd=3)
  plot(sf::st_geometry(out.shape),col=2,cex=2,add=T,pch=16)

  return(out.shape)
}


#' Make a sf network
#'
#' This function converts a sf shape with linestring or multilinestring geometries into
#' a sfnetwork object. Mostly, it is used internally.
#'
#' @param shape An sf object with LINESTRING or MULTILINESTRING geometry
#' @param attach.data Logical, whether data fields should be preserved
#' @param simple Logical, whether to use the fast and simple method, or the slower more thorough method
#'
#' @return A sfnetwork object
#' @export
#'
#' @examples
Makesfnetwork<-function(shape,attach.data=F,simple=T){

  shape.net<-sfnetworks::as_sfnetwork(sf::st_as_sf(shape),)
  shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape.net,"nodes"))
  shape.edges<-sf::st_as_sf(sfnetworks::activate(shape.net,"edges"))

  dosplit<-ifelse(simple,F,T)

  if(simple){
    connect.check<-sfnetworks::st_network_cost(shape.net,1,direction="all")[1,]
    dosplit<-any(is.infinite(connect.check))
    if(dosplit){print("Converting shape to sfnetwork using the simple method")}
  }else{
    print("Converting shape to sfnetwork using the thorough method")
    # now an advanced step that goes piece by piece to prevent R from running over memory limits
    int.mat<-sf::st_intersects(shape.edges)
    pb<-utils::txtProgressBar(min = 0, max = nrow(shape.edges), initial = 0, style=3)
    for(i in 1:length(int.mat)){

      int.shapes<-shape.edges[int.mat[[i]],]
      intersect.shapes<-sf::st_as_sf(sf::st_as_sfc(unique(sf::st_geometry(sf::st_intersection(int.shapes)))),crs=sf::st_crs(shape))

      int.points0<-intersect.shapes[sf::st_geometry_type(intersect.shapes)=="POINT",]
      int.lines0<-intersect.shapes[sf::st_geometry_type(intersect.shapes)=="LINESTRING",]
      int.mpoints<-sf::st_cast(intersect.shapes[sf::st_geometry_type(intersect.shapes)=="MULTIPOINT",],"POINT")
      int.mlines<-sf::st_cast(intersect.shapes[sf::st_geometry_type(intersect.shapes)=="MULTILINESTRING",],"LINESTRING")

      int.points<-c(sf::st_geometry(int.points0),sf::st_geometry(int.mpoints))
      int.lines<-c(sf::st_geometry(int.lines0),sf::st_geometry(int.mlines))

      # This bit deals with lines that may be slightly overlapping
      l.mat<-sf::st_equals(int.lines,shape.edges)
      lines.to.blend<-unlist(lapply(l.mat,FUN=function(x){return(length(x)==0)}))
      if(any(lines.to.blend)){
        int.line.pts<-sf::st_cast(int.lines[which(lines.to.blend)],"POINT")
        int.points<-c(int.points,int.line.pts)
      }

      int.points2<-sf::st_as_sf(sf::st_as_sfc(unique(int.points),crs=sf::st_crs(shape)))
      p.mat<-sf::st_equals(int.points2,shape.nodes)
      points.to.blend<-unlist(lapply(p.mat,FUN=function(x){return(length(x)==0)}))
      if(any(points.to.blend)){
        shape.nodes<-rbind(shape.nodes,int.points2[points.to.blend,])
      }
      utils::setTxtProgressBar(pb,i)
    }
    close(pb)
  }

  if(dosplit){
    # now to cut at those points
    shape.edges.combined<-sf::st_union(shape.edges)
    new.edges0<-lwgeom::st_split(x=shape.edges.combined,y=shape.nodes)
    new.edges<-sf::st_sf(sf::st_collection_extract(new.edges0,"LINESTRING"))
    sf::st_geometry(new.edges)<-attr(sf::st_as_sf(shape),"sf_column")

    if(attach.data){
      new.edges2<-AttachData(shape = new.edges,refshape = shape)
    }else{
      new.edges2<-new.edges
    }
  }else{
    new.edges2<-shape.edges
  }

  return(sfnetworks::as_sfnetwork(new.edges2))
}

# runs a quick format check that a given point is a proper sf point
#' Point Setup and checking
#'
#' Checks that a point is an sf object with the appropriate crs, or converts inputs into an sf
#' object with POINT geometry; mostly used internally to check user inputs
#'
#' @param point A sf, sfc, dataframe, or vector containing geometry or coordinates representing a point
#' @param crs A crs to apply to the point, if it is not an sf or sfc object
#'
#' @return A sf object with POINT geometry
#' @export
#'
#' @examples
PointSetup<-function(point,crs){
  # Setup the point
  if(class(point)[1]=="sf"){point.sf<-point}
  if(class(point)[1]=="sfc_POINT"){point.sf<-sf::st_as_sf(point)}
  if(class(point)[1]!="sf" & is.data.frame(point)){point.sf<-sf::st_as_sf(x=point,coords=1:2,crs=crs)}
  if(class(point)[1]%in%c("integer","numeric")){point.sf<-sf::st_as_sf(x=data.frame(X=point[1],Y=point[2]),
                                                                       coords=1:2,crs=crs)}
  if(class(point)[1]%in%c("sf","sfc_POINT","integer","numeric")==F & is.data.frame(point)==F){point.sf<-NA}
  sf::st_geometry(point.sf)<-"geometry"
  return(point.sf)
}

# runs a quick function that returns the root of a network
#' Locate the root of a network
#'
#' Locates the root of a network, assuming it has one. If the network has unconnected nodes
#' or does not have a valid root, returns an error
#'
#' @param network A sf network object with a valid root
#'
#' @return A sf object with POINT geometry representing the root of the network
#' @export
#'
#' @examples
LocateRoot<-function(network){

  if(inherits(network,"sf")){
    network<-suppressWarnings(sfnetworks::as_sfnetwork(network))
  }
  connection.fail<-any(is.infinite(sfnetworks::st_network_cost(network,from = 1,direction = "all")[1,]))
  if(connection.fail){
    warning("Not all nodes connected to the network! Root cannot be located!")
    return(NA)
  }

  network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
  network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))

  cost.mat<-sfnetworks::st_network_cost(network)
  root.index<-which(apply(cost.mat,MARGIN = 1,FUN=function(x){return(any(is.infinite(x))==F)}))
  if(length(root.index)==0){
    warning("Network does not have a root!")
    return(NA)
  }

  root.node<-network.nodes[root.index,]
  root.node$Index<-root.index
  return(root.node)
}

#' Adds or changes the root of a network
#'
#'  This function sets up the root and ensures that all nodes are connected to it, and
#'  preserves any old data. Any nodes not connected to the root are removed. By default, the
#'  function will root the network at the closest existing node to the specified root. Set tolerance = 0
#'  to instead cause the function to create a new node
#'
#' @param network A sfnetwork
#' @param root A sf point, or coordinates for the new root
#' @param root.crs a crs or string for the root, if it is different from the network
#' @param tolerance A distance that determines whether to use existing nodes or blend in a new one
#'
#' @return A sfnetwork with the new root
#' @export
#'
#' @examples
RootNetwork<-function(network,root,root.crs,tolerance=Inf){

  if(missing(root.crs)){
    root.sf<-PointSetup(point = root,crs = sf::st_crs(network))
  }else{
    root.sf0<-PointSetup(point = root,crs = root.crs)
    root.sf<-sf::st_transform(root.sf0,crs=sf::st_crs(network))
  }

  network.nodes0<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))
  network.edges0<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
  network0<-network
  root.node.index<-sf::st_nearest_feature(root.sf,network.nodes0)

  # If the specified root is too far from the nearest node, blend it into the network
  if(as.numeric(sf::st_distance(root.sf,network.nodes0[root.node.index,]))>tolerance){
    network0<-suppressWarnings(sfnetworks::st_network_blend(network,root.sf))
    network.nodes0<-sf::st_as_sf(sfnetworks::activate(network0,"nodes"))
    network.edges0<-sf::st_as_sf(sfnetworks::activate(network0,"edges"))

    root.node.index<-sf::st_nearest_feature(root.sf,network.nodes0)
  }
  # First, drop anything not connected to the root
  node.costs<-sfnetworks::st_network_cost(network0,from=root.node.index,to=network.nodes0,direction = "all")[1,]
  bad.nodes<-which(is.infinite(node.costs))
  if(length(bad.nodes)>0){
    network2<-suppressWarnings(sfnetworks::as_sfnetwork(network.edges0[network.edges0$to%in%bad.nodes==F &
                                                                         network.edges0$from%in%bad.nodes==F,]))
  }else{
    network2<-network0
  }

  network.nodes<-sf::st_as_sf(sfnetworks::activate(network2,"nodes"))
  network.edges<-sf::st_as_sf(sfnetworks::activate(network2,"edges"))
  root.node.index2<-sf::st_nearest_feature(root.sf,network.nodes)

  # first, check that all stream segments are going the right direction
  root.distances<-sfnetworks::st_network_cost(network2,from=root.node.index2,to=network.nodes,direction = "all")[1,]
  reversed<-which(root.distances[network.edges$from]>root.distances[network.edges$to])

  new.network.edges<-network.edges
  if(length(reversed>0)){
    for(i in 1:length(reversed)){
      index<-reversed[i]
      new.network.edges$from[index]<-network.edges$to[index]
      new.network.edges$to[index]<-network.edges$from[index]

      old.geo<-sf::st_geometry(network.edges)[index]
      old.coords<-sf::st_coordinates(old.geo)[,1:2]
      new.coords<-old.coords[nrow(old.coords):1,]
      new.geo<-sf::st_cast(sf::st_combine(sf::st_as_sf(as.data.frame(new.coords),coords=1:2,
                                                       crs=sf::st_crs(network))),"LINESTRING")
      sf::st_geometry(new.network.edges)[index]<-new.geo
    }
  }
  return(suppressWarnings(sfnetworks::as_sfnetwork(new.network.edges)))
}

#' Locate a feature on a shape or network
#'
#' Allows the user to select a number of shapes from a sf object or network by clicking
#' in the plot window. A shape will be converted to a sfnetwork first
#'
#' @param shape A sf object with LINESTRING geometry or a sfnetwork
#' @param type Either "POINT", "VERTEX" or "SEGMENT" which indicates the type of shape to select
#'
#' @return A sf object containing the POINTS or LINESTRINGS selected by the user
#' @export
#'
#' @examples
LocateFeature<-function(shape,type){

  if(toupper(type)%in%c("POINT","P","PT")){type<-"POINT"}
  if(toupper(type)%in%c("V","VER","VRTX","VERTEX")){type<-"VERTEX"}
  if(toupper(type)%in%c("S","SEG","SGMT","SEGMENT")){type<-"SEGMENT"}
  if(type%in%c("POINT","VERTEX","SEGMENT")==F){stop("Invalid 'type' argument!")}

  # Setup some plot parameters
  bounds<-sf::st_bbox(shape)
  xspan<-(bounds$xmax-bounds$xmin)
  yspan<-(bounds$ymax-bounds$ymin)
  midpoint<-c(X=mean(c(bounds$xmin,bounds$xmax)),Y=mean(c(bounds$ymin,bounds$ymax)))

  # set up a network
  if(class(shape)[1]=="sfnetwork"){
    shape.edges<-sf::st_as_sf(sfnetworks::activate(shape,"edges"))
    shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape,"nodes"))
    shape.net<-shape
  }else{
    if(type=="VERTEX"){warning("Vertex Indexing may be unstable if shape is not a sfnetwork")}
    shape.net<-suppressWarnings(sfnetworks::as_sfnetwork(shape))
    shape.edges<-sf::st_as_sf(sfnetworks::activate(shape.net,"edges"))
    shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape.net,"nodes"))
  }
  # output object
  if(type%in%c("POINT","VERTEX")){out.sf<-shape.nodes[0,]}
  if(type=="SEGMENT"){out.sf<-shape.edges[0,]}

  still.going<-T
  while(still.going){
    # start by plotting the network and show any updates
    plot(sf::st_geometry(shape.edges),xlim=c(midpoint[1]-xspan/2,midpoint[1]+xspan/2),
         ylim=c(midpoint[2]-yspan/2,midpoint[2]+yspan/2),pch=16)
    plot(sf::st_geometry(shape.nodes),pch=16,add=T,col=1)
    plot(sf::st_geometry(out.sf),add=T,col=4,pch=16,cex=1.25,lwd=1.25)

    # user chooses an action
    response<-readline(prompt="Choose an action : CENTER/zoomIN/zoomOUT/SELECT/FINISH")

    # adjust for spelling or abbreviations by the user
    if(toupper(response)%in%c("C","CEN","CENT","CNTR","CENTER")){resp<-"CENTER"}
    if(toupper(response)%in%c("S","SEL","SE","SLCT","SELECT")){resp<-"SELECT"}
    if(toupper(response)%in%c("F","FIN","FI","FINISH")){resp<-"FINISH"}
    if(toupper(response)%in%c("I","IN","ZOOMIN")){resp<-"IN"}
    if(toupper(response)%in%c("O","OUT","ZOOMOUT")){resp<-"OUT"}

    # should consider adding code to handle typos without an error


    # if finished, flag to not repeat the loop
    if(resp=="FINISH"){still.going<-F}

    # if center or zoom, adjust plot parameters
    if(resp%in%c("IN","OUT")){
      scale<-ifelse(resp=="IN",2,.5)
      xspan<-xspan/scale
      yspan<-yspan/scale
    }
    if(resp%in%c("C","CEN","CENTER")){
      print("Please click on the image to center.")
      midpoint<-unlist(graphics::locator(n=1))
    }

    # if selecting a feature, prompt the user to pick
    if(resp=="SELECT"){
      print("Please click on the image to select a feature")
      user.pick<-unlist(graphics::locator(n=1))
      user.sf<-sf::st_as_sf(data.frame(x=user.pick[1],y=user.pick[2]),coords=1:2,crs=sf::st_crs(shape))

      # because of CompSci shenanigans, this point is often not quite on the line, so fix that
      if(type=="POINT"){
        shape.line<-shape.edges[sf::st_nearest_feature(user.sf,shape.edges),]
        pick0<-sf::st_as_sf(sf::st_cast(sf::st_nearest_points(user.sf,shape.line),"POINT"))[2,]
        temp.net<-suppressWarnings(sfnetworks::st_network_blend(shape.net,pick0))
        temp.nodes<-sf::st_as_sf(sfnetworks::activate(temp.net,"nodes"))
        temp.edges<-sf::st_as_sf(sfnetworks::activate(temp.net,"edges"))
        pick.sf<-temp.nodes[sf::st_nearest_feature(pick0,temp.nodes),]
        plot(sf::st_geometry(pick.sf),add=T,pch=21,bg=6)
      }
      if(type=="VERTEX"){
        pick.index<-sf::st_nearest_feature(user.sf,shape.nodes)
        pick.sf<-shape.nodes[pick.index,]
        pick.sf$Index<-pick.index
        plot(sf::st_geometry(pick.sf),add=T,pch=21,bg=6)
      }
      if(type=="SEGMENT"){
        pick.index<-sf::st_nearest_feature(user.sf,shape.edges)
        pick.sf<-shape.edges[pick.index,]
        pick.sf$Index<-pick.index
        plot(sf::st_geometry(pick.sf),add=T,col=6,lwd=2)
      }
    }
    # Have the user confirm the correct feature was selected
    if(resp=="SELECT"){
      good.pick<-toupper(readline(prompt = "Is this correct (Y/N)?   "))
      if(toupper(good.pick)%in%c("Y","YES","T",T)){
        out.sf<-rbind(out.sf,pick.sf)
      }
    }
  }
  sf::st_geometry(out.sf)<-"geometry"
  return(out.sf)
}


#' Edit a shape or network in real time
#'
#' Allows the user to add points or remove features by selecting them from the plot window
#' and command line prompts
#'
#' @param shape A sf or sfnetworks object to be editted
#' @param root A sf POINT for the root (might be obsolete)
#'
#' @return A sf or sfnetworks object with the user specificed alterations
#' @export
#'
#' @examples
EditFeatures<-function(shape,root){

  if(class(shape)[1]!="sfnetwork"){
    shape.net<-Makesfnetwork(shape)
  }else{
    shape.net<-shape
  }
  shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape.net,"nodes"))
  shape.edges<-sf::st_as_sf(sfnetworks::activate(shape.net,"edges"))

  still.going<-T

  # Setup some plot parameters
  bounds<-sf::st_bbox(shape.edges)

  xspan<-(bounds$xmax-bounds$xmin)
  yspan<-(bounds$ymax-bounds$ymin)

  cost.mat<-sfnetworks::st_network_cost(shape.net)
  root.index<-which(apply(cost.mat,MARGIN = 1,FUN=function(x){return(any(is.infinite(x))==F)}))
  has.root<-length(root.index)>0
  if(has.root){root.sf<-shape.nodes[root.index,]}

  #make a copy of the shapes, so that we have original and updated version
  temp.net<-shape.net
  temp.edges<-shape.edges
  temp.nodes<-shape.nodes
  new.points<-shape.nodes[0,]
  midpoint<-c(X=mean(c(bounds$xmin,bounds$xmax)),Y=mean(c(bounds$ymin,bounds$ymax)))

  # user input loop
  while(still.going){

    # start by plotting the network and show any updates
    plot(sf::st_geometry(shape.edges),xlim=c(midpoint[1]-xspan/2,midpoint[1]+xspan/2),
         ylim=c(midpoint[2]-yspan/2,midpoint[2]+yspan/2))
    plot(sf::st_geometry(shape.edges),add=T,col=2)
    plot(sf::st_geometry(shape.nodes),pch=21,add=T,col=1,bg=2)
    plot(sf::st_geometry(temp.edges),add=T,col=1)
    plot(sf::st_geometry(temp.nodes),add=T,col=1,pch=16)
    if(nrow(new.points)>0){plot(sf::st_geometry(new.points),add=T,pch=21,col=1,bg=4)}

    # user chooses an action
    response<-readline(prompt="Choose an action : CENTER/zoomIN/zoomOUT/ADD/REMOVE/FINISH")

    # adjust for spelling or abbreviations by the user
    if(toupper(response)%in%c("C","CEN","CENT","CNTR","CENTER")){resp<-"CENTER"}
    if(toupper(response)%in%c("A","AD","ADD")){resp<-"ADD"}
    if(toupper(response)%in%c("R","RE","RMOVE","RM","REMOVE")){resp<-"REMOVE"}
    if(toupper(response)%in%c("S","SEL","SE","SLCT","SELECT")){resp<-"SELECT"}
    if(toupper(response)%in%c("F","FIN","FI","FINISH")){resp<-"FINISH"}
    if(toupper(response)%in%c("I","IN","ZOOMIN")){resp<-"IN"}
    if(toupper(response)%in%c("O","OUT","ZOOMOUT")){resp<-"OUT"}

    # if finished, flag to not repeat the loop
    if(resp=="FINISH"){still.going<-F}

    # if center or zoom, adjust plot parameters
    if(resp%in%c("IN","OUT")){
      scale<-ifelse(resp=="IN",2,.5)
      xspan<-xspan/scale
      yspan<-yspan/scale
    }
    if(resp%in%c("C","CEN","CENTER")){
      print("Please click on the image to center.")
      midpoint<-unlist(graphics::locator(n=1))
    }

    # Can add a point to the network
    if(resp=="ADD"){
      print("Please click on the image to add a point")
      user.pick<-unlist(graphics::locator(n=1))
      user.sf<-sf::st_as_sf(data.frame(x=user.pick[1],y=user.pick[2]),coords=1:2,crs=sf::st_crs(shape))

      shape.line<-temp.edges[sf::st_nearest_feature(user.sf,temp.edges),]
      pick0<-sf::st_as_sf(sf::st_cast(sf::st_nearest_points(user.sf,shape.line),"POINT"))[2,]

      # because of CompSci shenanigans, this point is often not quite on the line
      # but sf_networks is good at fixing the problem
      temp.net0<-suppressWarnings(sfnetworks::st_network_blend(temp.net,pick0))
      temp.nodes0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"nodes"))
      temp.edges0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"edges"))
      pick.sf<-temp.nodes0[sf::st_nearest_feature(pick0,temp.nodes0),]
      plot(sf::st_geometry(pick.sf),add=T,pch=21,bg=3)
    }

    # if removing features, the user must choose vertex or segment
    if(resp=="REMOVE"){
      type<-readline(prompt = "What type of feature would you like to select?  (VERTEX/SEGMENT)")
      if(toupper(type)%in%c("V","VER","VRTX","VERTEX","POINT","P","PT")){type<-"VERTEX"}
      if(toupper(type)%in%c("S","SEG","SGMT","SEGMENT")){type<-"SEGMENT"}

      print("Please click on the image to select a feature")
      user.pick<-unlist(graphics::locator(n=1))
      user.sf<-sf::st_as_sf(data.frame(x=user.pick[1],y=user.pick[2]),coords=1:2,crs=sf::st_crs(shape))

      # When removing a vertex, you can't remove the root
      if(type=="VERTEX"){
        pick.index<-sf::st_nearest_feature(user.sf,temp.nodes)

        if(has.root & pick.index==root.index){
          warning("Can not remove the root! Selection canceled!")
          resp<-"none"
        }

        connect.edges<-temp.edges[temp.edges$from==pick.index | temp.edges$to==pick.index,]

        #if the vertex is an endpoint, its safe to delete the segment as well
        if(nrow(connect.edges)==1){
          bad.edges<-which(temp.edges$to==pick.index)
          plot(sf::st_geometry(temp.edges[bad.edges,]),add=T,lwd=2,col=6)
          temp.net0<-sfnetworks::as_sfnetwork(temp.edges[-bad.edges,])
        }

        # if the vertex has exactly two edges, then it is internal, so delete the vertex and combine the edges
        if(nrow(connect.edges)==2){
          old.edges<-temp.edges[temp.edges$from!=pick.index & temp.edges$to!=pick.index,]
          new.edges<-sf::st_as_sf(c(sf::st_geometry(old.edges),sf::st_line_merge(sf::st_combine(connect.edges))))
          temp.net0<-sfnetworks::as_sfnetwork(new.edges)
        }

        # if the vertex has 3 or more edges, then it's important to have a root
        if(nrow(connect.edges)>=3){
          if(has.root==F){
            warning("Can not resolve network without root! Remove other features first!")
            resp<-"none"
            temp.net0<-temp.net
          }else{

            bad.edges<-connect.edges
            bad.nodes<-temp.nodes[unique(c(connect.edges$from,connect.edges$to)),]

            temp.net0<-sfnetworks::as_sfnetwork(temp.edges[temp.edges$from!=pick.index & temp.edges$to!=pick.index,])
            temp.edges0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"edges"))
            temp.nodes0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"nodes"))

            bad.nodes.index<-which(is.infinite(sfnetworks::st_network_cost(x=temp.net0,from = root.sf)[1,]))
            bad.edges.index<-which(temp.edges0$from%in%bad.nodes.index | temp.edges0$to%in%bad.nodes.index)

            bad.edges<-rbind(bad.edges,temp.edges0[bad.edges.index,])
            bad.nodes<-rbind(bad.nodes,temp.nodes0[bad.nodes.index,])

            plot(sf::st_geometry(bad.edges),add=T,lwd=2,col=6)
            plot(sf::st_geometry(bad.nodes),add=T,col=6,pch=16)
            temp.net0<-sfnetworks::as_sfnetwork(temp.edges0[-bad.edges.index,])
          }
        }
        temp.nodes0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"nodes"))
        temp.edges0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"edges"))

        plot(sf::st_geometry(temp.nodes[pick.index,]),add=T,pch=21,col=1,bg=6,cex=1.5)
      }

      # if removing segments, can remove just a segment or prune the tree
      if(type=="SEGMENT"){
        pick.index<-sf::st_nearest_feature(user.sf,temp.edges)

        temp.net0<-sfnetworks::as_sfnetwork(temp.edges[-pick.index,-c(1:2)])
        if(has.root){

          temp.nodes0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"nodes"))
          temp.edges0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"edges"))

          bad.nodes.index<-which(apply(sf::st_intersects(temp.nodes,temp.nodes0,sparse=F),MARGIN = 1,FUN = sum)==0)
          bad.nodes.index2<-which(is.infinite(sfnetworks::st_network_cost(temp.net0,root.sf)[1,]))
          bad.edges.index<-which(temp.edges0$from%in%bad.nodes.index2 | temp.edges0$to%in%bad.nodes.index2)
          if(length(bad.nodes.index)>0){plot(sf::st_geometry(temp.nodes[bad.nodes.index,]),add=T,col=6,pch=16)}
          if(length(bad.edges.index)>0){
            plot(sf::st_geometry(temp.edges0[bad.edges.index,]),add=T,lwd=2,col=6)
            temp.net0<-sfnetworks::as_sfnetwork(temp.edges0[-bad.edges.index,-c(1:2)])
            temp.edges0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"edges"))
            temp.nodes0<-sf::st_as_sf(sfnetworks::activate(temp.net0,"nodes"))
          }
        }
        plot(sf::st_geometry(temp.edges[pick.index,]),add=T,col=6,lwd=2)
      }
    }

    # Have the user confirm the correct feature was selected
    if(resp%in%c("REMOVE","ADD")){
      good.pick<-toupper(readline(prompt = "Is this correct (Y/N)?   "))
      if(toupper(good.pick)%in%c("Y","YES","T",T)){
        temp.net<-temp.net0
        temp.edges<-sf::st_as_sf(sfnetworks::activate(temp.net,"edges"))
        temp.nodes<-sf::st_as_sf(sfnetworks::activate(temp.net,"nodes"))

        if(has.root){
          cost.mat<-sfnetworks::st_network_cost(temp.net)
          root.index<-which(apply(cost.mat,MARGIN = 1,FUN=function(x){return(any(is.infinite(x))==F)}))
        }
        if(resp=="ADD"){new.points<-rbind(new.points,pick.sf)}
      }
    }
  }
  return(temp.net)
}


#' Add nodes or edges to a sfnetwork
#'
#' This function adds additional nodes or edge segments to a network. Before adding anything, it checks that the new features
#' are distinct from pre-existing features, and that they are relatively close to the existing network.
#' Note, this function has not been tested extensively.
#'
#' @param network a sfnetwork object
#' @param nodes a data frame with lon and lat columns, or sf object with POINT geometry
#' @param edges a sf object with LINESTRING geometry that includes edges to be added
#' @param crs a crs object, if nodes or edges has a different crs from network
#' @param tolerance a maximum distance for snapping a shape to the network, shapes outside this tolerance are removed
#' @param tolerance2 a maximum distance for snapping nodes together, new nodes within this tolerance are not distinct from existing nodes
#'
#' @return A sfnetwork with the indicated nodes added
#' @export
#'
#' @examples
AddFeatures<-function(network,nodes,edges,crs,tolerance=200,tolerance2=10){

  old.net<-network
  if(class(network)[1]!="sfnetwork"){
    print("Converting shape to sfnetwork")
    network<-Makesfnetwork(network)
  }
  network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
  network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))

  if(missing(crs)){crs<-sf::st_crs(network)}

  # do edges first
  if(missing(edges)==F){
    if(class(edges)[1]!="sf" | sf::st_geometry_type(edges)!="LINESTRING"){stop("Edges must be a sf object with LINESTRING geometry")}
    edges2<-sf::st_transform(edges,crs = crs)
    coords<-as.data.frame(sf::st_coordinates(edges2))
    nbad<-0
    for(i in 1:nrow(edges)){
      line.dat<-subset(coords,L1==i)

      start.node<-sf::st_as_sf(line.dat[1,],coords=1:2,crs=crs)
      end.node<-sf::st_as_sf(line.dat[nrow(line.dat),],coords=1:2,crs=crs)

      start.good<-any(as.numeric(sf::st_distance(start.node,network.edges))<=tolerance)
      end.good<-any(as.numeric(sf::st_distance(end.node,network.edges))<=tolerance)

      # now check for redundant nodes, and then update
      if(any(c(start.good,end.good))){
        if(start.good){
          add.start<-any(as.numeric(sf::st_distance(start.node,network.nodes))<=tolerance2)==F
          if(add.start){
            network<-sfnetworks::st_network_blend(x = network,y = sf::st_geometry(start.node))
          }
        }
        if(end.good){
          add.end<-any(as.numeric(sf::st_distance(end.node,network.nodes))<=tolerance2)==F
          if(add.end){network<-sfnetworks::st_network_blend(x = network,y = sf::st_geometry(end.node))}
        }
        # update the network
        network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
        network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))

        # now the new nodes are in there, so let's add the linestring; no fancy transformation
        line.dat2<-line.dat[,1:2]
        if(start.good){
          start.node2<-network.nodes[sf::st_nearest_feature(x = start.node,y=network.nodes),]
          coords.start2<-as.data.frame(sf::st_coordinates(start.node2))
          line.dat2[1,]<-coords.start2
        }
        if(end.good){
          end.node2<-network.nodes[sf::st_nearest_feature(x = end.node,y=network.nodes),]
          coords.end2<-as.data.frame(sf::st_coordinates(end.node2))
          line.dat2[nrow(line.dat2),]<-coords.end2
        }

        #update the network
        new.edges<-c(sf::st_geometry(network.edges),sf::st_cast(sf::st_combine(sf::st_as_sf(line.dat2,coords=1:2,crs=crs)),"LINESTRING"))
        network<-sfnetworks::as_sfnetwork(new.edges)
        network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
        network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))
      }else{
        nbad<-nbad+1
      }
    }
    if(nbad>0){warning("Unable to add ",nbad," out of ",nrow(edges)," lines to the network")}
  }

  # start by getting all the nodes set up
  if(missing(nodes)==F){
    if(class(nodes)[1]=="data.frame"){
      if(any((c("lon","lat")%in%names(nodes))==F)){stop("Invalid data frame for nodes. Must contain columns named 'lon' and 'lat'")}
      node.points<-sf::st_as_sf(nodes[,c("lon","lat")],coords=1:2,crs=crs)
      node.points<-sf::st_transform(node.points,crs = sf::st_crs(network))
    }
    if(class(nodes)[1]=="sf"){
      if(any(sf::st_geometry_type(nodes)!="POINT")){stop("Nodes must have POINT geometry")}
      node.points<-sf::st_transform(nodes,crs = sf::st_crs(network))
    }

    # check that the points are on the network
    dist.mat<-sf::st_distance(x = node.points,y = network.edges)
    good.points<-node.points[apply(dist.mat,MARGIN = 1,FUN = function(x){return(any(as.numeric(x)<=tolerance))}),]

    if(nrow(good.points)<nrow(node.points)){warning(paste0("Failed to snap ",nrow(node.points)-nrow(good.points)," to the network"))}

    # check if the nodes are redundant
    dist.mat2<-sf::st_distance(x = good.points,y = network.nodes)
    add.points<-good.points[apply(dist.mat2,MARGIN = 1,FUN = function(x){return(all(as.numeric(x)>tolerance2))}),]
    if(nrow(add.points)>0){
      network<-sfnetworks::st_network_blend(x = network,y = add.points)
    }
  }
  return(network)
}


#short function that removes internal nubs and any small terminal segments
#' Simplifies a sfnetwork
#'
#' This function removes extraneous nodes and small segments from a network that are often
#' leftover after other operations
#'
#' @param network A sfnetwork
#' @param root an optional sf POINT or other input capatible with RootNetwork; useful if the root gets removed
#' @param preserve an optional input that specifies nodes or segments that should not be altered
#' @param preserve.type either "nodes" or "edges", indicating which component preserve is referencing
#' @param tolerance A integer, segments with no connections below this size are removed
#' @param tolerance2 A fine scale tolerance for resolving overlapping edges; ought to be small
#' @param makeplot should the funciton output a plot showing before/after versions of the network
#'
#' @return A sfnetwork object
#' @export
#'
#' @examples
SimplifyNetwork<-function(network,root,preserve,preserve.type,tolerance=1000,tolerance2=1,makeplot=T){

  old.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
  old.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))

  #set aside nodes to preserve here, and blend them back in later
  if(missing(preserve)==F){

    # start with some checks
    preserve.class<-class(preserve)

    if(preserve.class=="logical"){
      if(missing(preserve.type)){
        if(length(preserve)==nrow(old.edges)){preserve.type<-"edges"}
        if(length(preserve)==nrow(old.nodes)){preserve.type<-"nodes"}
        if(length(preserve)%in%c(nrow(old.edges),nrow(old.nodes))==F){stop("Vector 'preserve' is wrong length")}
        if(nrow(old.edges)==nrow(old.nodes)){stop("Please specify 'preserve.type'.")}
      }
      if(preserve.type%in%c("nodes","edges")==F){stop("'preserve.type' must be either 'nodes' or 'edges'")}
      if(preserve.type=="nodes"){reserve.nodes<-old.nodes[which(preserve),]}
      if(preserve.type=="edges"){
        reserve.edges<-old.edges[which(preserve),]
        reserve.nodes<-old.nodes[sort(unique(c(reserve.edges$to,reserve.edges$from))),]
      }
    }
    if(preserve.class%in%c("numeric","integer")){
      if(missing(preserve.type)){
        warning("No input for preserve.type! Defaulting to 'edges'!")
        preserve.type<-"edges"
      }
      if(preserve.type=="nodes"){reserve.nodes<-old.nodes[preserve,]}
      if(preserve.type=="edges"){
        reserve.edges<-old.edges[preserve,]
        reserve.nodes<-old.nodes[sort(unique(c(reserve.edges$to,reserve.edges$from))),]
      }
    }
    if(preserve.class=="character"){
      col.matches<-apply(as.data.frame(old.edges),MARGIN = 2,FUN = function(x){return(sum(preserve%in%x))})
      b.col<-which.max(col.matches)
      if(length(preserve)>max(col.matches)){
        nomatch<-which(preserve%in%as.data.frame(shape)[,b.col]==F)
        warning("Could not find '",preserve[nomatch],"' using best column '",names(shape)[b.col],"'")

      }
      reserve.edges<-old.edges[which(as.data.frame(old.edges)[,b.col]%in%preserve),]
      reserve.nodes<-old.nodes[sort(unique(c(reserve.edges$to,reserve.edges$from))),]
    }
  }else{
    reserve.nodes<-old.nodes[0,]
  }


  # first simplify by combining segments
  new.edges<-sf::st_as_sf(sf::st_cast(sf::st_line_merge(sf::st_combine(old.edges)),"LINESTRING"))

  # but there may still be overlapping lines
  new.network0<-suppressWarnings(sfnetworks::as_sfnetwork(new.edges))
  new.network.edges0<-sf::st_as_sf(sfnetworks::activate(new.network0,"edges"))

  new.network.edges<-new.network.edges0[0,]
  already.done<-NULL
  for(i in 1:nrow(new.network.edges0)){
    active.edge<-new.network.edges0[i,]
    active.length<-as.numeric(sf::st_length(active.edge))

    matches<-which(new.network.edges0$from==active.edge$from & new.network.edges0$to==active.edge$to)
    matches<-matches[matches%in%already.done==F]

    if(length(matches)==1){
      new.network.edges<-rbind(new.network.edges,active.edge)
      already.done<-c(already.done,matches)
    }
    if(length(matches)>1){
      match.edges<-new.network.edges0[matches,]
      match.lengths<-as.numeric(sf::st_length(match.edges))

      # check whether the edges are distinct
      distinct.edges<-which(abs(match.lengths-active.length)>tolerance2)
      if(length(distinct.edges)==0){
        new.network.edges<-rbind(new.network.edges,active.edge)
        already.done<-c(already.done,matches)
      }
      if(length(distinct.edges)>0){
        new.network.edges<-rbind(new.network.edges,active.edge,match.edges[distinct.edges,])
        already.done<-c(already.done,matches)
      }
    }
  }
  new.network<-suppressWarnings(sfnetworks::as_sfnetwork(new.network.edges))
  new.network.edges<-sf::st_as_sf(sfnetworks::activate(new.network,"edges"))
  new.network.nodes<-sf::st_as_sf(sfnetworks::activate(new.network,"nodes"))

  # then check if there are any nubs to drop
  shorts<-which(as.numeric(sf::st_length(new.network.edges))<tolerance)
  ends<-which(new.network.edges$to%in%new.network.edges$from==F)
  drop<-shorts[shorts%in%ends]

  if(length(drop)>0){
    new.network2<-suppressWarnings(sfnetworks::as_sfnetwork(new.network.edges[-drop,]))
  }else{
    new.network2<-new.network
  }

  #simplify again
  new.network.edges2<-sf::st_as_sf(sfnetworks::activate(new.network2,"edges"))
  new.edges2<-sf::st_as_sf(sf::st_cast(sf::st_line_merge(sf::st_combine(new.network.edges2)),"LINESTRING"))
  out.network<-sfnetworks::as_sfnetwork(new.edges2)

  if(nrow(reserve.nodes)>0){
    out.network<-sfnetworks::st_network_blend(x = out.network,y = reserve.nodes)
  }

  if(missing(root)==F){
    out.network<-RootNetwork(network = out.network,root = root,tolerance = 0)
  }

  # plot it out
  if(makeplot){
    out.edges<-sf::st_as_sf(sfnetworks::activate(out.network,"edges"))
    out.nodes<-sf::st_as_sf(sfnetworks::activate(out.network,"nodes"))

    plot(network,col=2)
    plot(out.network,add=T)
  }
  return(out.network)
}


#' Extracts a basin or subbasin
#'
#' this function takes a shape and extracts the network above that node
#' Note that this function does not preserve information associated with the original shapefile
#' Future versions should allow the user to specify information to preserve, and allow for
#'
#' @param shape an sf object consisting of linestrings (not multilinestrings!)
#' @param root a sf point, or length 2 vector with xy coordinates
#' @param basins a vector of names or sf feature indicating the basin to be extracted
#' @param exclude vector indexes to be dropped, will attempt to match a single name
#' @param extract.inner logical; should nodes that are in between the basins and the root be extracted
#' @param attach.data logical; should old data fields be obtained; time consuming
#' @param upper a point or points to stop the extraction
#' @param upper.crs a crs string for the upper point, if different from shape
#'
#' @return A sfnetwork containing the specified subbasins
#' @export
#'
#' @examples
ExtractBasin<-function(shape,root,basins,exclude,extract.inner=F,attach.data=T,upper,upper.crs){

  if(inherits(shape,what = "sfnetwork")){shape<-sf::st_as_sf(sfnetworks::activate(shape,"edges"))}

  # check for some errors before we waste a lot of processing time
  if(missing(basins) & missing(upper)){stop("Must specify a basin to extract")}
  echeck<-ifelse(missing(exclude),F,is.character(exclude))
  if(echeck){
    exclude.matches<-apply(as.data.frame(shape),MARGIN = 2,FUN = function(x){return(sum(exclude%in%x))})
    if(sum(exclude.matches==max(exclude.matches))>1){stop("Multiple columns match exclude argument!")}
    if(sum(exclude.matches==max(exclude.matches))==0){stop("Zero columns match exclude argument!")}
  }

  # First off, drop anything in exclude
  if(missing(exclude)==F){

    if(is.logical(exclude)){
      drops<-which(exclude)
    }
    if(is.numeric(exclude) | is.integer(exclude)){
      drops<-exclude
    }
    if(is.character(exclude)){
      col.matches<-apply(as.data.frame(shape),MARGIN = 2,FUN = function(x){return(sum(exclude%in%x))})
      ex.col<-which.max(col.matches)
      drops<-which(as.data.frame(shape)[,ex.col]%in%exclude)
    }
    shape2<-shape[-drops,]
  }else{
    shape2<-shape
  }

  root.sf<-PointSetup(point = root,crs = sf::st_crs(shape))

  # Setup the network
  shape.net<-Makesfnetwork(shape = shape2,attach.data = F)
  shape.net2<-RootNetwork(network = shape.net,root=root.sf)

  shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape.net2,'nodes'))
  shape.edges<-sf::st_as_sf(sfnetworks::activate(shape.net2,'edges'))
  root.node.index<-which.min(sf::st_distance(x=root.sf,y=shape.nodes))

  # If upper points have been specified, drop everything above those points
  if(missing(upper)==F){
    if(missing(upper.crs)){upper.crs<-sf::st_crs(shape)}
    upper.points<-PointSetup(upper,crs=upper.crs)
    upper.points2<-sf::st_transform(x = upper.points,crs = sf::st_crs(shape))

    shape.net2<-sfnetworks::st_network_blend(x = shape.net2,y = upper.points2)
    shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape.net2,'nodes'))
    shape.edges<-sf::st_as_sf(sfnetworks::activate(shape.net2,'edges'))
    upper.index<-sf::st_nearest_feature(x = upper.points2,y = shape.nodes)
    all.paths<-sfnetworks::st_network_paths(x = shape.net2,from = root.sf)

    drop.nodes.index<-which(unlist(lapply(all.paths[[1]],FUN=function(x){return(any(x%in%upper.index))})))
    drop.nodes.index<-drop.nodes.index[drop.nodes.index%in%upper.index==F]

    drop.segs<-which(shape.edges$from%in%drop.nodes.index | shape.edges$to%in%drop.nodes.index)

    new.edges<-shape.edges[-drop.segs,]
    shape.net2<-sfnetworks::as_sfnetwork(new.edges)
    shape.nodes<-sf::st_as_sf(sfnetworks::activate(shape.net2,'nodes'))
    shape.edges<-sf::st_as_sf(sfnetworks::activate(shape.net2,'edges'))
    root.node.index<-which.min(sf::st_distance(x=root.sf,y=shape.nodes))
  }

  # The ideal method is to specify a basin to extract
  if(is.logical(basins)){
    b.segments<-shape[which(basins),]
  }
  if(is.numeric(basins) | is.integer(basins)){
    b.segments<-shape[basins,]
  }
  if(is.character(basins)){
    col.matches<-apply(as.data.frame(shape),MARGIN = 2,FUN = function(x){return(sum(basins%in%x))})
    b.col<-which.max(col.matches)
    if(length(basins)>max(col.matches)){
      nomatch<-which(basins%in%as.data.frame(shape)[,b.col]==F)
      warning("Could not find match for basin '",basins[nomatch],"' using best column '",names(shape)[b.col],"'")
      basins<-basin[-nomatch]
    }
    b.segments<-shape[which(as.data.frame(shape)[,b.col]%in%basins),]
  }

  # need to identify redundant segments
  b.roots.index<-rep(NA,nrow(b.segments))
  b.inner.index<-rep(NA,nrow(b.segments))

  # we locate the start and end point for each and determine which is the root
  for(b in 1:nrow(b.segments)){
    seg.points<-sf::st_cast(sf::st_geometry(b.segments[b,]),"POINT")
    pts<-c(seg.points[1],seg.points[length(seg.points)])
    pts.index<-sf::st_nearest_feature(x = pts,y = shape.nodes)

    p.costs<-sfnetworks::st_network_cost(x = shape.net,from = root.sf,to=pts,direction = "all")

    b.roots.index[b]<-pts.index[which.min(p.costs)]

    # now we go one up from the basin root
    pts.index2<-which(sf::st_geometry(shape.nodes)%in%sf::st_geometry(seg.points))
    p.costs2<-sfnetworks::st_network_cost(x = shape.net,from = root.sf,to=shape.nodes[pts.index2,],direction = "all")
    b.inner.index[b]<-pts.index2[order(p.costs2)[2]]
  }

  # we keep only the roots that are not also basin nodes
  b.roots.index<-b.roots.index[b.roots.index%in%b.inner.index==F]

  # we take all paths that link through inner nodes (if you include end nodes, you'll get the whole system)
  all.paths<-sfnetworks::st_network_paths(x = shape.net2,from = root.sf,to = shape.nodes)
  good.nodes.index<-which(unlist(lapply(all.paths[[1]],FUN=function(x){return(any(x%in%b.inner.index))})))

  # separately, we find the edges that connect the basins, and put that aside for a few steps
  basin.connectors<-sort(unique(unlist(all.paths[[2]][b.roots.index])))

  # then we construct a connector that links those roots and the main root
  if(extract.inner){
    connector.paths<-sfnetworks::st_network_paths(x = shape.net2,from = root.sf,to = b.roots.index)
    inner.connect.nodes0<-sort(unique(unlist(connector.paths[[1]])))
    inner.connect.nodes<-inner.connect.nodes0[inner.connect.nodes0%in%c(root.node.index,b.roots.index)==F]

    # function determines if a path hits an inner node before one of the root nodes
    good.connect.nodes<-which(unlist(lapply(all.paths[[1]],FUN=function(x){
      conn<-x[which(x%in%c(inner.connect.nodes,b.roots.index))]
      if(length(conn)==0){return(FALSE)}
      return(utils::tail(conn,1)%in%b.roots.index==F)})))

    good.nodes.index<-unique(c(good.nodes.index,good.connect.nodes))
  }

  # we extract anything that is attached to a valid inner node, and anything required to link the basin roots
  good.edges<-shape.edges[shape.edges$from%in%good.nodes.index | shape.edges$to%in%good.nodes.index | 1:nrow(shape.edges)%in%basin.connectors,]

  if(attach.data){
    good.edges<-AttachData(shape = good.edges,refshape = shape2)
  }

  plot(sf::st_geometry(good.edges))
  plot(sf::st_geometry(shape),add=T)
  plot(sf::st_geometry(good.edges),add=T,col=2,lwd=2)
  plot(sf::st_geometry(root.sf),add=T,pch=16,col=2)

  return(sfnetworks::as_sfnetwork(good.edges))
}


#' Attach Data to a network
#'
#' Take data contained in one sf LINESTRING object and attaches it to a network.
#' Values are averaged or combined where appropriate. Take care with ID numbers or coded data that may appear numeric,
#' but is not meant to be combined or averaged.
#'
#' @param shape A sfnetwork or sf object with LINESTRING geometry to attach data to
#' @param refshape A sf object with LINESTRING geometry that contains the desired data
#' @param weight.type interger indicating which weight type is to be used; see details
#' @param fields A character vector with the names of the columns to attach, or "all"
#' @param logical T/F determining how to handle logical values
#' @param tolerance A number used when matching lines that are not perfectly aligned
#'
#' @return A sfnetwork consisting of the original shape with data attached now
#' @export
#'
#' @examples
AttachData<-function(shape,
                     refshape,
                     weight.type=1,
                     fields="all",
                     logical=F,
                     tolerance=100){

  if(weight.type%in%c(1,2)==F){stop("Select weight.type '1' or '2'")}

  if(class(shape)[1]=="sfnetwork"){
    shape.net<-shape
    shape.sf<-sf::st_as_sf(sfnetworks::activate(shape,"edges"))
  }else{
    shape.net<-sfnetworks::as_sfnetwork(shape)
    shape.sf<-shape
  }
  shape.sf.union<-sf::st_union(shape.sf)

  if(class(refshape)[1]=="sfnetwork"){
    refshape.net<-refshape
    refshape.sf<-sf::st_as_sf(sfnetworks::activate(refshape,"edges"))
  }else{
    refshape.net<-sfnetworks::as_sfnetwork(refshape)
    refshape.sf<-refshape
  }

  if(sf::st_crs(shape)!=sf::st_crs(refshape)){
    refshape.sf<-sf::st_transform(refshape.sf,crs=sf::st_crs(shape))
  }
  refshape.sf<-refshape.sf[sf::st_intersects(x=refshape.sf,y=shape.sf.union,sparse = F),]

  shape.names<-names(shape.sf)
  ref.names<-names(refshape.sf)
  sf.name<-attr(refshape,"sf_column")

  if(fields[1]=="all"){
    transfer.names<-ref.names[tolower(ref.names)%in%tolower(shape.names)==F & ref.names!=sf.name]
  }else{
    transfer.names<-fields
    if(any(fields%in%ref.names==F)){
      badnames<-fields[fields%in%ref.names==F]
      stop(paste0("Column names not found: ",paste(badnames,collapse = ", ")))
    }
  }

  transfer.mat<-matrix(ncol=length(transfer.names),nrow=nrow(shape.sf),data=NA)
  colnames(transfer.mat)<-transfer.names
  transfer.dat<-as.data.frame(transfer.mat)

  #We'll also need to know the kind of data in each column
  ref.types<-unlist(sapply(as.data.frame(refshape.sf),class))[transfer.names]
  number.names<-names(ref.types)[which(ref.types%in%c("numeric","integer"))]
  char.names<-names(ref.types)[which(ref.types%in%c("character","factor"))]
  logic.names<-names(ref.types)[which(ref.types=="logical")]

  pb<-utils::txtProgressBar(min = 0, max = nrow(shape.sf), initial = 0, style=3)
  # start by pulling out each feature
  for(i in 1:nrow(shape.sf)){
    shape.points<-sf::st_as_sf(as.data.frame(sf::st_coordinates(shape.sf[i,])[,1:2]),
                               coords=1:2,crs=sf::st_crs(shape))

    end.points<-sf::st_coordinates(shape.sf[i,])[c(1,nrow(sf::st_coordinates(shape.sf[i,]))),]
    end.points.sf<-sf::st_as_sf(as.data.frame(end.points)[,1:2],coords=1:2,crs=sf::st_crs(shape))

    start.dists<-as.numeric(sf::st_length(sf::st_nearest_points(x=end.points.sf[1,],y=refshape.sf)))
    end.dists<-as.numeric(sf::st_length(sf::st_nearest_points(x=end.points.sf[2,],y=refshape.sf)))

    start.min<-min(start.dists)
    end.min<-min(end.dists)

    if(start.min<=tolerance & end.min<=tolerance){
      start.feature<-which(start.dists==start.min,)
      end.feature<-which(end.dists==end.min)

      # we use two types of weights: the first is the % of shape made up of each refshape feature
      # and the second is the % of each refshape covered by shape
      #the shape in contained in one feature then the wgt is one
      if(all(start.feature==end.feature) & (length(start.feature)==1 & length(end.feature)==1)){
        wgt.vec<-1
        alt.wgt<-as.numeric(sf::st_length(shape.sf[i,])/sf::st_length(refshape.sf[start.feature,]))
        all.features<-start.feature
      }else{
        all.features<-unique(sf::st_nearest_feature(shape.points,refshape.sf))

        ref.points<-sf::st_as_sf(as.data.frame(sf::st_coordinates(refshape.sf[all.features,])[,1:2]),
                                 coords=1:2,crs=sf::st_crs(shape))
        snap.points<-sf::st_nearest_feature(shape.points,ref.points)

        new.line<-sf::st_cast(sf::st_combine(ref.points[snap.points,]),"LINESTRING")
        if(length(unique(snap.points))>1){
          snap.intersect<-sf::st_intersection(new.line,refshape.sf[all.features,])
          wgt.vec<-as.numeric(sf::st_length(snap.intersect)/sum(sf::st_length(snap.intersect)))
          alt.wgt<-as.numeric(sf::st_length(snap.intersect))/as.numeric(sf::st_length(refshape.sf[all.features,]))
        }else{
          wgt.vec<-0  # these cases handle when there is only a single point of overlap
          alt.wgt<-0
        }
      }

      if(weight.type==1){use.wgt<-wgt.vec}
      if(weight.type==2){use.wgt<-alt.wgt}

      # numeric values are given a weighted average
      if(length(number.names)>0){
        number.data<-as.data.frame(refshape.sf)[all.features,number.names]
        transfer.dat[i,number.names]<-apply(X=as.data.frame(number.data)*use.wgt,FUN="sum",MARGIN=2)
      }

      # characters and factors are concatenated, depending on weight
      # Have tried to handle a variety of error types
      if(length(char.names)>0){
        char.data<-as.data.frame(as.data.frame(refshape.sf)[all.features,char.names])
        char.data$use.wgt<-use.wgt
        out.char<-vector(length = length(char.names))

        # determine if characters are unique
        for(j in 1:length(char.names)){
          unique.char<-unique(char.data[,j])

          if(length(unique.char)==0){unique.char<-NA}
          if(length(unique.char)==1){
            out.char<-unique.char
          }else{
            char.wgts<-stats::aggregate(char.data$use.wgt,by=list(as.factor(char.data[,j])),FUN=sum)
            out.char<-paste(char.wgts$Group.1[char.wgts$x>=.25],collapse = "*")
          }
          transfer.dat[i,char.names[j]]<-out.char
        }
      }
      #untested : logical values are true only if the entire reach is true
      if(length(logic.names)>0){
        logic.data<-as.data.frame(refshape.sf)[,logic.names]
        transfer.dat[i,logic.names]<-apply(X=as.integer(logic.data),MARGIN=2,FUN="floor")
      }
    }
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)
  out.shape.sf<-cbind(shape.sf,transfer.dat)

  # return the shape in the same format as it was given
  if(class(shape)[1]=="sfnetwork"){
    return(sfnetworks::as_sfnetwork(out.shape.sf))
  }else(
    return(out.shape.sf)
  )
}


#' Prunes a network based on various criteria
#'
#' New and improved version of prune network, it handles a larger variety of inputs and is better implemented
#' exclude and match can be used in a variety of ways, and the function can sometimes handle them at the same time
#' But it is often better to call the function multiple times, using different arguments
#'
#' @param network A sfnetwork to be pruned
#' @param root a sf point or coordinates for the root
#' @param exclude a logical vector or vector of values to exclude
#' @param match a sf object with points or lines
#' @param tolerance a tolerance for matching points to the network
#' @param plot logical, should a plot be made
#'
#' @return A sfnetwork pruned as specified
#' @export
#'
#' @examples
PruneNetwork<-function(network,root,exclude,match,tolerance=100,plot=T){

  root.sf<-PointSetup(point = root,crs = sf::st_crs(network))

  # check the network setup; may want to add more checks later
  if(sfnetworks::is.sfnetwork(network)){
    network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
    network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))
  }else{
    if(class(network)[1]%in%c("sf","sfc_LINESTRING")){
      if(all(sf::st_geometry_type(network)=="LINESTRING")){
        network.edges<-sf::st_as_sf(network)
      }else{
        stop("Network argument must be a sfnetwork or sf with LINESTRING geometry")
      }
    }else{
      stop("Network argument must be a sfnetwork or sf with LINESTRING geometry")
    }
  }

  # Flexible, but complicated to handle different data types
  if(missing(exclude)==F){
    exclude.class<-class(exclude)[1]
  }else{
    exclude.class<-"none"
  }
  if(missing(match)==F){
    match.class<-class(match)[1]
    if(strsplit(match.class,split="_")[[1]][1]=="sfc"){
      match<-sf::st_as_sf(match)
      match.class<-"sf"
    }

  }else{
    match.class<-"none"
  }
  if(exclude.class=="sf"){stop("SF objects are invalid for exclude")}
  # If not using sf objects, then exclude and match are redundant, and can use the same code with a sign change
  switcheroo<-F
  if(match.class%in%c("sf","none")==F){
    if(exclude.class!="none"){stop("Redundant criteria between exclude and match arguments; please choose one")}
    if(exclude.class=="none"){
      switcheroo<-T
      exclude.class<-match.class
      match.class<-"none"
      exclude<-match
    }
  }

  # this section handles exclusions and matches by logical vectors, character string, or index number
  if(exclude.class=="logical"){
    if(length(exclude)!=nrow(network.edges)){stop("Exclude does not match length of network")}
    if(any(is.na(exclude))){stop("Invalid value (NA) in exclude")}
    drop<-which(exclude)
  }
  if(exclude.class=="character"){
    ex.no<-apply(as.data.frame(network.edges),MARGIN = 2,FUN = function(x){return(sum(exclude%in%x))})
    if(max(ex.no)==0){stop("Could not find exclude criterion; Check spelling?")}
    ex.col<-which.max(ex.no)
    drop<-which(as.data.frame(network.edges)[,ex.col]%in%exclude)
  }
  if(exclude.class%in%c("integer","numeric")){
    if(any(exclude>nrow(network.edges)) | any(exclude<1)){stop("Invalid rows selected")}
    drop<-exclude
  }
  if(exclude.class=="none"){
    good.edges1<-network.edges
  }else{
    if(switcheroo){
      good.edges1<-network.edges[drop,]
    }else{
      good.edges1<-network.edges[-drop,]
    }
  }

  # this section prunes based on an sf object
  if(match.class=="sf"){
    if(any(sf::st_geometry_type(match)%in%c("POINT","LINESTRING")==F)){stop("Invalid geometry type in match argument")}
    #if(any(sf::st_geometry_type(match)))){stop("Multiple geometry types in match argument")}

    network2<-Makesfnetwork(good.edges1,attach.data = F)
    network2<-RootNetwork(network2,root = root.sf)
    network.edges2<-sf::st_as_sf(sfnetworks::activate(network2,"edges"))
    network.nodes2<-sf::st_as_sf(sfnetworks::activate(network2,"nodes"))
    network.points<-sf::st_sf(sf::st_cast(sf::st_geometry(network.edges2),"POINT"))

    # Let's generalize it to handle geometry collections
    for(i in 1:nrow(match)){

      # turn each feature into an endpoint
      if(sf::st_geometry_type(match)[i]=="POINT"){
        net.pts0<-network.points[sf::st_nearest_feature(x = match[i,],network.points),]
        net.pts<-net.pts0[as.numeric(sf::st_distance(match[i,],net.pts0,by_element = T))<=tolerance,]
      }
      if(sf::st_geometry_type(match)[i]=="LINESTRING"){
        coords<-unique(sf::st_coordinates(match[i,])[,1:2])
        coords.sf<-sf::st_as_sf(as.data.frame(coords),coords=1:2,crs=sf::st_crs(network))

        net.pts0<-network.points[sf::st_nearest_feature(x = coords.sf,network.points),]
        net.pts<-net.pts0[as.numeric(sf::st_distance(x=coords.sf,y=net.pts0,by_element = T))<=tolerance,]
      }

      if(i==1){net.points<-net.pts}
      if(i>1){net.points<-rbind(net.points,net.pts)}
    }

    # # start by looking for endpoints
    # if(sf::st_geometry_type(match)[1]=="POINT"){
    #   net.points0<-network.points[sf::st_nearest_feature(x = match,network.points),]
    #   net.points<-net.points0[as.numeric(sf::st_distance(match,net.points0,by_element = T))<=tolerance,]
    # }
    # if(sf::st_geometry_type(match)[1]=="LINESTRING"){
    #   coords<-unique(sf::st_coordinates(match)[,1:2])
    #   coords.sf<-sf::st_as_sf(as.data.frame(coords),coords=1:2,crs=sf::st_crs(network))
    #
    #   net.points0<-network.points[sf::st_nearest_feature(x = coords.sf,network.points),]
    #   net.points<-net.points0[as.numeric(sf::st_distance(x=coords.sf,y=net.points0,by_element = T))<=tolerance,]
    # }

    # now that we have endpoints, we'll start by making a new network and blending the points in
    blend.points<-unique(sf::st_as_sf(c(sf::st_geometry(net.points),sf::st_geometry(network.nodes2))))
    network3<-sfnetworks::st_network_blend(network2,blend.points)

    # in order to accomodate loops, we will need to add additional code here
    # code will check for loops, find the subroot of any given segment, and ensure that
    # stream direction is set up appropriately
    network3<-RootNetwork(network3,root.sf)

    edge.paths<-sfnetworks::st_network_paths(network3,root,net.points)[[2]]
    good.edges2<-sf::st_as_sf(sfnetworks::activate(network3,"edges"))[unique(unlist(edge.paths)),]
  }else{
    good.edges2<-good.edges1
  }

  # removed the simplification step, I think its better if the user calls it separately
  # good.edges3<-sf::st_cast(sf::st_line_merge(sf::st_combine(good.edges2)),"LINESTRING")
  out.network<-RootNetwork(network = sfnetworks::as_sfnetwork(good.edges2),root = root.sf,tolerance = 0)

  if(plot){
    plot(sf::st_geometry(network.edges))
    plot(sf::st_geometry(good.edges2),col=2,add=T)
  }

  return(out.network)
}


#' Divide a network into reaches
#'
#' This function takes a network and divides it into distinct reaches based on a target size
#' if you want only set boundaries, use setnodes and set targetsize very high
#' The setnodes argument can be used to ensure certain reach boundaries (not implemented yet)
#'
#' @param network A sfnetwork
#' @param targetsize a desired size to use for reaches, set to NA to use only set nodes
#' @param habitat a logical vector that identifies areas of habitat (T) vs non-habitat (F)
#' @param setnodes A sf POINTS object that contains nodes to add to the shape
#'
#' @return A new sfnetwork object divided into reaches of approximately the correct size
#' @export
#'
#' @examples
AssignReaches<-function(network,targetsize=NA,habitat,setnodes){

  if(inherits(network,"sfnetwork")==F){
    print("Converting shape to sfnetwork")
    network<-Makesfnetwork(network)
  }
  network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
  network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))

  # Set up the habitat distinctions
  if(missing(habitat)){
    habitat2<-rep(TRUE,nrow(network.edges))
  }else{
    if(is.logical(habitat)){
      if(length(habitat)!=nrow(network.edges)){stop("Length of 'habitat' does not match # of edges")}
      habitat2<-habitat
    }
    if(is.integer(habitat)){
      if(any(habitat>nrow(network.edges)) | any(habitat<1)){stop("Invalid index in 'habitat'")}
      habitat2<-rep(FALSE,nrow(network.edges))
      habitat2[habitat]<-TRUE
    }
    if(is.logical(habitat)==F & is.integer(habitat)==F){stop("Invalid input for 'habitat'")}
  }
  network.edges$habitat<-habitat2

  # If using set reach nodes, blend them in here
  if(missing(setnodes)==F){

    # start by checking nodes
    if(inherits(setnodes,"sf")==F){stop("Setnodes argument must be an sf object with POINT geometry")}
    if(any(sf::st_geometry_type(setnodes)!="POINT")){stop("Setnodes argument must be an sf object with POINT geometry")}
    if(sf::st_crs(setnodes)!=sf::st_crs(network.edges)){stop("CRS for setnodes does not match network")}

    set.edges0<-sfnetworks::as_sfnetwork(network.edges)
    set.edges<-sfnetworks::st_network_blend(x = set.edges0,y = setnodes)
    network.edges<-sf::st_as_sf(sfnetworks::activate(set.edges,"edges"))
  }

  if(is.na(targetsize)){
    out.edges<-network.edges
  }else{
    # we want to identify the best size before diving into the loop
    reach.lengths<-sf::st_length(network.edges)
    size.range<-targetsize*90:110/100
    targetsize.vec<-rep(size.range,each=nrow(network.edges))

    mod.mat<-matrix(nrow=21,ncol=nrow(network.edges),data=reach.lengths%%targetsize.vec,byrow=T)
    if(any(network.edges$habitat==F)){mod.mat[,which(network.edges$habitat==F)]<-NA}
    bestsize<-size.range[which.min(apply(mod.mat,1,sum,na.rm=T))]
    nreaches.vec<-as.integer(round(reach.lengths/bestsize))
    nreaches.vec[nreaches.vec==0]<-1

    print(paste0("Using reach size of ",bestsize))
    for(i in 1:nrow(network.edges)){

      reach<-network.edges[i,]

      if(nreaches.vec[i]==1 | network.edges$habitat[i]==F){
        out.edges0<-reach
      }else{
        reach.nodes0<-sf::st_line_sample(x = reach,sample = c(0,1/nreaches.vec[i]*1:(nreaches.vec[i]-1),1))
        reach.nodes<-sf::st_as_sf(sf::st_cast(reach.nodes0,"POINT"))

        reach.net<-suppressWarnings(sfnetworks::st_network_blend(sfnetworks::as_sfnetwork(reach),reach.nodes))
        out.edges0<-sf::st_as_sf(sfnetworks::activate(reach.net,"edges"))
      }

      sf::st_geometry(out.edges0)<-"geometry"
      if(i==1){
        out.edges<-out.edges0
      }else{
        out.edges<-rbind(out.edges,out.edges0)
      }
    }
  }
  out.edges$Reach<-1L:nrow(out.edges)

  return(sfnetworks::as_sfnetwork(out.edges))
}


#' Check a network and correct common problems
#'
#' this function checks that the direction of the stream network is away from the root
#' and fixes any problems, also adds a column named "parent" to the shape and the distance
#' from the midpoint of the parent to the midpoint of the child
#'
#' @param network A sfnetwork object
#' @param root A sf POINT or vector with coordinates for the root
#'
#' @return A sfnetwork with some additional data attached
#' @export
#'
#' @examples
CheckNetwork<-function(network,root){

  network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))
  root.sf<-PointSetup(point = root,crs = sf::st_crs(network))

  # first, check that all stream segments are going the right direction
  root.node.index<-sf::st_nearest_feature(root.sf,sf::st_as_sf(sfnetworks::activate(network,"nodes")))
  root.distances<-as.vector(sfnetworks::st_network_cost(network,root.node.index,direction="all"))
  reversed<-which(root.distances[network.edges$from]>root.distances[network.edges$to])

  old.network.edges<-network.edges
  if(length(reversed)>0){
    old.network.edges$from[reversed]<-network.edges$to[reversed]
    old.network.edges$to[reversed]<-network.edges$from[reversed]

    for(i in 1:length(reversed)){
      old.line<-old.network.edges$geometry[reversed[i]]
      old.network.edges$geometry[reversed[i]]<-sf::st_reverse(old.line)
    }
  }

  # assign the root
  old.network.edges$root<-old.network.edges$from==root.node.index

  #Now, let renumber the Reachs so that they procede in an easy to interpret order
  new.network<-sfnetworks::as_sfnetwork(old.network.edges)
  new.edges<-sf::st_as_sf(sfnetworks::activate(new.network,"edges"))
  new.nodes<-sf::st_as_sf(sfnetworks::activate(new.network,"nodes"))

  root.index<-new.edges$from[which(new.edges$root)]
  end.nodes<-new.edges$to[new.edges$to%in%new.edges$from==F]

  #need to make a table that converts the old order to the new order, then we can grab values from old data frame
  network.table<-data.frame(oldid=new.edges$Reach,newids=0)
  costs<-sfnetworks::st_network_cost(x = new.network,from = root.index,to = end.nodes)
  costs2<-sort(costs,decreasing=T)
  end.nodes2<-end.nodes[order(costs,decreasing=T)]

  for(i in 1:length(end.nodes2)){
    start.val<-max(network.table$newids)+1

    path<-sfnetworks::st_network_paths(x = new.network,from = root.index,to = end.nodes2[i])$node_paths[[1]]
    path2<-path[path!=root.index]

    match.vec<-match(path2,new.edges$to)
    match.vec2<-match.vec[network.table$newids[match.vec]==0]

    network.table$newids[match.vec2]<-start.val:(start.val+length(match.vec2)-1)
  }
  new.edges$Reach<-network.table$newid[match(new.edges$Reach,network.table$oldid)]
  new.edges<-new.edges[order(new.edges$Reach),]

  # now can compute the parents, though this is going to be complicated by habitat/nonhabitat
  new.edges$parent<-new.edges$Reach[match(new.edges$from,new.edges$to)]

  new.parent.vec<-rep(NA,nrow(new.edges))
  parent.distance.vec<-rep(NA,nrow(new.edges))

  # the root must be marked habitat
  new.edges$habitat[which(new.edges$root)]<-TRUE

  # now account for habitat
  for(i in 1:nrow(new.edges)){
    if(new.edges$root[i]==F){
      parent<-new.edges$parent[i]
      # if the parent is habitat, things are easy
      if(new.edges$habitat[parent]){
        new.parent.vec[i]<-parent
        parent.distance.vec[i]<-.5*sf::st_length(new.edges)[i]+.5*sf::st_length(new.edges)[parent]
      }else{
        # loop through to parent is easier than converting to network
        bad.parent<-T
        path.segs<-parent
        active.seg<-parent
        while(bad.parent){
          new.parent<-new.edges$parent[active.seg]
          if(new.edges$habitat[new.parent]){
            bad.parent<-F
          }else{
            path.segs<-c(path.segs,new.parent)
            active.seg<-new.parent
          }
        }
        new.parent.vec[i]<-new.parent
        parent.distance.vec[i]<-.5*sf::st_length(new.edges)[i]+.5*sf::st_length(new.edges)[new.parent]+
          sum(sf::st_length(new.edges)[path.segs])
      }
    }
  }

  new.edges$parent<-new.parent.vec
  new.edges$parent.distance<-as.numeric(parent.distance.vec)

  return(sfnetworks::as_sfnetwork(new.edges))
}


#This is a new function that separates out the survey tracks
#' Makes sf LINESTRINGs representing survey transects
#'
#' This function uses coordinates to locate the start and stop points for a survey
#' and then draws a line that follows the stream contours between those points.
#' If multiple paths are available, it will choose the shortest, so it's best to eliminate
#' multichannel sections ahead of time
#'
#' @param shape A sf network or sf object with LINESTRING geometry
#' @param surveys A data frame with coordinates and any other desired information
#' @param surveys.crs a crs object or code for the survey coordinates
#' @param save.col character vector with names of any columns to save, or "all"
#' @param maxdist a maximum distance for when the coordinates are compared to the shape
#' @param survey.coords vector of names or column numbers; order matters
#'
#' @return A sf object with LINESTRINGS and data representing each survey
#' @export
#'
#' @examples
MakeSurveyTracks<-function(shape, surveys,surveys.crs="wgs84",save.col="all",
                           maxdist=25, survey.coords=c("startlon","startlat","endlon","endlat")){

  if(class(shape)[1]=="sfnetwork"){
    shape.net<-shape
    shape.edges<-sf::st_as_sf(sfnetworks::activate(shape,"edges"))
  }else{
    shape.net<-Makesfnetwork(shape,attach.data=F)
    shape.edges<-sf::st_as_sf(shape)
  }

  if(save.col[1]=="all"){
    save.col<-names(surveys)[names(surveys)%in%survey.coords==F]
  }

  # need to turn the surveys in line objects
  shape.union<-sf::st_geometry(shape.edges)
  if(nrow(shape.edges)>1 | any(sf::st_geometry_type(shape.edges)=="MULTILINESTRING")){
    shape.union<-sf::st_line_merge(sf::st_union(shape.edges))
  }

  if(sf::st_geometry_type(shape.union)=="MULTILINESTRING"){
    shape.union<-sf::st_cast(shape.union,"LINESTRING")
  }

  start.points0<-sf::st_as_sf(surveys,coords=survey.coords[1:2],crs=surveys.crs)
  start.points<-sf::st_transform(start.points0,crs=sf::st_crs(shape))
  end.points0<-sf::st_as_sf(surveys,coords=survey.coords[3:4],crs=surveys.crs)
  end.points<-sf::st_transform(end.points0,crs=sf::st_crs(shape))

  # Setup the output objects
  out.shape<-sf::st_as_sf(shape.union[0])
  sf::st_geometry(out.shape)<-"geometry"
  good.surveys<-NULL
  bad.coords<-NULL

  # loop through and make a feature collection of the survey lines
  pb<-utils::txtProgressBar(min = 0, max = nrow(start.points), initial = 0, style=3)
  for(i in 1:nrow(start.points)){
    test.point<-rbind(start.points[i,"geometry"],end.points[i,"geometry"])

    # first check that its on the network
    test.dists<-c(min(sf::st_distance(test.point[1,],shape.union)),
                  min(sf::st_distance(test.point[2,],shape.union)))

    # now structured so that if either point is out of area, discard it
    if(all(as.numeric(test.dists)<=maxdist)){

      start<-NA
      end<-NA
      # snap the points to the network, the careful way
      nearest1<-sf::st_nearest_points(test.point[1,],shape.union)
      pt1<-nearest1[which.min(sf::st_length(nearest1))]
      if(as.numeric(sf::st_length(pt1))<=maxdist){
        start<-sf::st_as_sf(as.data.frame(sf::st_coordinates(pt1))[2,1:2],coords=1:2,
                            crs=sf::st_crs(shape.union))
        start.seg<-shape.union[which.min(sf::st_length(nearest1))]
      }
      nearest2<-sf::st_nearest_points(test.point[2,],shape.union)
      pt2<-nearest2[which.min(sf::st_length(nearest2))]
      if(as.numeric(sf::st_length(pt2))<=maxdist){
        end<-sf::st_as_sf(as.data.frame(sf::st_coordinates(pt2))[2,1:2],coords=1:2,
                          crs=sf::st_crs(shape.union))
        end.seg<-shape.union[which.min(sf::st_length(nearest2))]
      }

      # check for that very annoying rounding error
      start.ok<-sf::st_intersects(start,start.seg,sparse = F)[1,1]
      end.ok<-sf::st_intersects(end,end.seg,sparse = F)[1,1]

      # if present, fix it the error
      start2<-start
      end2<-end
      if(!start.ok){
        start.seg.points<-sf::st_cast(start.seg,"POINT")
        start.dists<-sf::st_distance(start,start.seg.points)
        if(as.numeric(min(start.dists))<=maxdist){
          start2<-sf::st_as_sf(start.seg.points[which.min(start.dists)])
        }else{
          # its possible that the fixed point is further away than the desired tolerance
          # in which case, we can resample the line segment to find a closer point
          start.seg.points2<-sf::st_cast(sf::st_line_sample(start.seg,density=1/maxdist),"POINT")
          start.dists2<-sf::st_distance(start,start.seg.points2)
          start2<-sf::st_as_sf(start.seg.points2[which.min(start.dists2)])
        }
      }
      if(!end.ok){
        end.seg.points<-sf::st_cast(end.seg,"POINT")
        end.dists<-sf::st_distance(end,end.seg.points)
        if(as.numeric(min(end.dists))<=maxdist){
          end2<-sf::st_as_sf(end.seg.points[which.min(end.dists)])
        }else{
          end.seg.points2<-sf::st_cast(sf::st_line_sample(end.seg,density=1/maxdist),"POINT")
          end.dists2<-sf::st_distance(end,end.seg.points2)
          end2<-sf::st_as_sf(end.seg.points2[which.min(end.dists2)])
        }
      }

      # fix a problem where the geometry names are different, (there's prob a better way to do this)
      sf::st_geometry(start2)<-"geometry"
      sf::st_geometry(end2)<-"geometry"

      # now that we are sure both points are exactly on the stream line, we can add them to the network
      shape.net2<-sfnetworks::st_network_blend(x=shape.net,rbind(sf::st_as_sf(start2),sf::st_as_sf(end2)))
      edge.path<-sfnetworks::st_network_paths(x=shape.net2,from=sf::st_as_sf(start),to=sf::st_as_sf(end),mode="all")
      new.shape<-sf::st_as_sf(sfnetworks::activate(shape.net2,"edges"))[edge.path$edge_paths[[1]],]

      if(nrow(new.shape)>0){
        new.shape<-sf::st_as_sf(sf::st_line_merge(sf::st_combine(new.shape)))
        sf::st_geometry(new.shape)<-"geometry"
        out.shape<-rbind(out.shape,new.shape)
        good.surveys<-c(good.surveys,i)
      }else{
        bad.coords<-c(bad.coords,i)
      }
    }
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)

  if(length(bad.coords)>0){
    print(paste0("Warning: Removed ",length(bad.coords), " surveys due to bad coordinates!"))
    print(paste0("Effected survey indices : ",paste(bad.coords,collapse = ", ")))
  }

  # Attach the other data and return
  out.shape<-cbind(out.shape,surveys[good.surveys,save.col])
  sf::st_geometry(out.shape)<-"geometry"
  return(out.shape)
}


#' Guides the user in how to fix loops
#'
#' This function will guide the user through eliminating any loops in the network
#' using the command line. Note that for vast, a rooted tree is required,
#'so no reach can have multiple parents
#'
#' @param network A sfnetwork
#' @param root A sf object or coordinates for the root
#' @param guides optional list of sf objects that may assist decision making
#'
#' @return A sfnetwork with user specified sections removed
#' @export
#'
#' @examples
FixLoops<-function(network,root,guides=NULL){

  network.nodes<-sf::st_as_sf(sfnetworks::activate(network,"nodes"))
  network.edges<-sf::st_as_sf(sfnetworks::activate(network,"edges"))

  if(is.null(guides)==F){
    if(is.null(names(guides))){names(guides)<-paste0("Guide ",1:length(guides))}
    if(class(guides)!="list"){stop("The guides argument must be a list")}
  }
  root.sf<-PointSetup(point = root,crs = sf::st_crs(network))

  print("Identifying Loops")
  loop.vec<-rep(F,nrow(network.edges))
  pb<-utils::txtProgressBar(min = 0, max = nrow(network.edges), initial = 0, style=3)
  for(i in 1:nrow(network.edges)){

    network.edges2<-network.edges[-i,]

    #first, check if a point is still on the network
    intersect.mat<-sf::st_intersects(network.nodes,network.edges2)
    all.points<-!any(unlist(lapply(intersect.mat,sum))==0)

    # now check if every point still has a valid path
    if(all.points){
      loop.vec[i]<-!any(is.infinite(sfnetworks::st_network_cost(x = sfnetworks::as_sfnetwork(network.edges2),
                                                                from = root.sf,direction = "all")))
    }
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)
  if(all(loop.vec==F)){stop("No loops detected")}

  # now we need to group the loops together
  loop.index<-which(loop.vec)

  bounds<-sf::st_bbox(network.edges)
  xspan<-(bounds$xmax-bounds$xmin)*.2
  yspan<-(bounds$ymax-bounds$ymin)*.2
  temp.edges<-network.edges
  temp.edges$ID<-1:nrow(temp.edges)

  dists<-sf::st_distance(network.edges[loop.index,],root.sf)

  needs.doing<-temp.edges$ID[loop.index[order(dists)]]
  while(length(needs.doing>0)){

    active<-needs.doing[1]
    midpoint<-sf::st_coordinates(sf::st_line_sample(x = network.edges[active,],sample=.5))

    taking.action<-T
    while(taking.action){

      plot(sf::st_geometry(network.edges)[active],xlim=c(midpoint[1]-xspan/2,midpoint[1]+xspan/2),
           ylim=c(midpoint[2]-yspan/2,midpoint[2]+yspan/2))
      plot(sf::st_geometry(network.edges),add=T,col=7,lwd=2)
      plot(sf::st_geometry(temp.edges),add=T,col="gray70",lwd=2)
      plot(sf::st_geometry(network.edges)[needs.doing],add=T,lwd = 3)
      plot(sf::st_geometry(network.edges)[active],add=T,lwd = 4,col=3)

      if(length(guides)>0){
        for(g in 1:length(guides)){
          plot(sf::st_geometry(guides[[g]]),add=T,col=c(2,6,5,7)[[g]],lwd=.5)
        }
        guide.col<-c("red","purple","blue")[1:length(guides)]
      }else{
        guide.col<-NULL
      }
      graphics::legend("topright",legend=c("Basin","Loop Segments","Active Segment", "Cut Segments",names(guides)),lty=1,
                       col = c("gray70","black","green","yellow",guide.col))

      action0<-readline("Choose an action for the red highlighted segment: KEEP/REMOVE/TRUNCATE/zoom'IN'/zoom'OUT'")
      if(toupper(action0)%in%c("K","KEEP")){action<-"KEEP"}
      if(toupper(action0)%in%c("R","RE","RMOVE","RM","REMOVE")){action<-"REMOVE"}
      if(toupper(action0)%in%c("T","TR","TRUNC","TRUN","TRUNCATE")){action<-"TRUNCATE"}
      if(toupper(action0)%in%c("I","IN","ZOOMIN")){action<-"IN"}
      if(toupper(action0)%in%c("O","OUT","ZOOMOUT")){action<-"OUT"}

      if(action%in%c("IN","OUT")){
        scale<-ifelse(action=="IN",2,.5)
        xspan<-xspan/scale
        yspan<-yspan/scale
      }
      if(action=="KEEP"){
        needs.doing<-needs.doing[-1]
        taking.action<-F
      }
      if(action=="REMOVE"){
        needs.doing<-needs.doing[-1]
        taking.action<-F
        temp.edges<-temp.edges[temp.edges$ID!=active,]
      }
      if(action=="TRUNCATE"){
        needs.doing<-needs.doing[-1]
        taking.action<-F
        old.coords<-sf::st_coordinates(temp.edges[temp.edges$ID==active,])[,1:2]
        new.coords<-old.coords[1:(nrow(old.coords)*.9),]

        new.edge0<-sf::st_as_sf(as.data.frame(new.coords),coords=1:2,crs=sf::st_crs(temp.edges))
        sf::st_cast(sf::st_combine(new.edge0),"LINESTRING")
        sf::st_geometry(temp.edges)[temp.edges$ID==active]<-sf::st_cast(sf::st_combine(new.edge0),"LINESTRING")
      }
    }

    # if we changed the network, run some checks,
    # need to trim off any parts that are now cut from the root
    if(action%in%c("TRUNCATE","REMOVE")){
      temp.network<-sfnetworks::as_sfnetwork(temp.edges)
      temp.edges2<-sf::st_as_sf(sfnetworks::activate(temp.network,"edges"))

      temp.costs<-sfnetworks::st_network_cost(x = temp.network,from = root.sf,direction = "all")[1,]
      bad.nodes<-which(is.infinite(temp.costs))
      if(length(bad.nodes)>0){
        bad.ids<-temp.edges2$ID[temp.edges2$from%in%bad.nodes & temp.edges2$to%in%bad.nodes]

        temp.edges<-temp.edges[temp.edges$ID%in%bad.ids==F,]
        needs.doing<-needs.doing[needs.doing%in%bad.ids==F]
      }
    }
  }

  out.network<-Makesfnetwork(temp.edges)
  out.network2<-RootNetwork(out.network,root.sf)
  return(out.network)
}


# New version assumes that the MakeSurveyTracks function has already been used, makes things simpler
#' Combines survey tracks, counts, and a network into a data frame
#'
#' This function takes an sfnetwork or sf shape, a set of surveys, and optionally a set of
#' georeferenced redds observations, and returns a data set with the effort assigned to each reach
#' and the count types associated with that observation date.
#'
#' @param shape A sf object with LINESTRING geometries, or a sfnetwork
#' @param georedds An optional dataframe with coordinates for each redd observation
#' @param georedds.type An optional name for the column with classification of redds
#' @param georedds.coords A character vector with the column names for lon/lat in georedds
#' @param reddcrs A crs for the coordinates in georedds
#' @param tolerance.redds A maximum distance for snapping redds to the shape
#' @param tolerance.surveys A maximum distance for snapping surveys to the shape (should be small)
#' @param surveys A sf object with LINESTRINGs representing surveys
#' @param survey.redds optional name of a column with survey counts
#'
#' @return a dataframe with dates, counts, and effort for each surveyed location
#' @export
#'
#' @examples
AssembleReddData<-function(shape,                 # shape or network organized into reaches
                           georedds=NA,           # data frame with a redd in each row
                           georedds.type=NA,      # name of a column to sort redds by, counts are then kept separately
                           georedds.coords=NA,    # vector of length two with the names of the lon/lat columns
                           reddcrs="wgs84",       # crs for the redd coordinates
                           tolerance.redds=250,   # maximum distance for matching redds to a reach
                           tolerance.surveys=10,  # maximum distance for matching surveys to a reach; should be small
                           surveys,               # data frame containing
                           survey.redds=NA){      # name for the column with redd counts

  if(class(shape)[1]=="sfnetwork"){
    shape.edges<-sf::st_as_sf(sfnetworks::activate(shape,"edges"))
  }else(
    shape.edges<-shape
  )

  #remove some surveys
  if(is.na(survey.redds)==F){
    good.surveys<-surveys[is.na(as.data.frame(surveys)[,survey.redds])==F,]
  }else{
    good.surveys<-surveys
  }
  good.surveys$tempID<-1:nrow(good.surveys)

  # some setup
  surveys.dat<-as.data.frame(good.surveys)
  day.col<-which(tolower(names(surveys.dat))=="day")
  year.col<-which(tolower(names(surveys.dat))=="year")
  if(length(year.col)!=1 | length(day.col)!=1){stop("Please ensure your survey data contains columns named 'Year' and 'Day'")}

  # Setup any georeferenced redds, should improve handling of dates so that it's not just assumed to be in
  # correct format
  if(is.data.frame(georedds)){
    print("Formatting georeferenced counts")
    if(inherits(georedds,"sf")){
      georedds.sf<-georedds
    }else{
      if(is.na(georedds.coords)[1]){
        lon.col<-which(tolower(names(georedds))%in%c("longitude","lon","long","x"))[1]
        lat.col<-which(tolower(names(georedds))%in%c("latitude","lat","y"))[1]
      }else{
        lon.col<-georedds.coords[1]
        lat.col<-georedds.coords[2]
      }

      badrow<-which(is.na(georedds[,lon.col]) | is.na(georedds[,lat.col]))
      if(length(badrow)>0){
        georedds.sf<-sf::st_as_sf(georedds[-badrow,],coords=c(lon.col,lat.col),crs=reddcrs)
      }else{
        georedds.sf<-sf::st_as_sf(georedds,coords=c(lon.col,lat.col),crs=reddcrs)
      }
    }
    georedds.sf2<-sf::st_transform(georedds.sf,crs=sf::st_crs(shape))

    closest.reach<-shape.edges$Reach[sf::st_nearest_feature(x=georedds.sf2,shape.edges)]
    closest.dists<-sf::st_distance(georedds.sf2,shape.edges[match(closest.reach,shape.edges$Reach),],by_element = T)

    georedds.sf2$closest<-closest.reach
    georedds.sf3<-georedds.sf2[which(as.integer(closest.dists)<=tolerance.redds),]

    redd.day.col<-which(tolower(names(georedds.sf3))=="day")
    redd.year.col<-which(tolower(names(georedds.sf3))=="year")

    names(georedds.sf3)[redd.day.col]<-"Day"
    names(georedds.sf3)[redd.year.col]<-"Year"

    georedds.sf3$tempID<-NA
    # need to match redds to surveys
    pb<-utils::txtProgressBar(min = 0, max = nrow(georedds.sf3), initial = 0, style=3)
    for(i in 1:nrow(georedds.sf3)){
      good.surveys2<-good.surveys[surveys.dat[,year.col]==as.data.frame(georedds.sf3)[i,redd.year.col] &
                                    surveys.dat[,day.col]==as.data.frame(georedds.sf3)[i,redd.day.col],]
      if(nrow(good.surveys2)>0){
        survey.dists<-sf::st_length(sf::st_nearest_points(georedds.sf3[i,],good.surveys2))

        georedds.sf3$tempID[i]<-good.surveys2$tempID[which.min(survey.dists)]
      }
      utils::setTxtProgressBar(pb,i)
    }
    close(pb)
    sf::st_geometry(georedds.sf3)<-"geometry"
    georedds.sf3<-georedds.sf3[is.na(georedds.sf3$tempID)==F,c("Year","Day","closest","tempID")]
  }else{
    georedds.sf3<-NULL
  }

  # There are a lot of potential ways to assign non-georeferenced redds, but for now
  # we will just assume even spacing along the survey line
  # need to improve method of handling dates
  # Also need to adjust to ensure it functions when there are no survey redds
  print("Geo-locating Survey Redds w/o GPS")
  if(is.na(survey.redds)==F){
    has.redds<-which(as.data.frame(good.surveys)[,survey.redds]>0)

    pb<-utils::txtProgressBar(min = 0, max = length(has.redds), initial = 0, style=3)

    survey.redds.sf<-sf::st_as_sf(data.frame(Year=vector(),Day=vector(),closest=vector(),tempID=vector(),
                                             x=vector(),y=vector()),coords=c("x","y"),crs=sf::st_crs(shape))
    for(i in 1:length(has.redds)){
      active.survey<-good.surveys[has.redds[i],]
      redds<-as.data.frame(active.survey)[,survey.redds]

      if(inherits(georedds.sf3,"sf")){
        accounted.for<-georedds.sf3[georedds.sf3$tempID==active.survey$tempID,]
      }else{
        accounted.for<-data.frame(x=vector())
      }

      no.unaccounted<-redds-nrow(accounted.for)
      if(no.unaccounted>0){

        redds.sf<-sf::st_as_sf(sf::st_cast(sf::st_line_sample(active.survey,
                                                              sample = (1:no.unaccounted)/(no.unaccounted+1)),"POINT"))
        redds.sf$Year<-surveys.dat[has.redds[i],year.col]
        redds.sf$Day<-surveys.dat[has.redds[i],day.col]

        closest.reach<-shape.edges$Reach[sf::st_nearest_feature(x=redds.sf,shape.edges)]
        closest.dists<-sf::st_distance(redds.sf,shape.edges[match(closest.reach,
                                                                  shape.edges$Reach),],by_element = T)
        redds.sf$closest<-closest.reach
        redds.sf$tempID<-active.survey$tempID

        survey.redds.sf<-rbind(survey.redds.sf,redds.sf[as.numeric(closest.dists)<=tolerance.redds,])
      }
      utils::setTxtProgressBar(pb,i)
    }
    sf::st_geometry(survey.redds.sf)<-"geometry"
    close(pb)
  }else{
    survey.redds.sf<-NULL
  }
  good.redds<-rbind(survey.redds.sf,georedds.sf3)

  pb<-utils::txtProgressBar(min = 0, max = nrow(good.surveys), initial = 0, style=3)
  print("Resolving Survey Effort & assigning Redds")
  first.iter<-T
  for(i in 1:nrow(good.surveys)){

    overlapping<-which(sf::st_intersects(good.surveys[i,],shape.edges,sparse=F))

    if(length(overlapping)>0){
      for(j in 1:length(overlapping)){

        # because of imperfect snapping, we're going to do this the harder way
        # check which points are within the distance, and assume all inbetween points
        # ought to be included, even if the tracks are slightly off
        reach.shape<-shape.edges[overlapping[j],]
        reach.points<-suppressWarnings(sf::st_cast(reach.shape,"POINT"))

        good.points<-which(as.numeric(sf::st_length(sf::st_nearest_points(reach.points,good.surveys[i,])))<tolerance.surveys)

        # assume that any gaps are an error in snapping the shapes together
        # might cause problems if we ever start working with multi-channel streams, but ok for now
        line.points<-min(good.points):max(good.points)
        reach.intersection<-sf::st_cast(sf::st_combine(reach.points[line.points,]),"LINESTRING")
        intersection.length<-sf::st_length(reach.intersection)

        # drop very small intersections,
        if(as.numeric(intersection.length/sf::st_length(reach.shape))>.02){
          out.dat<-data.frame(Year=surveys.dat[i,year.col],
                              Day=surveys.dat[i,day.col],
                              Reach=overlapping[j],
                              Effort=intersection.length)

          active.redds<-good.redds[good.redds$closest==overlapping[j] &
                                     good.redds$tempID==good.surveys$tempID[i],]

          if(is.na(georedds.type)){
            Redd.counts<-data.frame(Redds=nrow(active.redds))
            redd.codes<-NA
          }else{
            redd.codes<-unique(as.data.frame(good.redds)[,georedds.type])
            Redd.counts<-as.data.frame(matrix(nrow=1,ncol=length(redd.codes),data=0))
            names(Redd.counts)<-paste("Redds",redd.codes,sep="_")
            for(t in 1:length(redd.codes)){
              Redd.counts[1,t]<-length(which(as.data.frame(active.redds)[,georedds.type]==redd.codes[t]))
            }
          }
          if(first.iter){
            out.data<-cbind(out.dat,Redd.counts)
            first.iter<-F
          }else{
            out.data<-rbind(out.data,cbind(out.dat,Redd.counts))
          }
        }
      }
    }
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)

  # Need to fix cases where the same reach is covered partially by two surveys
  # in future, may wish to preserve this difference, would require attaching a survey id to the data
  # Will need updating when it is time to work with surveys that don't have georeferenced redds
  doubles<-names(which(table(paste(out.data$Year,out.data$Day,out.data$Reach,sep="_"))>1))
  if(length(doubles)>0){
    for(i in 1:length(doubles)){
      y<-strsplit(doubles[i],split="_")[[1]][1]
      d<-strsplit(doubles[i],split="_")[[1]][2]
      r<-strsplit(doubles[i],split="_")[[1]][3]

      d.rows<-which(out.data$Year==y & out.data$Day==d & out.data$Reach==r)
      d.data<-out.data[d.rows,]

      new.data<-data.frame(Year=as.integer(y),Day=as.integer(d),Reach=as.integer(r),
                           Effort=sum(d.data$Effort))
      if(length(redd.codes)==1){new.redds<-data.frame(Redds=sum(d.data$Redds))}
      if(length(redd.codes)>1){new.redds<-t(as.data.frame(apply(d.data[,-(1:4)],MARGIN=2,FUN=sum)))}
      out.data<-rbind(out.data[-d.rows,],cbind(new.data,new.redds))
    }
  }
  # carry over any covariates
  var.col<-names(shape.edges)[names(shape.edges)%in%c("from","to","Reach","geometry")==F]
  reach.match<-match(out.data$Reach,shape.edges$Reach)
  out.data<-cbind(out.data,as.data.frame(shape.edges)[reach.match,var.col])
  rownames(out.data)<-NULL

  return(out.data)
}

#' Converts river measures and Lon/Lat coordinates
#'
#' This tool take a LINESTRING and either returns coordinates based on a set of river measures
#' or river measures based on a set of given coordinates. Note that a single, unbranched LINESTRING
#' should be used for shape, which may require calling the function multiple times for complex systems
#' The shape should be constructed in the right direction as well
#' (this should happen automatically if using the other functions).
#'
#' @param shape A sf LINESTRING corresponding stream
#' @param measures a vector of river measures, in the same units as shape
#' @param points a dataframe of coordinates or sf object with POINT geometries
#' @param crs A crs for the points
#'
#' @return a sf object with POINT geometries and matched coordinates and river measure data
#' @export
#'
#' @examples
RiverMeasureLL<-function(shape,measures,points,crs){

  root.sf<-LocateRoot(shape)

  if(missing(crs)){crs<-sf::st_crs(shape)}
  if(missing(measures)==F & missing(points)==F){stop("Doesn't make sense to use both points and river measure!")}
  if(is.na(root.sf[1])){stop("Root could not be located!")}

  if(missing(points)){
    riverline<-sf::st_line_merge(sf::st_combine(shape))
    if(sf::st_geometry_type(riverline)=="MULTILINESTRING"){
      river.coords<-sf::st_as_sf(as.data.frame(sf::st_coordinates(riverline)),coords=1:2,
                                 crs=sf::st_crs(shape))
      riverline<-sf::st_cast(sf::st_combine(river.coords),"LINESTRING")
    }
    points0<-sf::st_line_sample(riverline,sample = as.numeric(measures/sf::st_length(riverline)))
    points.sf2<-sf::st_as_sf(as.data.frame(sf::st_coordinates(points0))[,1:2],coords=1:2,
                             crs=sf::st_crs(shape))
    points.sf2$RiverDistance<-measures
  }
  if(missing(measures)){
    network<-suppressWarnings(sfnetworks::as_sfnetwork(shape))

    points.sf<-PointSetup(point=points,crs=sf::st_crs(crs))
    points.sf<-sf::st_transform(points.sf,sf::st_crs(shape))
    save.col<-which(names(points.sf)!=attr(points.sf,"sf_column"))

    network2<-suppressWarnings(sfnetworks::st_network_blend(network,sf::st_geometry(points.sf)))
    network2.nodes<-sf::st_as_sf(sfnetworks::activate(network2,"nodes"))

    points.sf2<-network2.nodes[sf::st_nearest_feature(points.sf,network2.nodes),]
    dists<-sfnetworks::st_network_cost(network2,from=root.sf,to = points.sf2)
    points.sf2<-cbind(points.sf2,as.data.frame(points.sf)[,save.col])
    points.sf2$RiverDistance<-as.numeric(sfnetworks::st_network_cost(network2,from=root.sf,to = points.sf2)[1,])
  }
  # Attach coordinates as data columns
  coordsLL<-sf::st_coordinates(sf::st_transform(points.sf2,crs = "wgs84"))
  coordsepgs<-sf::st_coordinates(sf::st_transform(points.sf2,crs = sf::st_crs(shape)))
  points.sf2<-cbind(points.sf2,data.frame(Lon=coordsLL[,1],Lat=coordsLL[,2],
                                          Eastings=coordsepgs[,1],Northings=coordsepgs[,2]))
  return(points.sf2)
}



#' Crop an sf shape - A graphical wrapper for the st_crop function
#'
#' @param shape A sf object
#' @param guides Any sf object you would like highlighted to help make your selections
#'
#' @return A sf object cropped to the user inputs
#' @export
#'
#' @examples
CropArea<-function(shape,guides,xmin,xmax,ymin,ymax){

  if(inherits(shape,"sfnetwork")){
    shape2<-sf::st_as_sf(sfnetworks::activate(shape,"edges"))
    makenet<-T
  }else{
    shape2<-shape
    makenet<-F
  }

  start<-as.numeric(sf::st_bbox(shape2))
  start2<-as.numeric(sf::st_bbox(sf::st_transform(shape2,crs="wgs84")))
  if(missing(guides)){guides<-list()}
  if(inherits(guides,"list")==F){stop("Guides argument must be a list")}

  if(any(c(missing(xmin),missing(xmax),missing(ymin),missing(ymax)))){
    map<-ggplot2::ggplot()+ggplot2::geom_sf(data=sf::st_geometry(shape2))
    if(length(guides)>0){
      for(i in 1:length(guides)){
        highlight<-sf::st_transform(guides[[i]],crs="wgs84")
        map<-map+ggplot2::geom_sf(data=highlight,col=(2:7)[i])
      }
    }
    map<-map+ggplot2::coord_sf(xlim=start[c(1,3)],ylim=start[c(2,4)])
    print(map)
  }
  print("Please input new coordinates for the bounding box; entering a blank will keep the current setting")
  if(missing(xmin)){xmin<-as.numeric(readline(prompt = "Minimum Longitude ="))}
  if(missing(xmax)){xmax<-as.numeric(readline(prompt = "Maximum Longitude ="))}
  if(missing(ymin)){ymin<-as.numeric(readline(prompt = "Minimum Latitude ="))}
  if(missing(ymax)){ymax<-as.numeric(readline(prompt = "Maximum Latitude ="))}

  if(is.na(xmin)){xmin<-start2[1]}
  if(is.na(xmax)){xmax<-start2[3]}
  if(is.na(ymin)){ymin<-start2[2]}
  if(is.na(ymax)){ymax<-start2[4]}

  bound.points0<-sf::st_as_sf(data.frame(x=c(xmin,xmax),y=c(ymin,ymax)),coords=1:2,crs="wgs84")
  bound.points<-sf::st_transform(bound.points0,crs=sf::st_crs(shape2))

  new.bounds<-as.numeric(sf::st_bbox(bound.points))

  new.shape<-sf::st_crop(shape2,bound.points)
  newmap<-ggplot2::ggplot()+ggplot2::geom_sf(data=new.shape)
  if(length(guides)>0){
    for(i in 1:length(guides)){
      highlight<-sf::st_transform(guides[[i]],crs="wgs84")
      newmap<-newmap+ggplot2::geom_sf(data=highlight,col=(2:7)[i])
    }
  }
  newmap<-newmap+ggplot2::coord_sf(xlim=new.bounds[c(1,3)],ylim=new.bounds[c(2,4)])
  print(newmap)

  # Sometimes this process breaks the geometries, so put in a step to try and fix
  broke<-which(sf::st_geometry_type(new.shape)!="LINESTRING")
  if(length(broke)>0){
    out.shape.ok<-new.shape[-broke,]
    fixed<-sf::st_cast(new.shape[broke,],"LINESTRING")
    out.shape<-rbind(out.shape.ok,fixed)
  }else{
    out.shape<-new.shape
  }

  if(makenet){
    out.shape<-sfnetworks::as_sfnetwork(out.shape)
  }

  return(out.shape)
}


