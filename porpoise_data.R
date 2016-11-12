setwd("/home/fbachl/devel/git/porpoise")

####### SETTINGS ########################

  # target.p4s = "+proj=laea +init=epsg:3035 +units=km"
  target.p4s = "+proj=laea +units=km"

####### SURVEY DATA #######################
  
  porpoise = read.csv("data/sightings_NEODAAS.csv")
  
  ## Transform data equal area projection
  
  coordinates(porpoise) =  ~ Longitude + Latitude
  proj4string(porpoise) <- CRS("+proj=longlat")
  porpoise = spTransform(porpoise, CRS(target.p4s))
  coordnames(porpoise) = c("x","y")
  
  # Set effort
  porpoise$eff = porpoise$Effort_Tot*(porpoise$Strip_Widt/1000)
  

  
  
  ##### EXACT DATA #########################
  
  lns = read.csv("data/Survey_Lines_XY_Points.csv")
  pts = read.csv("data/All_Sightings.csv")
  
  # POINTS
  coordinates(pts) = c("Longitude","Latitude")
  proj4string(pts) = "+proj=longlat"
  coordnames(pts) = c("lon","lat")
  
  
  # LINES
  nlines = nrow(lns)
  sp = lns[,c("StartX","StartY")]
  ep = lns[,c("EndX","EndY")]
  colnames(sp) = c("lon","lat")
  colnames(ep) = colnames(sp)
  lilist = lapply(1:nrow(sp), function(k) { Lines(list(Line(rbind(sp[k,],ep[k,]))), ID = k) } )
  splines = SpatialLines(lilist, proj4string = CRS(proj4string(pts)))
  splines = SpatialLinesDataFrame(SpatialLines(lilist, proj4string = CRS(proj4string(pts))), 
                         data = data.frame(weight = lns$Strip_Width/1000))
  
  pts = spTransform(pts, CRS(target.p4s))
  splines = spTransform(splines, CRS(target.p4s))
  
  
  ####### MAKE MESH #######################

  coastline<-read.csv("data/Coastline_updated_for_INLAbru.csv", header=TRUE)
  coordinates(coastline) = c("Long","Lat")
  proj4string(coastline) = "+proj=longlat"
  coastline = spTransform(coastline, CRS(target.p4s))
  plot(coastline)
  bnd.coast <- inla.mesh.segment(coordinates(coastline), is.bnd=1)
  
  tr.co = do.call(rbind, lapply(coordinates(splines), function(s) s[[1]]))
  coords_data = rbind(tr.co, coordinates(pts))
  
  bnd1 <- inla.nonconvex.hull(coords_data, 9)
  bnd2 <- inla.nonconvex.hull(coords_data, 27)
  mesh <- inla.mesh.2d(boundary=list(bnd1, list(bnd2, bnd.coast)), min.angle = 25, max.edge=c(4, 10), cutoff=3.8)
  plot(mesh)
  mesh$n
  points(coastline)
  points(coords_data)
  lines(splines)

  mesh$crs = CRS(target.p4s)
  
  ###### SSAVE EVERYTHING ################
  porpoise = list(mesh = mesh, samplers = splines, points = pts, griddata = porpoise)
  save("porpoise", file = "porpoise.RData")
  
  