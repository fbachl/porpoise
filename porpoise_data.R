setwd("/home/fbachl/devel/r/porpoise_scotland")

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
  
  data1 = read.csv("data/sightings_NEODAAS.csv")
  coords_data = as.matrix(data1[,c("Longitude","Latitude")])
  plot(data1[,c("Longitude","Latitude")])
  
  # data2<-read.csv("coords_sightings_no_out.csv")
  # points(data2[,c("Longitude","Latitude")], col = "red")
  
  coastline<-read.csv("data/coastline_updated.csv", header=TRUE)
  bnd.coast <- inla.mesh.segment(cbind(coastline$Long, coastline$Lat), is.bnd=1)
  
  #plot(coastline, type = "l")
  
  bnd1 <- inla.nonconvex.hull(coords_data, 0.05, resolution=100)
  bnd2 <- inla.nonconvex.hull(coords_data, 0.3)
  mesh <- inla.mesh.2d(boundary=list(bnd1, list(bnd2, bnd.coast)), max.edge=c(0.05, 0.5), cutoff=0.04)
  
  # Transform CRS
  mesh$loc[,c(1,2)] = coordinates(spTransform(SpatialPoints(mesh$loc[,c(1,2)], CRS("+proj=longlat")), CRS(target.p4s)))
  ggplot() + gg.mesh(mesh) + geom_point(data = data.frame(porpoise), aes(x=x, y=y, size=Size/eff))
  
  # Make new mesh since Laura's mesh does not cover all transects
  
      #' Boundary locations
      bnd.loc = porpoise$mesh$loc[mesh$segm$bnd$idx,c(1,2)]
      
      #' refined mesh
      rmesh = mesh.refine(mesh, refine = list(max.edge = 3))
      
      #' Sighting and effort locations
      smp.loc = do.call(rbind, lapply(coordinates(splines), function(x) x[[1]]))
      smp.loc = rbind(smp.loc, coordinates(pts))
      smp.hull = inla.nonconvex.hull(smp.loc, convex = -0.01)

      #' Vertices of new mesh
      rm.loc = rmesh$loc[,c(1,2)]
      bnd.new = inla.nonconvex.hull(rbind(rm.loc, smp.loc), convex = -0.007)
      mesh = inla.mesh.2d(boundary = bnd.new, max.edge = 5)
      # plot(mesh)
      # lines(porpoise$samplers)
      # all(is.inside(mesh, smp.loc))
      
      mesh$crs = CRS(target.p4s)
  
  ####### MAKE DATA SET ####
  
  porpoise = list(mesh = mesh, samplers = splines, points = pts, griddata = porpoise)
  save("porpoise", file = "porpoise.RData")
  
  