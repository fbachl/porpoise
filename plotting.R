
laura.plot = function(mesh, values, ng = 100, ...) {
  

  pgrid <- inla.mesh.projector(mesh, dims = c(ng,ng))
  pred<- inla.mesh.project(pgrid, values)
  
  data0 <- matrix(NA, ncol=3, nrow=length(pgrid$x)*length(pgrid$y))
  data0[,1]<-rep(pgrid$y,each=length(pgrid$x))
  data0[,2]<-rep(pgrid$x,length(pgrid$y))
  data0[,3]<-rapply(list(pred), c)
  data0<-data.frame(data0)
  names(data0)<-c("Lat","Lon","LatentFieldMean")
  
  library(fields)
  quilt.plot(data0$Lon,data0$Lat,data0$LatentFieldMean, nx = ng, ny = ng, ...)  # this works
  
}



