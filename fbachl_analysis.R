#'---
#' title: "Harbour porpoise Scotland"
#' author: Fabian E. Bachl
#' output:
#'   html_document:
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'       smooth_scroll: false
#'---

#' Install and load inlabru
#'=================================================

#' Use this to install inlabru:
#+eval=FALSE
library(devtools)
install_git("https://github.com/fbachl/inlabru.git")

#' Assuming that inlabru is installed:
library(inlabru)

#' Regression analysis 
#'=================================================

#' Data plots
#'----------------------------

#' Load data

  setwd("/home/fbachl/devel/r/porpoise_scotland/")
  load("porpoise.RData")
  griddata = porpoise$griddata
  
#' Plot the data (currently deactived as my rgdal package is not functioning)
#' 
#+eval=FALSE

  library(ggmap)
  gg.map(griddata) + gg.point(griddata, CRS = CRS("+proj=longlat"))

#' Non-discretized data
  ggplot() + gg.mesh(porpoise$mesh) + gg.segment(porpoise$samplers) +gg.point(porpoise$points)
  
#' Detections per unit effort and effort per grid cell. It is easy to see that there 
#' are some cells with an extremely large number of animals per effort.    
#' 
#+out.width='50%'
  
  ggplot() + gg.mesh(porpoise$mesh) + 
    geom_point(data = data.frame(griddata), aes(x=x, y=y, color=Size/eff), size = 3) +
    scale_color_gradientn(colors = topo.colors(100)) ; ggplot() + gg.mesh(porpoise$mesh) + 
    geom_point(data = data.frame(griddata), aes(x=x, y=y, color=eff), size = 5) +
    scale_color_gradientn(colors = topo.colors(100))


#' Poisson Likelihood
#'----------------------------

#' Run inferece using default model (rhs of formula), lhs determines count and exposure columns in data frame
#' 
#+warning=FALSE,results='hide'

  rp = poiss(griddata, Size + eff ~ ., mesh = porpoise$mesh)

#'### Results

#' Predict spatial intensity pattern. The default model has SPDE plus Intercept, which we combine
#' and then apply $exp$ as inverse link. Also predict log intensity.

  int = predict(rp, coordinates ~ exp(Intercept + spde))
  log.int = predict(rp, coordinates ~ Intercept + spde)

#' At a first glance the result is not very informative. The intensity field seems flat with some very local
#' peaks. In the log scale it becomes clear that the small local areas of high density  exactly reflect 
#' the grid cells with high animal counts per effort area.
#' 
#+out.width='50%'

  plot(int); plot(log.int) + geom_point(data = data.frame(griddata), aes(x=x, y=y, color=Size/eff), size = 1) +
    scale_color_gradientn(colors = topo.colors(100))

#' Summary: When using the Poisson likelihood the SPDE will pick up the local animal clusters. Hence,
#' we obtain a low range and a field that is not very smooth on the larger spatial scale. 
#' 

  
#' Negative Binomial Likelihood
#'--------------------------------

#' The negative Binomial likelihood can mitigate the local clustering effect by allowing some additional
#' degree of overdispersion. In what follows we will compare the effect of using PCpriors.  
#'
#' Model 1: PCprior
  
  fml1 = Size + eff ~ g(spat, model=inla.spde2.pcmatern(porpoise$mesh, prior.range=c(20, 0.1), prior.sigma=c(100, 0.1)), mesh = porpoise$mesh)
  
#' Model 2: Default SPDE parameterization
  
  fml2 = Size + eff ~ g(spat, model=inla.spde2.matern(porpoise$mesh), mesh = porpoise$mesh)

#' Run inference
#' 
#+warning=FALSE,results='hide'
  
  r1 = poiss(griddata, fml1, family = "nbinomial")
  r2 = poiss(griddata, fml2, family = "nbinomial")

  
#'### Results

#' Hyper parameter posteriors

  sr1 = inla.spde.result(r1, "spat", inla.spde2.pcmatern(porpoise$mesh, prior.range=c(20, 0.1), prior.sigma=c(100, 0.1)))  
  sr2 = inla.spde.result(r2, "spat", inla.spde2.matern(porpoise$mesh))

#' Not much difference!

  exp(rbind(sr1$summary.log.range.nominal, sr1$summary.log.variance.nominal))
  exp(rbind(sr2$summary.log.range.nominal, sr2$summary.log.variance.nominal))
  
      
#' Spatial intensity predictions
#' 
#+warning=FALSE

  int1 = predict(r1, coordinates ~ exp(Intercept + spat))
  int2 = predict(r2, coordinates ~ exp(Intercept + spat))

#' Plot them side by side and on the same scale
#' 
#+out.width='50%',warning=FALSE,message=FALSE,results='hide'

  plot(int1) + scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"), limits = c(0,0.6)); plot(int2) + scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"), limits = c(0,0.6))
  
#' Plot using Laura's method (don't evaluate by default) 
#+eval=FALSE
  
  source("fbachl_plotting.R")
  laura.plot(porpoise$mesh, exp(r1$summary.random$spat$mean + r1$summary.fixed["Intercept","mean"])) ; 
  laura.plot(porpoise$mesh, exp(r2$summary.random$spat$mean + r2$summary.fixed["Intercept","mean"]))
  
#' Intensity integrated over space (just as an example)
#' 
#+warning=FALSE
  
  tint = predict(r1, ~ exp(Intercept + spat), integrate = "coordinates")
  plot(tint)

  
  
#' Covariates
#'--------------------------------
#'
#' This is a simple example on how to use covariates. Let's have a look at the effect of depth. 
#' We need to remove the mean of depth to avoid confounding the effect of depth with the intercept

  griddata$DepthCentered = griddata$Depth - mean(griddata$Depth)

#' Set up the model and run inference
#'   
#+warning=FALSE,results='hide'

  fml = Size + eff ~ g(spat, model=inla.spde2.matern(porpoise$mesh), mesh = porpoise$mesh) + DepthCentered
  rc = poiss(griddata, fml, family = "nbinomial")
  
#' Use `plot.marginal` to inspect the posterior of the depth effect. The Predictive interval covers 0,
#' so the effect is not significant:

  plot.marginal(rc, "DepthCentered")

#' Let's predict the combined intensity over space. Nothing fancy going on here.
  
  pr = predict(rc, coordinates ~ exp(DepthCentered + spat + Intercept))
  plot(pr)

#' What if we only predict the exponential depth effect over space? The values are around 1
#' everywhere. This makes sense as we figured out that the depth effect is not significant.

  pr = predict(rc, coordinates ~ exp(DepthCentered))
  plot(pr)
  