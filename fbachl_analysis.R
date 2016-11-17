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

#' more stuff
library(rgdal)
library(ggmap)

#' The data 
#'=================================================

#' Load data

  setwd("/home/fbachl/devel/git/porpoise/")
  load("porpoise.RData")
  griddata = porpoise$griddata
  
#' Survey data 
#'--------------------

ovg = ggplot() + gg(porpoise$mesh) + gg(porpoise$samplers) + gg(porpoise$points)
ovm = gg.map(porpoise$points, zoom=7) + gm(porpoise$mesh) + gm(porpoise$samplers) + gm(porpoise$points)

#+out.width='50%'

ovg ; ovm


#' Grid aggregation
#'--------------------

  
  gd.cnt = ggplot() + gg.mesh(porpoise$mesh) + 
    geom_point(data = data.frame(griddata), aes(x=x, y=y, color=Size), size = 3) +
    scale_color_gradientn(colors = topo.colors(100), limits = c(0,20)) ; 
  
  gd.eff = ggplot() + gg.mesh(porpoise$mesh) + 
    geom_point(data = data.frame(griddata), aes(x=x, y=y, color=eff), size = 5) +
    scale_color_gradientn(colors = topo.colors(100))
  
  gd.dns = ggplot() + gg.mesh(porpoise$mesh) + 
    geom_point(data = data.frame(griddata), aes(x=x, y=y, color=Size/eff), size = 3) +
    scale_color_gradientn(colors = topo.colors(100)) ; 
  
  gd.ldns = ggplot() + gg.mesh(porpoise$mesh) + 
    geom_point(data = data.frame(griddata), aes(x=x, y=y, color=log(Size/eff)), size = 3) +
    scale_color_gradientn(colors = topo.colors(100)) ; 

#' Animal count (left) and effort (right) per grid cell
#+out.width='50%'
  gd.cnt ; gd.eff 

#' Density (count/effort) and lod density per cell
#+out.width='50%'
  gd.dns ; gd.ldns


#' Mesh aggregation 
#'--------------------


  md.cnt = ggplot() + gg(porpoise$mesh) + 
    geom_point(data.frame(porpoise$meshdata), mapping = aes(x,y,color=count), size = 3) + 
    scale_color_gradientn(colors = topo.colors(100), limits = c(0,20))
  md.eff = ggplot() + gg(porpoise$mesh) + 
    geom_point(data.frame(porpoise$meshdata), mapping = aes(x,y,color=area), size = 3) + 
    scale_color_gradientn(colors = topo.colors(100))
  md.dns = ggplot() + gg(porpoise$mesh) + 
    geom_point(data.frame(porpoise$meshdata), mapping = aes(x,y,color=count/area), size = 3) + 
    scale_color_gradientn(colors = topo.colors(100))
  md.ldns = ggplot() + gg(porpoise$mesh) + 
    geom_point(data.frame(porpoise$meshdata), mapping = aes(x,y,color=log(count/area)), size = 3) +
    scale_color_gradientn(colors = topo.colors(100))

#+out.width='50%'
md.cnt ; md.eff
#+out.width='50%'
md.dns ;md.ldns
  
#' Comparison
#'--------------------

#+out.width='50%'
gd.cnt ; md.cnt

#+out.width='50%'
gd.eff ; md.eff  



#' Grid regression 
#'=================================================

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
  
  source("plotting.R")
  laura.plot(porpoise$mesh, exp(r1$summary.random$spat$mean + r1$summary.fixed["Intercept","mean"])) ; 
  laura.plot(porpoise$mesh, exp(r2$summary.random$spat$mean + r2$summary.fixed["Intercept","mean"]))
  
#' Intensity integrated over space (just as an example)
#' 
#+warning=FALSE
  
  tint = predict(r1, ~ exp(Intercept + spat), integrate = "coordinates")
  plot(tint)


#' Comparison
#'--------------------------------
  int.g.pois = plot(int) + scale_fill_gradientn(colours = brewer.pal(7,"YlOrRd"),  limits = c(0,0.5))
  int.g.nbin1 = plot(int1) + scale_fill_gradientn(colours = brewer.pal(7,"YlOrRd"),  limits = c(0,0.5))
  int.g.nbin2 = plot(int2) + scale_fill_gradientn(colours = brewer.pal(7,"YlOrRd"),  limits = c(0,0.5))

#+out.width='30%'  
  int.g.pois ; int.g.nbin1 ; int.g.nbin2

  #+out.width='30%'  
  int.g.pois ; int.g.nbin1 ; int.g.nbin2
  
  
  
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

#' In order to make predictions we need to be able to evaluate the covariate everywhere.
#' For that purpose we use 'covariate' to interpolate the data. 'evaluator' will then 
#' give you a function that can be applied to coordinates x,y and will return the depth
#' at that position

  depth.cv = covariate(griddata, predictor = DepthCentered, mesh = porpoise$mesh)
  depth.ev = evaluator(depth.cv)
  
#' Predict depth over space (actuall we are not predicting depth here, this is just to check
#' if the covariate does the right thing)

  pr = predict(rc, coordinates ~ depth.ev(x,y), n = 2)
  plot(pr)

#' Predict the effect of the covariate. Note that the naming is a bit confusing. depth.ev is
#' the actual depth at some location. DepthCentered is the name of the depth effect, so just
#' some factor we are looking for.

  pr = predict(rc, coordinates ~ DepthCentered * depth.ev(x,y), n = 250)
  plot(pr)

#' Let's predict the spatial intensity with and without the covariate. As the depth effect
#' is not significant it is not surprising that the patterns are basically identical.

  pr1 = predict(rc, coordinates ~ exp(Intercept + spat + DepthCentered * depth.ev(x,y)))
  pr2 = predict(rc, coordinates ~ exp(Intercept + spat))
  
#+out.width='50%'    
  
  plot(pr1) ; plot(pr2)
  
    
#' Mesh regression 
#'=================================================
#'
#' 
#' `porpoise$meshdata` is a projection of the exact surey data to mesh vertices. This
#' is the triangular equivalent to a grid-projection but makes use of the assumption that
#' the underlying intensity is linear within each triangle.
#' 
#' Poisson regression on this data should give results that are similar to LGCP modeling.
#' However, instead of Poisson regression we can also perform neg-Binomial regression which
#' we can compare to the respective result using gridded data.
#' 
#' The statistical interpretation of this is not clear!!!!!
#' 

#' Let's run Poisson and neg-Binomial regression.

r.p.pois = poiss(porpoise$meshdata, model = count + area ~ ., mesh = porpoise$mesh, family = "poisson")
r.p.nbin = poiss(porpoise$meshdata, model = count + area ~ ., mesh = porpoise$mesh, family = "nbinomial")


int.p.pois = predict(r.p.pois, coordinates ~ exp(Intercept + spde))
int.p.nbin = predict(r.p.nbin, coordinates ~ exp(Intercept + spde))

#' Comparison
#'* Poisson on mesh projection
#'* nBin on mesh projection
#'* nBin on grid
#'
#+out.width='30%'

p1 = plot(int.p.pois) + scale_fill_gradientn(colours = brewer.pal(7,"YlOrRd"),  limits = c(0,2.2))
p2 = plot(int.p.nbin) + scale_fill_gradientn(colours = brewer.pal(7,"YlOrRd"),  limits = c(0,2.2))
p3 = plot(int2) + scale_fill_gradientn(colours = brewer.pal(7,"YlOrRd"),  limits = c(0,2.2))

p1 ; p2 ; p3

#' Point process analysis 
#'=================================================
#'
#' Run LGCP inference
  
  r = lgcp(porpoise$points, samplers = porpoise$samplers, mesh = porpoise$mesh)
  
#' Plot integration weights
  
  ggplot() + geom_point(data = data.frame(r$ips), aes(x=X1, y=X2, color = weight) , size = 3)
  
#' Predict spatial (log) intensity
  
  spint = predict(r, coordinates ~ exp(spde + Intercept))
  lspint = predict(r, coordinates ~ spde + Intercept)
  
#+out.width='50%'
  
  plot(spint) + gg(porpoise$points, size = 1) ; plot(lspint) + gg(porpoise$points, size = 1)
  
