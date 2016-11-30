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

  rp = bru(griddata, predictor = Size ~ spde + Intercept, mesh = porpoise$mesh, E = griddata$eff, linear = TRUE)

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
  
  mdl1 = ~ spat(model=inla.spde2.pcmatern(porpoise$mesh, prior.range=c(20, 0.1), prior.sigma=c(100, 0.1)), mesh = porpoise$mesh) +
           Intercept - 1
  
#' Model 2: Default SPDE parameterization
  
  mdl2 = ~ spat(model=inla.spde2.matern(porpoise$mesh), mesh = porpoise$mesh) +
           Intercept - 1

#' Run inference
#' 
#+warning=FALSE,results='hide'
  
  r1 = bru(griddata, model = mdl1, predictor = Size ~ spat + Intercept, E = griddata$eff, linear = TRUE, family = "nbinomial")
  r2 = bru(griddata, model = mdl2, predictor = Size ~ spat + Intercept, E = griddata$eff, linear = TRUE, family = "nbinomial")

  
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

  mdl = ~ spat(model=inla.spde2.matern(porpoise$mesh), mesh = porpoise$mesh) + DepthCentered + Intercept - 1
  rc = bru(griddata, model = mdl, predictor = Size ~ spat + Intercept + DepthCentered, E = griddata$eff, linear = TRUE, family = "nbinomial")
  
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

r.p.pois = bru(porpoise$meshdata, predictor = count ~ spde + Intercept, mesh = porpoise$mesh, 
               E = porpoise$meshdata$area, linear = TRUE, family = "poisson")
r.p.nbin = bru(porpoise$meshdata, predictor = count ~ spde + Intercept, mesh = porpoise$mesh, 
               E = porpoise$meshdata$area, linear = TRUE, family = "nbinomial")


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
  
  
#'### Covariates
#'
#' For this example I will use the covariates as stored in the `griddata` object. Ideally, I would
#' have the original data and apply the following procedure. Let's start with interpolating
#' the SST data.

  SSTdata = covariate(griddata, predictor = SST, mesh = porpoise$mesh)
  
#+out.width='50%'
  pl1 = ggplot() + gg(porpoise$mesh) + gg(griddata, mapping = aes(x,y,color=SST), size = 4)
  pl2 = plot(SSTdata) + gg(porpoise$mesh) + gg(porpoise$points) + gg(porpoise$samplers)
  pl1 ; pl2
  
#' The interpolant can now be used to create functions that return the SST value at any point
#' in space.

  local.sst = evaluator(SSTdata)
  mean.sst = evaluator(SSTdata, get.smean) # Returns spatial mean
  centered.sst = function(...) local.sst(...) - mean.sst(...)
  # For example
  centered.sst(x=-150, y=6000)
  
#' Unfortunately there is currently a problem with using covariates within the predictor. It works but
#' it is extremely slow. A workaround is to evaluate the covariates at the detection and integration points
#' beforehand. For that purpose we need a covariate function that acts on data frames:
   
  centered.sst2 = function(df) centered.sst(x=df$x, y=df$y)
  
#' Now, let's build a model!
  
  mdl = ~ sst.eff + 
    g(spat, model = inla.spde2.matern(porpoise$mesh), mesh = porpoise$mesh) + 
    Intercept - 1
  
#' ... and a predictor. Note how `csst` is used instead of the more straight forward `centered.sst(x,y)` syntax.  
  prd = coordinates ~ sst.eff * csst + spat + Intercept -1
  
#' Via the `append` command we let `lgcp` know that it should evaluate centered.sst2 for all points and
#' append the respective values to the data frame used during inference. This makes `csst` a valid
#' expression in the predictor. 

  r = lgcp(points = porpoise$points, 
           model = mdl, 
           predictor = prd, 
           mesh = porpoise$mesh, 
           n = 1,
           append = list(csst = centered.sst2))

#' Let's see if the effect is significant: Nope!
  plot.marginal(r, "sst.eff")
  r$summary.fixed
  
  
#' Predictions work as usual. This will predict the joint intensity:
  pr = predict(r, coordinates ~ exp(sst.eff * centered.sst(x,y) + spat + Intercept))
  plot(pr)

#' And this is only SPDE plus intercept  
  pr2 = predict(r, coordinates ~ exp(spat + Intercept))
  plot(pr2)

  