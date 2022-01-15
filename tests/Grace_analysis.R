library(ngBap)
library(multiedge)
plot.dat <- read.csv("data/grace_plot_level.csv")
site.dat <- read.csv("data/grace_site_level.csv")



plot.sem.dat <- with(plot.dat, data.frame(psitecode))
names(plot.sem.dat)[names(plot.sem.dat)=="psitecode"] <- "PlotSiteCode"
plot.sem.dat$PlotRichRaw    <- plot.dat$prich
plot.sem.dat$PlotRich       <- plot.dat$ln.prich
plot.sem.dat$SiteRich       <- plot.dat$ln.site.rich
plot.sem.dat$PlotProdRaw    <- plot.dat$pprod
plot.sem.dat$PlotProd       <- plot.dat$ln.pprod
plot.sem.dat$SiteProd       <- plot.dat$ln.site.prod
plot.sem.dat$PlotBiomass    <- plot.dat$ln.ptotmass
plot.sem.dat$SiteBiomass    <- plot.dat$ln.site.totmass
plot.sem.dat$PlotShade      <- plot.dat$ln.pshade
plot.sem.dat$PlotSoilSuit   <- plot.dat$SoilSuitability
plot.sem.dat$SoilWithShade  <- plot.dat$SoilWithShade #control variable for soil influences on effectiveness of shading
plot.sem.dat$PlotN          <- plot.dat$pn
plot.sem.dat$PlotN.sqr      <- plot.dat$sqr.pn
plot.sem.dat$PlotN.ln       <- plot.dat$ln.pn
plot.sem.dat$PlotP.ln       <- plot.dat$ln.pp
plot.sem.dat$PlotC.ln       <- plot.dat$ln.pc

Y <- cbind(plot.sem.dat[,c("PlotRich", "PlotSoilSuit", "SiteRich", "PlotShade", "PlotBiomass",
                           "PlotProd","SiteBiomass", "SiteProd")])#, SiteSoilSuit = site.dat$SoilSuitability[match(plot.sem.dat$PlotSiteCode, site.dat$site.code)])


Y <- scale(Y)
out_01 <- bang(Y, 4, level = .01, verbose = T)
out_05 <- bang(Y, 4, level = .05, verbose = T)
out_001 <- bang(Y, 4, level = .001, verbose = T)


out_01a <- bang(Y, 4, level = .01, verbose = T, restrict  = 2)
out_05a <- bang(Y, 4, level = .05, verbose = T, restrict = 2)
out_001a <- bang(Y, 4, level = .001, verbose = T, restrict = 2)



out_01_M=MBANG(t(Y),out_01)
out_05_M=MBANG(t(Y),out_05)
out_001_M=MBANG(t(Y),out_001)
out_01a_M=MBANG(t(Y),out_01a)
out_05a_M=MBANG(t(Y),out_05a)
out_001a_M=MBANG(t(Y),out_001a)

print(out_01_M)
#[[1]]
#[1] 2 6 7

#[[2]]
#[1] 2 7 8

#[[3]]
#[1] 3 7 8
