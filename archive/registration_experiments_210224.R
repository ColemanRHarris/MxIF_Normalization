require(fda)
var = 'Median_Cell_VIMENTIN_log10'
cb_atl = readRDS("~/../../media/disk2/atlas_mxif/combat/handling_zeroes/method9_0122/atl_with_all_methods_updated_0122.rds")
len = 512

## ----- starts here
x = split(cb_atl,factor(cb_atl$SlideID))
rang = range(cb_atl[ which(cb_atl[,var]!=0),var])

## setup basis functions
fdobj_basis = create.bspline.basis(rangeval = rang,
                                   norder = 4)
wbasis = create.bspline.basis(rangeval = rang, norder = 2)
Wfd0   <- fd(matrix(0,wbasis$nbasis,1),wbasis)
WfdPar <- fdPar(Wfd0,
                Lfdobj = int2Lfd(0))

#indices = round(seq(1,length(vals),length.out = 1024))
argvals = seq(rang[1],rang[2],len=1024)

## density
densY = sapply(1:length(x), function(i){
  vals = x[[i]][,var][which(x[[i]][,var]>0)]
  z = ecdf(vals)
  z(argvals)
})

# test_dat = data.frame(cbind(unlist(densYp),rep(argvals,length(x))))
# test_dat$SlideID = rep(levels(factor(cb_atl$SlideID)), each=1024)
# ggplot(test_dat) + geom_point(aes(x=X2,y=X1,color=SlideID),alpha=0.25) + theme(legend.position = "None")

## capture the shape of the histogram
fdobj   <- smooth.basis(argvals, densY, fdobj_basis, ## estimates the density using bsplines
                           fdnames = c("x", "samples", "density"))$fd

fdobj = deriv.fd(fdobj)

regDens = register.fd(yfd=fdobj, 
                      WfdParobj = WfdPar,
                      dbglev = 0)

#dfdobj <- deriv.fd(fdobj)

y0s = mean.fd(fdobj)
y0s$coefs = do.call(cbind,
                    rep(list(y0s$coefs),
                        ncol(fdobj$coefs)))
regDens = register.fd(y0fd = fdobj, 
                      yfd=y0s, 
                      WfdParobj = WfdPar,
                      dbglev = 0)
# 
# 
# slides = levels(factor(test_dat$SlideID))
# x1 = lapply(1:length(slides),
#             function(i){
#               test_dat[test_dat$SlideID == slides[i],'X2p'] = eval.fd(test_dat[test_dat$SlideID == slides[i],'X2'], regDens$warpfd[i])
#               test_dat[test_dat$SlideID == slides[i],]
#             })
# 
# pdat1 = data.table::rbindlist(x1)
# ggplot(pdat1) + geom_point(aes(x=X2p,y=X1,color=SlideID),alpha=0.25) + theme(legend.position = "None")
# 
# s1 = pdat1[pdat1$SlideID == unique(pdat1$SlideID)[1],]
# #ggplot(s1) + geom_point(aes(x=X2p,y=X1,color=SlideID),alpha=0.25) + theme(legend.position = "None")
# 
# s2 = pdat1[pdat1$SlideID == unique(pdat1$SlideID)[2],]
# ggplot(s2) + geom_point(aes(x=X2p,y=X1,color=SlideID),alpha=0.25) + theme(legend.position = "None")
# 
# finite.differences <- function(x, y) {
#   if (length(x) != length(y)) {
#     stop('x and y vectors must have equal length')
#   }
#   
#   n <- length(x)
#   
#   # Initialize a vector of length n to enter the derivative approximations
#   fdx <- vector(length = n)
#   
#   # Iterate through the values using the forward differencing method
#   for (i in 2:n) {
#     fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
#   }
#   
#   # For the last value, since we are unable to perform the forward differencing method 
#   # as only the first n values are known, we use the backward differencing approach
#   # instead. Note this will essentially give the same value as the last iteration 
#   # in the forward differencing method, but it is used as an approximation as we 
#   # don't have any more information
#   fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
#   
#   return(fdx)
# }
# 
# fin1 = finite.differences(s1$X2p, s1$X1)
# fin2 = finite.differences(s2$X2p, s2$X1)
# 
# original_form = lapply(X=1:length(slides),
#        FUN=function(i){
#          yp = finite.differences(pdat1[pdat1$SlideID == unique(pdat1$SlideID)[i],]$X2p,
#                             pdat1[pdat1$SlideID == unique(pdat1$SlideID)[i],]$X1)
#          df = data.frame(cbind(yp,
#                                pdat1[pdat1$SlideID == unique(pdat1$SlideID)[i],]$X2p))
#          df$SlideID = slides[i]
#          df
#        })
# 
# plot(s1$X2p, fin1,col="red")
# points(s2$X2p, fin2,col="blue")

y0s = mean.fd(dfdobj)
y0s$coefs = do.call(cbind,
                    rep(list(y0s$coefs),
                        ncol(dfdobj$coefs)))
regDens = register.fd(y0fd = dfdobj, 
                      yfd=y0s, 
                      WfdParobj = WfdPar,
                      dbglev = 0)

#deriv.fd(fdobj)

## transform back into data
normedVar = paste0(var, '_fda_registered')

x = lapply(1:length(x), 
           function(ind){ 
             x[[ind]][,normedVar] = 0;
             x[[ind]][ which(x[[ind]][,var]>0),normedVar] = eval.fd(x[[ind]][which(x[[ind]][,var]>0), var], regDens$warpfd[ind]);
             x[[ind]]})

pdat = data.table::rbindlist(x)
cb_atl[,normedVar] = data.frame(pdat)[,normedVar]
rm(pdat)

## ------ ends here
## unadjusted densities

ggplot(cb_atl) + geom_density(aes(x=Median_Cell_VIMENTIN_log10,col=SlideID)) + theme(legend.position = 'None')
ggplot(pdat) + geom_density(aes(x=Median_Cell_VIMENTIN_log10_fda_registered,col=SlideID)) + theme(legend.position = 'None')
