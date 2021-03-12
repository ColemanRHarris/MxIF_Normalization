require(mclust)

## load in data
#atl <- readRDS("~/../../media/disk2/atlas_mxif/colon_map_20201209.rds")
#cb_atl = readRDS('~/../../media/disk2/atlas_mxif/combat_201214/colon_map_20201209_combat_adj.rds')

## function to correctly align classification
get_class = function(dat,mod){
  sub_dat = dat[mod$classification == 1,]
  mv = mean(sub_dat$Median_Cell_VIMENTIN)
  mo = mean(unlist(sub_dat[,c('Median_Cell_NAKATPASE','Median_Cell_PANCK')]))
  
  if(mv > mo){
    return(ifelse(mod$classification==2,1,2))
  }
  return(mod$classification)
}

## function to fit clustering algos 
get_clusters = function(cb_atl, 
                        clus_vars,
                        clus_type){
  c1 = paste0('cluster_within_slide_',clus_type)
  c2 = paste0('cluster_across_slide_',clus_type)
  
  cb_atl[,c1] = -1.0
  cb_atl[,c2] = -1.0
  
  for(s in 1:length(all_slides)){
    print(all_slides[s])
    ssub = cb_atl[cb_atl$SlideID == all_slides[s],]
    smod = Mclust(as.data.frame(ssub)[,clus_vars], 
                  G=2, 
                  prior=priorControl(), 
                  modelNames='EEI')
    
    cb_atl[cb_atl$SlideID == all_slides[s],c1] = get_class(ssub,smod)
  }
  
  
  ## calculate clusters across slides
  fullmod = Mclust(as.data.frame(cb_atl)[,clus_vars], 
                   G=2, 
                   prior=priorControl(), 
                   modelNames='EEI')
  
  cb_atl[,c2] = get_class(cb_atl, fullmod)
  return(cb_atl)
}

clus_vars =  c('Median_Cell_NAKATPASE',
               'Median_Cell_PANCK',
               'Median_Cell_VIMENTIN')
               #,'Median_Cell_COLLAGEN')
all_clus = list(clus_vars,
                paste0(clus_vars,"_log10"),
                paste0(clus_vars,"_simple_adjusted"),
                paste0(clus_vars,"_combat_slide_adjusted")
                )
all_types = c("raw","log","simple","cb")

for(i in 1:length(all_clus)){
  cb_atl = get_clusters(cb_atl,
                        clus_vars=all_clus[[i]],
                        clus_type=all_types[i])
}

for(c in colnames(cb_atl)[grepl("cluster",colnames(cb_atl))]){
  cb_atl[,c] = as.factor(as.data.frame(cb_atl)[,c])
}

## just one adjustment
all_slides = unique(cb_atl$SlideID)

## VIM: 3.554297
## PANCK: 3.586874
## NAKATPASE:  3.468282

## cb_atl[cb_atl$SlideID == "MAP02112_0000_02_02",]$Median_Cell_VIMENTIN_log10_simon = cb_atl[cb_atl$SlideID == "MAP02112_0000_02_02",]$Median_Cell_VIMENTIN_log10_simon - 3.554297
## cb_atl[cb_atl$SlideID == "MAP02112_0000_02_02",]$Median_Cell_PANCK_log10_simon = cb_atl[cb_atl$SlideID == "MAP02112_0000_02_02",]$Median_Cell_PANCK_log10_simon - 3.586874
## cb_atl[cb_atl$SlideID == "MAP02112_0000_02_02",]$Median_Cell_NAKATPASE_log10_simon = cb_atl[cb_atl$SlideID == "MAP02112_0000_02_02",]$Median_Cell_NAKATPASE_log10_simon - 3.468282

require(dbscan)

#x = hdbscan(cb_atl[,paste0(clus_vars,'_log10')],minPts = 50)
cb_atl$Median_Cell_COLLAGEN_log10 = log10(cb_atl$Median_Cell_COLLAGEN+1)

cb_atl = get_clusters(cb_atl,
                      clus_vars = paste0(clus_vars,"_log10_fda_registered10"),
                      clus_type = "fda_registered1")

cb_atl$cluster_across_slide_fda_registered1 = as.factor(cb_atl$cluster_across_slide_fda_registered1)
cb_atl$cluster_within_slide_fda_registered1 = as.factor(cb_atl$cluster_within_slide_fda_registered1)


confusionMatrix(data = cb_atl$cluster_across_slide_log,
                reference = cb_atl$cluster_within_slide_log)

cmat = confusionMatrix(data = cb_atl$cluster_across_slide_fda_registered1,
                reference = cb_atl$cluster_within_slide_log)

## use k-means to compare clusters across methods
## bronze standard: k-means clusters fit within each slide on log10
## simple comparison: k-means clusters fit across slides on log10
## experimental: k-means clusters fit across slides on fda registered

## first pass: univariate clustering with otsu
## after normalization/registration, otsu thresholds will be closer than unnormalized data
## given threshold, compare confusion matrices
## otsu(as.array(image), range = c(0,255)), range is range of the data

## require(EBImage)
## require(autothresholdr)
## kmeans()


## kmeans analysis



X = cb_atl[,paste0(clus_vars,"_log10")]

l = list()
for(s in 1:length(all_slides)){
  print(s)
  l[[s]] = confusionMatrix(data = cb_atl[cb_atl$SlideID == all_slides[s],]$cluster_across_slide_fda_registered1,
                  reference = cb_atl[cb_atl$SlideID == all_slides[s],]$cluster_within_slide_log)$overall['Accuracy']
}

l1 = list()
for(s in 1:length(all_slides)){
  print(s)
  l1[[s]] = confusionMatrix(data = cb_atl[cb_atl$SlideID == all_slides[s],]$cluster_across_slide_log,
                           reference = cb_atl[cb_atl$SlideID == all_slides[s],]$cluster_within_slide_log)$overall['Accuracy']
}

xs = c(1:32)
y_fdas = unlist(l)
y_logs = unlist(l1)

plot(xs,y_fdas,col='red')
points(xs,y_logs, col='blue')

ggpairs(cb_atl[cb_atl$SlideID == "MAP01391_0000_01_01",],
        columns=paste0(clus_vars,"_log10"),
        lower=list(continuous = wrap("points",alpha=0.25)),
        ggplot2::aes(color=cluster_across_slide_log)) +
        theme_minimal()



saveRDS(cb_atl,
        "~/../../media/disk2/atlas_mxif/combat/handling_zeroes/method9_0122/atl_with_all_methods_updated_0122.rds")

### ------- deprecated vvvvvvvv -------

## separate data by slide
all_slides = unique(cb_atl$SlideID)
clus_vars = c('Median_Cell_NAKATPASE',
              'Median_Cell_PANCK',
              'Median_Cell_VIMENTIN')
clus_vars_cb = c('Median_Cell_NAKATPASE_Adjusted',
              'Median_Cell_PANCK_Adjusted',
              'Median_Cell_VIMENTIN_Adjusted')

## calculate clusters within slide
cb_atl$cluster_within_slide_raw = -1.0
cb_atl$cluster_within_slide_cb = -1.0

get_class = function(dat,mod){
  sub_dat = dat[mod$classification == 1,]
  mv = mean(sub_dat$Median_Cell_VIMENTIN)
  mo = mean(unlist(sub_dat[,c('Median_Cell_NAKATPASE','Median_Cell_PANCK')]))

  if(mv > mo){
    return(ifelse(mod$classification==2,1,2))
  }
  return(mod$classification)
}

for(s in 1:length(all_slides)){
  print(all_slides[s])
  ssub = cb_atl[cb_atl$SlideID == all_slides[s],]
  smod = Mclust(ssub[,clus_vars], 
                G=2, 
                prior=priorControl(), 
                modelNames='EEI')
  
  smod_cb = Mclust(ssub[,clus_vars_cb], 
                   G=2, 
                   prior=priorControl(), 
                   modelNames='EEI')
  
  ## c1 = ssub[(smod$classification == 1), ]
  ## c1cb = ssub[(smod_cb$classification == 1), ]
  ## 
  ## ## c1 :: VIM downregulated
  ## m1 = mean(c1$Median_Cell_VIMENTIN)
  ## m2 = mean(unlist(c1[,c('Median_Cell_NAKATPASE','Median_Cell_PANCK')]))
  ## if(m1 > m2){
  ##   print("cluster 1 seems weird")
  ##   print(paste0("Mean vim: ",m1))
  ##   print(paste0("Mean nakatpase,panck: ",m2))
  ## }
  ## 
  ## m1cb = mean(c1cb$Median_Cell_VIMENTIN)
  ## m2cb = mean(unlist(c1cb[,c('Median_Cell_NAKATPASE','Median_Cell_PANCK')]))
  ## if(m1cb > m2cb){
  ##   print("cluster 1 cb seems weird")
  ##   print(paste0("Mean vim cb: ",m1cb))
  ##   print(paste0("Mean nakatpase,panck cb: ",m2cb))
  ## }
  
  cb_atl[cb_atl$SlideID == all_slides[s],]$cluster_within_slide_raw = get_class(ssub,smod)
  cb_atl[cb_atl$SlideID == all_slides[s],]$cluster_within_slide_cb = get_class(ssub,smod_cb)
}

## calculate clusters across slides
fullmod = Mclust(cb_atl[,clus_vars], 
                 G=2, 
                 prior=priorControl(), 
                 modelNames='EEI')

fullmod_cb = Mclust(cb_atl[,clus_vars_cb], 
                 G=2, 
                 prior=priorControl(), 
                 modelNames='EEI')

cb_atl$cluster_across_slide_raw = get_class(cb_atl, fullmod)
cb_atl$cluster_across_slide_cb = get_class(cb_atl, fullmod_cb)

## load data
saveRDS(cb_atl,
        '~/../../media/disk2/atlas_mxif/combat/handling_zeroes/method3/colon_map_20201217_combat_adj_with_clusters.rds')
#saveRDS(cb_atl,'~/../../media/disk2/atlas_mxif/combat_201214/colon_map_20201214_combat_adj_with_clusters.rds')
#saveRDS(cb_atl,'~/../../media/disk2/atlas_mxif/combat_201214/colon_map_20201214_combat_adj_with_corrected_clusters.rds')
#cb_atl = readRDS('~/../../media/disk2/atlas_mxif/combat_201214/colon_map_20201214_combat_adj_with_clusters.rds')

## plotting
cb_atl$cluster_within_slide_raw = as.factor(cb_atl$cluster_within_slide_raw)
cb_atl$cluster_within_slide_cb = as.factor(cb_atl$cluster_within_slide_cb)
cb_atl$cluster_across_slide_raw = as.factor(cb_atl$cluster_across_slide_raw)
cb_atl$cluster_across_slide_cb = as.factor(cb_atl$cluster_across_slide_cb)

# plot_slide = "MAP01391_0000_01_01"
# plot_atl = cb_atl[cb_atl$SlideID == plot_slide,]
# plot_atl_raw = plot_atl[,c( clus_vars,
#                             clus_vars_cb,
#                            'cluster_within_slide_raw',
#                            'cluster_across_slide_raw',
#                            'cluster_within_slide_cb',
#                            'cluster_across_slide_cb')]
# 
# require(GGally)
# 
# 
# ggpairs(plot_atl_raw,
#         columns=1:3,
#         ggplot2::aes(color=cluster_across_slide_cb))
