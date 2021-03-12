cb_atl = readRDS("~/../../media/disk2/atlas_mxif/combat/handling_zeroes/method7_0107/atl_with_all_methods_updated_0107.rds")

## LIMMA
## BiocManager::install("limma")
# require(limma)
# y = cb_atl$Median_Cell_VIMENTIN_log10
# batch = cb_atl$SlideID
# y2 = removeBatchEffect(y, batch)
# 
# design <- model.matrix(~ -1 + cb_atl$SlideID)
# fit <- lmFit(cb_atl[,paste0("Median_Cell_",c('VIMENTIN'),'_log10')], design)
# fit <- eBayes(fit)
# topTable(fit)

## SVA
# require(sva)
# mod = model.matrix(~as.factor(SlideID),data=cb_atl)
# mod0 = model.matrix(~1,data=cb_atl)
# svobj =sva(cb_atl[,c("Median_Cell_VIMENTIN_log10","SlideID")],mod,mod0)

## Python
## https://github.com/ushaham/BatchEffectRemoval2018

## MATLAB
## https://gitlab.com/Chang_Lab/cycif_int_norm

## R
## https://github.com/julia-wrobel/mica

require(mica)
require(ggplot2)
require(ggridges)

mica_dat = cb_atl[,c('SlideID','Pos','Median_Cell_VIMENTIN_log10')]
mica_dat$channel = "Median_Cell_VIMENTIN_log10"
colnames(mica_dat)[3] = "intensity"

#mica_dat$image = paste0(mica_dat$SlideID,"+",mica_dat$Pos)
#mica_dat$mica_id = paste(c("MedianCellVIMENTINlog10","Pos","SlideID"),collapse = "_")
#mica_dat_cdf = calculate_cdf(mica_dat)
#hist(mica_dat_cdf$cdf,breaks=1000)

## log10
ggplot(mica_dat) +
  geom_density_ridges(aes(x=intensity,y=SlideID,fill=SlideID)) +
  theme_minimal() +
  theme(legend.position = "None") +
  ggtitle("log10 intensity of VIM across slides")

## template for the channel

## warping function for each image?
mica::inverse_warps()


make_intensity_df = function(channels, subj_scan_scanner, ...) 
{
  intensities = map2_dfr(channels, subj_scan_scanner, vectorize_channel)
  as_tibble(intensities)
}

vectorize_channel <- function(channel, subj_scan_scanner, ...){
  df = data.frame(intensity = mica_dat[mica_dat$channel == channel,]$intensity)
  
  df$subject = strsplit(subj_scan_scanner, "_")[[1]][1]
  df$scan_id = strsplit(subj_scan_scanner, "_")[[1]][2]
  df$scanner = strsplit(subj_scan_scanner, "_")[[1]][3]
  df$voxel_position = row.names(df)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  mode = getmode(df$intensity)
  df = dplyr::mutate(df, mode = mode)
  df = dplyr::filter(df, intensity != mode)
  
  df
}

channels = c("Median_Cell_VIMENTIN_log10")
subj_scan_scanner = "channel_Pos_SlideID"

intensity_df = make_intensity_df(channels = c("Median_Cell_VIMENTIN_log10"), subj_scan_scanner = "channel_Pos_SlideID")

estimate_cdf <- function(intensity_df,
                         grid_length = 1000, ...){
  
  intensity_df = tidyr::nest(intensity_df, data = c(intensity, voxel_position))
  intensity_df = mutate(intensity_df, data = map(data, calculate_cdf))
  
  # downsample cdf to smaller, regular grid
  cdf_mat = intensity_mat = matrix(0, nrow = grid_length, ncol = dim(intensity_df)[1])
  
  for(i in 1:dim(intensity_df)[1]){
    
    intensity_minimum = intensity_df$intensity_minimum[i]
    intensity_maximum = intensity_df$intensity_maximum[i]
    intensity_grid = seq(intensity_minimum, intensity_maximum, length.out = grid_length)
    
    cdf_mat[, i] = approx(intensity_df$data[[i]]$intensity,
                          intensity_df$data[[i]]$cdf,
                          xout = intensity_grid, rule = 2, ties = mean)$y
    
    intensity_mat[, i] = intensity_grid
  }
  ls = list(intensity_df = intensity_df, cdf_mat = cdf_mat, intensity_mat = intensity_mat)
  
  return(ls)
}

## https://bioconductor.org/books/release/OSCA/messmer-hesc.html#batch-correction-1
## https://bioconductor.org/books/release/OSCA/merged-pancreas.html
## https://bioconductor.org/books/release/OSCA/lun-416b-cell-line-smart-seq2.html#batch-correction
## https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#1_Introduction
## https://bioconductor.org/packages/release/bioc/vignettes/pvca/inst/doc/pvca.pdf
## https://bioconductor.org/packages/3.12/bioc/vignettes/batchelor/inst/doc/correction.html

