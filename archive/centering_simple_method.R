cb_atl = readRDS("~/../../media/disk2/atlas_mxif/combat/handling_zeroes/method5_1229/atl_with_clusters1.rds")

require(ggplot2)
require(ggridges)

ggplot(cb_atl1) +
  #geom_density(aes(Median_Cell_VIMENTIN,color="raw")) +
  geom_density(aes(Median_Cell_VIMENTIN_log10,color="log10",group=SlideID)) +
  geom_density(aes(Median_Cell_VIMENTIN_simple_adjusted,color="simple",group=SlideID)) +
  geom_density(aes(Median_Cell_VIMENTIN_simple_adjusted_centered,color="centered",group=SlideID)) +
  #geom_density(aes((log(Median_Cell_VIMENTIN_simple_adjusted+1) + 2.506974)/1.01929,color="simple",group=SlideID)) +
  #geom_density(aes(Median_Cell_VIMENTIN_combat_slide_adjusted,color="combat")) +
  scale_x_continuous(limits = c(0,10)) +
  scale_y_continuous(limits = c(0,2)) +
  theme_minimal()

## calc avg slide mean
slide_mean = mean(aggregate(cb_atl$Median_Cell_VIMENTIN_log10, list(cb_atl$SlideID), mean)$x)

## calc avg slide variance
slide_var = mean(aggregate(cb_atl$Median_Cell_VIMENTIN_log10, list(cb_atl$SlideID), var)$x)


cb_atl$Median_Cell_VIMENTIN_simple_adjusted_centered = (cb_atl$Median_Cell_VIMENTIN_simple_adjusted + 
  mean(aggregate(cb_atl$Median_Cell_VIMENTIN_log10, list(cb_atl$SlideID), mean)$x))

cb_atl$Median_Cell_PANCK_simple_adjusted_centered = (cb_atl$Median_Cell_PANCK_simple_adjusted + 
                                                       mean(aggregate(cb_atl$Median_Cell_PANCK_log10, list(cb_atl$SlideID), mean)$x))

cb_atl$Median_Cell_NAKATPASE_simple_adjusted_centered = (cb_atl$Median_Cell_NAKATPASE_simple_adjusted + 
                                                 mean(aggregate(cb_atl$Median_Cell_NAKATPASE_log10, list(cb_atl$SlideID), mean)$x))



cb_atl1 = readRDS("~/../../media/disk2/atlas_mxif/combat/handling_zeroes/method5_1229/atl_with_clusters_and_centered_1231.rds")

ggplot(cb_atl1) +
  #geom_density(aes(Median_Cell_VIMENTIN,color="raw")) +
  geom_density(aes(Median_Cell_VIMENTIN_log10,color="log10",group=SlideID)) +
  geom_density(aes(Median_Cell_VIMENTIN_simple_adjusted,color="simple",group=SlideID)) +
  geom_density(aes(Median_Cell_VIMENTIN_simple_adjusted_centered,color="centered",group=SlideID)) +
  scale_x_continuous(limits = c(0,10)) +
  scale_y_continuous(limits = c(0,2)) +
  theme_minimal()

ggplot(cb_atl1) +
  geom_density_ridges(aes(Median_Cell_VIMENTIN_log10,SlideID,fill="log10")) +
  scale_x_continuous(limits = c(0,5))

ggplot(cb_atl1) +
  geom_density_ridges(aes(Median_Cell_VIMENTIN_simple_adjusted,SlideID,fill="simple adjusted")) +
  scale_x_continuous(limits = c(0,5))

ggplot(cb_atl1) +
  geom_density_ridges(aes(Median_Cell_VIMENTIN_combat_slide_adjusted,SlideID,fill="combat")) +
  scale_x_continuous(limits = c(0,5))
