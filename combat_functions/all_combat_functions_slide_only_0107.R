### --- Script to adjust raw data with ComBat algorithm

### -------SETUP-------
require(dplyr)
require(tidyr)

## ---- COMBAT functions

## data adjustment function
remove_single_images = function(chan, image_var){
  ## count cells by images
  sub_chan = chan %>% group_by_at(image_var) %>% count()
  sub_chan$bool = sub_chan$n <= 1
  
  ## mark cells that are n-of-1 in an image
  nof1s = sub_chan[sub_chan$bool == TRUE,image_var]
  
  ## return dataset with no N-of-1s
  return(chan[!(chan$Pos %in% nof1s$Pos),])
}

## internal function for delta functions
sqerr = function(x){sum((x - mean(x))^2)}

## update each iteration of the algo
update_gamma = function(batch_chan, gamma_c, tau_c,channel, slide_var){
  ## create numerator value
  batch_chan$gamma_num = (batch_chan[,channel] - 
                            batch_chan$alpha_c)/batch_chan$delta_ijc
  gamma_num = batch_chan %>%
    group_by_at(slide_var) %>%
    summarise(avg = mean(gamma_num),.groups='drop')
  gamma_num$avg = gamma_num$avg + gamma_c/tau_c
  
  ## create denominator value
  gamma_denom = batch_chan %>%
    group_by_at(slide_var) %>%
    summarise(avg = mean(delta_ijc_inv),.groups='drop')
  gamma_denom$avg = gamma_denom$avg + (1/tau_c)
  
  gamma_ic_star = gamma_num
  gamma_ic_star$avg = gamma_ic_star$avg/gamma_denom$avg 
  
  ## returns zero if only one slide
  if(is.na(gamma_ic_star$avg[1])){gamma_ic_star$avg<-0}
  return(gamma_ic_star)
}
# update_lambda = function(batch_chan, lambda_c, eta_c,channel,image_var){
#   ##create numerator value
#   batch_chan$lambda_num = (batch_chan[,channel] - 
#                              batch_chan$alpha_c - 
#                              batch_chan$gamma_ic)/batch_chan$delta_ijc
#   lambda_num = batch_chan %>%
#     group_by_at(image_var) %>%
#     summarise(avg = mean(lambda_num),.groups='drop')
#   lambda_num$avg = lambda_num$avg + lambda_c/eta_c
#   
#   ## create denominator value
#   lambda_denom = batch_chan %>%
#     group_by_at(image_var) %>%
#     summarise(avg = mean(delta_ijc_inv),.groups='drop')
#   lambda_denom$avg = lambda_denom$avg + (1/eta_c)
#   
#   lambda_ijc_star = lambda_num
#   lambda_ijc_star$avg = lambda_ijc_star$avg/lambda_denom$avg
#   return(lambda_ijc_star)
# }
update_delta = function(batch_chan, beta_c,omega_c,channel,image_var){
  batch_chan$delta_num = beta_c + (batch_chan[,channel] - 
                                     batch_chan$alpha_c - 
                                     batch_chan$gamma_ic -
                                     batch_chan$lambda_ijc)^2
  delta_num = batch_chan %>%
    group_by_at(image_var) %>%
    summarise(avg = mean(delta_num),.groups='drop')
  
  delta_denom = batch_chan %>%
    group_by_at(image_var) %>%
    count()
  delta_denom$n = delta_denom$n/2 + omega_c - 1
  
  delta_ijc_star = delta_num
  delta_ijc_star$avg = delta_ijc_star$avg/delta_denom$n
  return(delta_ijc_star)
}
update_delta2 = function(batch_chan, beta_c,omega_c,channel,image_var){
  batch_chan$delta_num = (batch_chan[,channel] - 
                            batch_chan$alpha_c - 
                            batch_chan$gamma_ic)
  delta_num = batch_chan %>%
    group_by_at(image_var) %>%
    summarise(avg = sqerr(delta_num),.groups='drop')
  
  delta_denom = batch_chan %>%
    group_by_at(image_var) %>%
    count()
  
  delta_denom$n = delta_denom$n/2 + omega_c - 1
  delta_num$avg = 0.5*delta_num$avg + beta_c
  
  delta_ijc_star = delta_num
  delta_ijc_star$avg = delta_ijc_star$avg/delta_denom$n
  
  delta_ijc_star[is.na(delta_ijc_star$avg),]$avg = 0.00001
  
  return(delta_ijc_star)
}

## checking convergence
gamma_conv = function(batch_chan, gamma_stars,slide_var){
  gams = batch_chan[,c(slide_var,'gamma_ic')] %>% distinct()
  return(mean(abs(gams[match(unlist(gamma_stars[,slide_var]),
                             gams[,slide_var]),]$gamma_ic - gamma_stars$avg))) ## MAE
}
# lambda_conv = function(batch_chan, lambda_stars,image_var){
#   lambs = batch_chan[,c(image_var,'lambda_ijc')] %>% distinct()
#   return(mean(abs(lambs[match(unlist(lambda_stars[,image_var]),
#                               lambs[,image_var]),]$lambda_ijc - lambda_stars$avg))) ## MAE
# }
delta_conv = function(batch_chan, delta_stars,image_var){
  dels = batch_chan[,c(image_var,'delta_ijc')] %>% distinct()
  return(mean(abs(dels[match(unlist(delta_stars[,image_var]),
                             dels[,image_var]),]$delta_ijc - delta_stars$avg))) ## MAE
}


## function to combat-adjust for one channel
adjust_vals = function(channel,slide_var,image_var,uid_var,h,remove_zeroes=TRUE,
                       tol = 0.0001){
  print(channel)
  
  ### ---- Subset the data for the ComBat analysis
  chan = as.data.frame(h[,c(uid_var,slide_var,image_var,channel)])
  chan[,channel] = log10(chan[,channel]+1)
  
  if(remove_zeroes){
    ## remove zeroes if needed
    leftover = chan[chan[,channel] <=0,]
    chan = chan[(chan[,channel] > 0),]
    
    ## take ln
    chan[,channel] = log(chan[,channel])
  }
  
  
  ## fix n=1
  #chan = remove_single_images(chan, image_var)
  
  ### -------COMBAT EMPIRICAL VALUES-------
  
  ## get alpha (grand mean)
  chan$alpha_c = mean(chan[,channel])
  
  ## get gammas (slide means)
  gamma_ic = chan %>% 
    group_by_at(slide_var) %>% 
    summarise(avg=mean(get(channel)), .groups = 'drop')
  
  chan$gamma_ic = gamma_ic[match(chan[,slide_var],unlist(gamma_ic[,slide_var])),]$avg
  chan$gamma_ic = chan$gamma_ic - chan$alpha_c
  
  ## get lambdas (image means)
  # lambda_ijc = chan %>% 
  #   group_by_at(image_var) %>% 
  #   summarise(avg=mean(get(channel)), .groups = 'drop')
  # 
  # chan$lambda_ijc = lambda_ijc[match(chan[,image_var],unlist(lambda_ijc[,image_var])),]$avg
  # chan$lambda_ijc = chan$lambda_ijc - chan$alpha_c - chan$gamma_ic
  
  ## get deltas (image variances)
  #chan$delta_ijc = (chan[,channel] - chan$alpha_c - chan$gamma_ic - chan$lambda_ijc)^2
  chan$delta_ijc = (chan[,channel] - chan$alpha_c - chan$gamma_ic)
  ## delta_ijc = chan %>%
  ##   group_by_at(image_var) %>%
  ##   summarise(v=var(delta_ijc), .groups = 'drop')
  
  delta_ijc = chan %>%
    group_by_at(image_var) %>%
    summarise(v=sqerr(delta_ijc), .groups='drop')
  
  #delta_ijc[is.na(delta_ijc$v),]$v = 0.0001 ## there are images with only one value
  chan$delta_ijc = (delta_ijc[match(chan[,image_var],unlist(delta_ijc[,image_var])),]$v)
  
  ### -------COMBAT HYPERPARAMETERS-------
  ## slide level mean
  gamma_c = mean(chan$gamma_ic)
  tau_c = var(chan$gamma_ic)
  
  ## image level mean
  # lambda_c = mean(chan$lambda_ijc)
  # eta_c = var(chan$lambda_ijc)
  
  ## image level variances
  M_c = mean(chan$delta_ijc)
  S_c = var(chan$delta_ijc)
  
  ## is this correct?
  omega_c = (M_c + 2*S_c)/S_c
  beta_c = (M_c^3 + M_c*S_c)/S_c
  
  ### -------CALLING COMBAT BATCH EFFECTS FUNCTIONS-------
  batch_chan = chan ## duplicate the dataframe to iterate
  batch_chan$delta_ijc_inv = 1/batch_chan$delta_ijc
  
  # gamma_c; tau_c
  # lambda_c; eta_c
  # M_c; S_c
  # omega_c; beta_c
  
  ### -------COMBAT BATCH EFFECT ADJUSTMENT-------
  
  ## run a single iteration
  ## run delta first
  #delta_stars = update_delta(batch_chan, beta_c, omega_c,channel,image_var=image_var)
  delta_stars = update_delta2(batch_chan, beta_c, omega_c,channel,image_var)
  check_delta_conv = delta_conv(batch_chan, delta_stars,image_var)
  batch_chan$delta_ijc = (delta_stars[match(batch_chan[,image_var],unlist(delta_stars[,image_var])),]$avg)
  batch_chan$delta_ijc_inv = 1/batch_chan$delta_ijc
  
  ## now update gamma
  gamma_stars = update_gamma(batch_chan, gamma_c, tau_c,channel,slide_var=slide_var)
  check_gamma_conv = gamma_conv(batch_chan, gamma_stars,slide_var=slide_var)
  batch_chan$gamma_ic = gamma_stars[match(batch_chan[,slide_var],unlist(gamma_stars[,slide_var])),]$avg
  
  ## now update lambda
  # lambda_stars = update_lambda(batch_chan, lambda_c, eta_c,channel,image_var=image_var)
  # check_lambda_conv = lambda_conv(batch_chan, lambda_stars,image_var=image_var)
  # batch_chan$lambda_ijc = lambda_stars[match(batch_chan[,image_var],unlist(lambda_stars[,image_var])),]$avg
  
  total_mae = sum(check_gamma_conv,check_delta_conv)
  iterations = 1
  ## first check of MAE
  print(paste0('Total MAE after ', iterations,' iterations: ', round(total_mae,8)))
  
  ## run until convergence 
  while(total_mae > tol){ 
    ## run delta first
    #delta_stars = update_delta(batch_chan, beta_c, omega_c,channel,image_var=image_var)
    delta_stars = update_delta2(batch_chan, beta_c, omega_c,channel,image_var)
    check_delta_conv = delta_conv(batch_chan, delta_stars,image_var)
    batch_chan$delta_ijc = (delta_stars[match(batch_chan[,image_var],unlist(delta_stars[,image_var])),]$avg)
    batch_chan$delta_ijc_inv = 1/batch_chan$delta_ijc
    
    ## now update gamma
    gamma_stars = update_gamma(batch_chan, gamma_c, tau_c,channel,slide_var=slide_var)
    check_gamma_conv = gamma_conv(batch_chan, gamma_stars,slide_var=slide_var)
    batch_chan$gamma_ic = gamma_stars[match(batch_chan[,slide_var],unlist(gamma_stars[,slide_var])),]$avg
    
    ## now update lambda
    # lambda_stars = update_lambda(batch_chan, lambda_c, eta_c,channel,image_var=image_var)
    # check_lambda_conv = lambda_conv(batch_chan, lambda_stars,image_var=image_var)
    # batch_chan$lambda_ijc = lambda_stars[match(batch_chan[,image_var],unlist(lambda_stars[,image_var])),]$avg
    
    total_mae = sum(check_gamma_conv,check_delta_conv)
    iterations = iterations + 1
    ## final check of MAE
    print(paste0('Total MAE after ', iterations,' iterations: ', round(total_mae,4)))
  }
  
  ### -------COMBAT BATCH EFFECT RESULTS-------
  
  ## NOW actually adjust for the batch effects
  batch_chan$Y_ijc_star = (batch_chan[,channel] -
                             batch_chan$alpha_c -
                             batch_chan$gamma_ic)/batch_chan$delta_ijc

  
  
  ## exponential after natural log
  #batch_chan$Y_ijc_star = exp(batch_chan$Y_ijc_star)
  
  ## add zeroes back in if needed
  if(remove_zeroes){
    ## add back in zeroes
    leftover$Y_ijc_star = 0
    leftover[,colnames(batch_chan)[!(colnames(batch_chan) %in% colnames(leftover))]] = NA
    batch_chan = rbind(batch_chan,leftover)
  }
  
  return(batch_chan)
}


### NOTES ###
## dataset    | SARDANA | mouse
## -------    | ------- | -----
## slide_var  | SlideID | slideID
## image_var  | image   | view
## fov_var    | Pos     | Pos

run_full_combat = function(data,
                           save_path,
                           vars_to_adjust,
                           slide_var,
                           image_var,
                           uid_var,
                           remove_zeroes,
                           tol=0.001){
  ## create combat adjustment dir if necessary
  if(!dir.exists(paste0(save_path,'ComBat_adjustment_files/'))){
    dir.create(paste0(save_path,'ComBat_adjustment_files/'))
  }
  
  ## remove n-of-1s in the full dataset
  h_cb = data
  #h_cb = remove_single_images(h_cb,image_var)
  
  #alphas = c(); gammas = c(); deltas = c()
  
  ## adjust within each channel
  for (i in 1:length(vars_to_adjust)){
    ## adjust for the channel
    chan_i = adjust_vals(channel=vars_to_adjust[i],
                         slide_var = slide_var,
                         image_var = image_var,
                         uid_var = uid_var,
                         tol=tol,
                         h = data,
                         remove_zeroes = FALSE)
    
    ## save the dataframe
    saveRDS(chan_i, paste0(save_path,'ComBat_adjustment_files/',vars_to_adjust[i],'.rds'))
    
    ## replace the combined data with the combat-adjusted values
    h_cb[,paste0(vars_to_adjust[i],'_ComBat_Adjusted_slide_only_image_variance')] = chan_i$Y_ijc_star
    
    ## save channel vals
    #als = unique(chan_i$alpha_c)
    #alphas = c(alphas, als[!is.na(als)])
    
    #gams = unique(chan_i$gamma_ic)
    #gammas = c(gammas, gams[!is.na(gams)])
    
    #dels = unique(chan_i$delta_ijc)
    #deltas = c(deltas, dels[!is.na(dels)])
  }
  
  ## scale to give data natural scale
  # grand_mean = exp(mean(alphas))
  # grand_var = exp(mean(sqrt(deltas)))
  
  # for(v in paste0(vars_to_adjust,'_Adjusted')){
  #   h_cb[h_cb[,v] != 0,v] = (h_cb[h_cb[,v] != 0,v] + grand_mean)/grand_var
  # }
  
  ## save the adjusted dataset
  return(h_cb)
  
}

## ---- END FUNCTIONS
