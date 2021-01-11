library(readxl)
library(magrittr)
library(ggplot2)
library(dplyr)

# pre_intervention enviroment
pre_iv = new.env()
# post_intervention enviroment
post_iv = new.env()

## load scan of pre and post intervention into two enviroments
load("data/cor_scan1.RData", envir = pre_iv)
load("data/cor_scan2.RData", envir = post_iv)

## load baseline data
baseline_1 = read_xlsx("data/PPA Baseline Database_MasterDatabase.xlsx",sheet = 1)
baseline_2 = read_xlsx("data/PPA Baseline Database_MasterDatabase.xlsx",sheet = 2)

## get subjects that have scaned both pre and post intervention 
pre_sub = pre_iv$subject[!sapply(pre_iv$tc.cor, is.null)]
post_sub = post_iv$subject[!sapply(post_iv$tc.cor, is.null)]

subjects = intersect(pre_sub, post_sub)


## pick a random subject and then print out the 
## density of correlations for pre and post intervention

sample.density = function(){
  subj = sample(subjects,1)
  n = dim(pre_iv$tc.cor[[subj]])[1]
  
  ## only upper triangular part of correlation
  corr_pre = pre_iv$tc.cor[[subj]][upper.tri(pre_iv$tc.cor[[subj]])]
  corr_post = post_iv$tc.cor[[subj]][upper.tri(post_iv$tc.cor[[subj]])]
  L = length(corr_pre)
  
  ## Trained items before and after phase 1
  tb1 = baseline_1$`Trained Items Before Phase 1`[baseline_1$Participant == subj]
  ta1 = baseline_1$`After 1..56`[baseline_1$Participant == subj]
  
  ## plot out densities
  data.frame(correlation = c(corr_pre, corr_post), 
             scan = c(rep('pre intervention',L), 
                      rep("post intervention",L))) %>%
    ggplot(aes(x = correlation, color = scan, fill = scan)) +
    geom_density(alpha = 0.2) + 
    ggtitle(paste("Density plot for correlatin coefficients\nSubject: ", subj, 
                  "\t\t Before:", round(tb1,2), "\t After:", round(ta1,2))) +
    theme_bw()
}


####
ROI = c('IFG_opercularis_L', 'IFG_opercularis_R',
        'IFG_orbitalis_L', 'IFG_orbitalis_R',
        'IFG_triangularis_L', 'IFG_triangularis_R',
        'SMG_L', 'FuG_L', 'STG_L', 'STG_L_pole',
        'MTG_L', 'MTG_L_pole', 'ITG_L')
DMN = c('MFG_DPFC_L', 'MFG_DPFC_R', 'AG_L', 'AG_R', 'PCC_L', 'PCC_R')
RES = setdiff(colnames(post_iv$tc.cor[[1]]), c(ROI,DMN))

##### baseline variables

features = baseline_1 %>% rename(Language_Severity = `Language Severity`, Age = `Age at Start of Therapy`,
                                 First_condition = `First-period condition`,
                                 Trained_items_1 = `Trained Items Before Phase 1`,
                                 Untrained_items_1 = `Untrained Items Before Phase 1`,
                                 Trained_After_1 = `After 1..56`,
                                 Untrained_After_1 = `After 1..64`) %>%
  select(Participant, First_condition, Variant, Gender, Language_Severity, Age, 
         Trained_items_1, Untrained_items_1, Trained_After_1, Untrained_After_1) %>%
  mutate(Participant = toupper(Participant),
         trained_diff = Trained_After_1-Trained_items_1, 
         trained_rdiff = (Trained_After_1-Trained_items_1)/(100-Trained_items_1)*100,
         untrained_diff = Untrained_After_1-Untrained_items_1, 
         untrained_rdiff = (Untrained_After_1-Untrained_items_1)/(100-Untrained_items_1)*100)

features$First_condition[startsWith(features$First_condition, "s")] = "s"
features$First_condition[startsWith(features$First_condition, "t")] = "t"
features$Variant[startsWith(features$Variant, "n")] = "n"

### post intervention data
data_post = merge(data.frame(Participant = post_sub), features, by = 'Participant')
data_post$Participant = as.character(data_post$Participant)

# data_post = data_post %>% dplyr::filter(Participant %in% pre_sub)

corr_post = c(); corr_z_post = c(); post_flag = c()
for(subj in data_post$Participant){
  mi = post_iv$tc.cor[[subj]]
  corr_i = mi[upper.tri(mi)]
  corr_i = corr_i[!is.na(corr_i)]
  post_flag = c(post_flag, length(corr_i) == 3003)
  if(length(corr_i) < 3003){
    corr_i = c(corr_i, sample(corr_i,3003-length(corr_i)))
  }
  zi = post_iv$tc.cor.z[[subj]]
  corr_z_i = zi[upper.tri(zi)]
  corr_z_i = corr_z_i[!is.na(corr_z_i)]
  if(length(corr_z_i) < 3003){
    corr_z_i = c(corr_z_i, sample(corr_z_i,3003-length(corr_z_i)))
  }
  corr_post = rbind(corr_post, corr_i)
  corr_z_post = rbind(corr_z_post, corr_z_i)
}
rownames(corr_post) = data_post$Participant
rownames(corr_z_post) = data_post$Participant

Apost = as.numeric(data_post$First_condition == 't')
Xpost = model.matrix(~Variant+Gender+Age-1, data_post)
Xpost = Xpost[,2:5]

###
rm(baseline_1,baseline_2)

### bootstrap

bootstrap.sample = function(Y,X1,X2,X3,sub=0,split=0){
  ind = sample(length(Y), length(Y), replace=TRUE)
  ind2 = sample(length(Y), length(Y), replace=TRUE)
  if(split > 0){
    ind = sample(length(Y), floor(split*length(Y)), replace=TRUE)
    ind2 = -ind
  }
  Y_train = Y[ind]; X1_train = X1[ind,]
  X2_train = X2[ind,]; X3_train = X3[ind,]
  Y_val = Y[ind2]; X1_val = X1[ind2,]
  X2_val = X2[ind2,]; X3_val = X3[ind2,]
  if(sub > 0){
    for(i in 1:length(Y_train)){
      X2_train[i,] = sample(X2_train[i,], dim(X2)[2], replace=TRUE)
      X3_train[i,] = sample(X3_train[i,], dim(X3)[2], replace=TRUE)
    }
    for(i in 1:length(Y_val)){
      X2_val[i,] = sample(X2_val[i,], dim(X2)[2], replace=TRUE)
      X3_val[i,] = sample(X3_val[i,], dim(X3)[2], replace=TRUE)
    }
    X2_train = X2_train[,1:sub]; X3_train = X3_train[,1:sub]
    X2_val = X2_val[,1:sub]; X3_val = X3_val[,1:sub]
  }
  list(train = list(Y=Y_train,X1=X1_train,X2=X2_train,X3=X3_train),
       val = list(Y=Y_val,X1=X1_val,X2=X2_val,X3=X3_val))
}

### plot functions
library(fda)

plot.f = function(xs, ys, nbasis=21, xlab='correlation', ylab='signal'){
  basisobj = create.bspline.basis(c(min(xs),max(xs)), nbasis)
  if(is.null(dim(ys))){
    xs = matrix(xs, 1, length(xs))
    ys = matrix(ys, 1, length(ys))
  }
  else if(is.null(dim(xs))){
    xs = matrix(rep(xs,dim(ys)[1]),dim(ys)[1],dim(ys)[2],byrow = TRUE)
  }
  para = c()
  for(i in 1:dim(ys)[1]){
    para = cbind(para,smooth.basis(xs[i,],ys[i,],basisobj)$fd$coef)
  }
  plot(fd(para, basisobj),xlab=xlab,ylab=ylab)
}

plot.ds = function(corr_trans, np=256, nb=21, qd=FALSE){
  xs = c(); ys = c()
  for(i in 1:nrow(corr_trans)){
    ds_i = density(corr_trans[i,], kernel="gaussian", 
                   from=min(corr_trans) * 1.01, 
                   to=max(corr_trans) * 1.01, n=np)
    if(qd){
      ds_i$y = exp(-dens2lqd(ds_i$y, ds_i$x))
      ds_i$x = seq(0,1,length.out=np)
    }
    xs = rbind(xs, ds_i$x)
    ys = rbind(ys, ds_i$y)
  }
  plot.f(xs, ys, nb)
  list(xs=xs,ys=ys)
}

standardize = function(X){
  Xnew = matrix(NA, nrow(X), ncol(X))
  for(i in 1:nrow(X)){
    Xnew[i,] = (X[i,]-mean(X[i,])) / sqrt(var(X[i,]))
  }
  Xnew
}

rnorm.mul = function(mui, Sig){
  ei = eigen(Sig)
  C = ei$vectors %*% diag(sqrt(ei$values)) %*% t(ei$vectors)
  c(C %*% rnorm(length(mui)) + mui)
}

nb = 19
basis1 = create.bspline.basis(c(-1,1), nb)
basis2 = create.bspline.basis(c(0,1), nb)
T0 = function(corr_i){
  # di = density(corr_i, kernel="gaussian", from=-1, to=1, n=512)
  densf = density.fd(corr_i, fd(matrix(0,nb,1), basis1))
  di = list(x=seq(-1,1,length.out=512),
            y=c(exp(eval.fd(seq(-1,1,length.out=512),densf$Wfdobj))/densf$C))
  list(x=di$x, y=di$y)
}
Tl = function(corr_i){
  # di = density(corr_i, kernel="gaussian", from=-1, to=1, n=512)
  densf = density.fd(corr_i, fd(matrix(0,nb,1), basis1))
  di = list(x=seq(-1,1,length.out=512),
            y=c(eval.fd(seq(-1,1,length.out=512),densf$Wfdobj)-log(densf$C)))
  # list(x=di$x, y=log(1e-7+di$y))
  list(x=di$x, y=di$y)
}
Tstd = function(corr_i){
  corr_std = (corr_i - mean(corr_i)) / sqrt(var(corr_i)) 
  di = density(corr_std, kernel="gaussian", from=-4, to=4, n=512)
  list(x=di$x, y=di$y)
}
Tq = function(corr_i){
  x = seq(0,1,length.out=512)
  y = quantile(corr_i, x)
  list(x=x, y=y)
}
Tqd = function(corr_i){
  di = density(corr_i, kernel="gaussian", from=min(corr_i)-0.1, to=1, n=512)
  y = exp(-dens2lqd(di$y, di$x))
  x = seq(0,1,length.out=512)
  list(x=x, y=y)
}
Tlqd = function(corr_i){
  di = density(corr_i, kernel="gaussian", from=min(corr_i)-0.1, to=1, n=512)
  # densf = density.fd(corr_i, fd(matrix(0,nb,1), basis1))
  # di = list(x=seq(-1,1,length.out=1024),
  #           y=c(exp(eval.fd(seq(-1,1,length.out=1024),densf$Wfdobj))/densf$C))
  # ind1 = min(which(di$y > 1e-3)); ind2 = max(which(di$y > 1e-3))
  # di = list(x = di$x[ind1:ind2], y = di$y[ind1:ind2])
  tryCatch({
     y = -dens2lqd(di$y, di$x, N=512)
     x = seq(0,1,length.out=512)
     list(x=x, y=y)
    },
    error = function(e){
      di$y = di$y + 1e-7
      y = -dens2lqd(di$y, di$x, N=512)
      x = seq(0,1,length.out=512)
      list(x=x, y=y)
    }
  )
}

FC = ecdf(corr_post)
Tgq = function(corr_i){
  corr_i = FC(corr_i)
  T0(corr_i)
}

## generate T for running model
T.generate = function(Tc, Ds, Ps){
  function(corr_i){
    corr_i = Tc(corr_i)
    di = Ds(corr_i)
    Ps(di$x, di$y)
  }
}

## show pairs selected
show.pairs = function(net, fnames){
  net[upper.tri(net)] = FALSE
  outer(fnames, fnames, function(a,b){paste(a,"<-->",b)})[net]
}

a = 3
a = tryCatch({44}, error = function(e){a = 4})
