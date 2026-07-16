# =============================================================================
# Assumption Tests
# =============================================================================

# Linear correlation assumption tests - Linearity, Homoscedasticity, Norm of residuals
library(crqa)
library(R.matlab)
library(lmtest)
setwd('')
data=readMat('Experimental_R.mat')

# Windowed versions of:
# - Pearson correlation p-value
# - Breusch–Pagan p-value (Homoscedasticity)
# - Shapiro p-value (normality of residuals)
# and returns:
# (a) per-window results (one row per dyad × window)
# (b) aggregated violation rates

library(lmtest)

windowed_tests <- function(x, y, window_size = 50L, overlap = 0L) {
  stopifnot(is.numeric(x), is.numeric(y))
  stopifnot(length(x) == length(y))
  window_size <- as.integer(window_size)
  overlap <- as.integer(overlap)
  if (window_size < 3L) stop("window_size must be >= 3 (needed for lm + tests).")
  if (overlap < 0L || overlap >= window_size) stop("overlap must be in [0, window_size-1].")
  
  n <- length(x)
  step <- window_size - overlap
  starts <- seq.int(1L, n - window_size + 1L, by = step)
  ends <- starts + window_size - 1L
  
  out <- vector("list", length(starts))
  for (w in seq_along(starts)) {
    idx <- starts[w]:ends[w]
    xi <- x[idx]
    yi <- y[idx]
    
    # Guard against constant series inside a window (cor / lm can misbehave)
    if (stats::sd(xi) == 0 || stats::sd(yi) == 0) {
      out[[w]] <- data.frame(
        start = starts[w], end = ends[w], n = window_size,
        cor_p = NA_real_, bp_p = NA_real_, shapiro_p = NA_real_
      )
      next
    }
    
    df <- data.frame(xi = xi, yi = yi)
    fit <- stats::lm(xi ~ yi, data = df)
    
    cor_p <- tryCatch(stats::cor.test(xi, yi)$p.value, error = function(e) NA_real_)
    bp_p  <- tryCatch(lmtest::bptest(fit)$p.value, error = function(e) NA_real_)
    sh_p  <- tryCatch(stats::shapiro.test(stats::residuals(fit))$p.value, error = function(e) NA_real_)
    
    out[[w]] <- data.frame(
      start = starts[w], end = ends[w], n = window_size,
      cor_p = cor_p, bp_p = bp_p, shapiro_p = sh_p
    )
  }
  
  do.call(rbind, out)
}

library(crqa)
library(R.matlab)

setwd('')
data <- readMat('Experimental_R.mat')
WScnormBPM <- list(); WSchomBPM <- list(); WSclinBPM <- list(); WSvioBPM <- list()


numTS=49
WS = c(5, 10, 20, 40, 80, 160, 320)

for (window in WS){
  all_BPM <- list()
  all_RR  <- list()
  
  it <- 0L
  for (i in c(1:numTS)){
    if (is.numeric(data$dyadBPMmin[[i+2*numTS]][[1]][1,])&is.numeric(data$dyadBPMmin[[i+3*numTS]][[1]][1,])){
      
      ts1BPM=scale(t(data$dyadBPMmin[[i+2*numTS]][[1]]))[,1]
      ts2BPM=scale(t(data$dyadBPMmin[[i+3*numTS]][[1]]))[,1]
      # Ensure BPM lengths match too (just in case)
      mb <- min(length(ts1BPM), length(ts2BPM))
      ts1BPM <- ts1BPM[1:mb]
      ts2BPM <- ts2BPM[1:mb]
      
      wBPM <- windowed_tests(ts1BPM, ts2BPM, window_size = window, overlap = window/2)
      
      wBPM$dyad <- i
      all_BPM[[i]] <- wBPM
      
      it <- it + 1L
      print(it)
    }}
  
  BPM_df <- do.call(rbind, all_BPM)
  
  # --- Aggregate violation rates across ALL windows 
  
  cnormBPM <- mean(BPM_df$shapiro_p < 0.05, na.rm = TRUE) # Normality of residuals
  
  chomBPM  <- mean(BPM_df$bp_p < 0.05, na.rm = TRUE) # Homoscedasticity
  
  clinBPM  <- mean(BPM_df$cor_p > 0.05, na.rm = TRUE)  # not significant correlation
  
  # Same structure as your last line (note: your original mixed BPM norm with RR hom by accident)
  vioBPM <- mean((BPM_df$shapiro_p < 0.05) | (BPM_df$bp_p < 0.05) | (BPM_df$cor_p > 0.05), na.rm = TRUE)
  
  # Results:
  WScnormBPM = c(WScnormBPM, cnormBPM); WSchomBPM = c(WSchomBPM, chomBPM); WSclinBPM = c(WSclinBPM, clinBPM); WSvioBPM = c(WSvioBPM, vioBPM)
  
  head(BPM_df)}

# Non windowed linear assumption tests


hom_pBPM=c();lin_pBPM=c();norm_pBPM=c();hom_pRR=c();lin_pRR=c();norm_pRR=c()
it=0
for (i in c(1:49)){
  if (is.numeric(data$dyadBPMmin[[i+2*numTS]][[1]][1,])&is.numeric(data$dyadBPMmin[[i+3*numTS]][[1]][1,])){
    ts1BPM=scale(t(data$dyadBPMmin[[i+2*numTS]][[1]]))[,1]
    ts2BPM=scale(t(data$dyadBPMmin[[i+3*numTS]][[1]]))[,1]
    tsBPM=cbind(ts1BPM,ts2BPM)
    
    lmBPM=lm(ts1BPM~ts2BPM,data=data.frame(tsBPM))
    homBPM=bptest(lmBPM)
    hom_pBPM=cbind(hom_pBPM,homBPM$p.value)
    linBPM=cor.test(ts1BPM,ts2BPM)
    lin_pBPM=cbind(lin_pBPM,linBPM$p.value)
    normBPM=shapiro.test(residuals(lmBPM))
    norm_pBPM=cbind(norm_pBPM,normBPM$p.value)
    it=it+1
    print(it)
  }}

cnormBPM=sum(norm_pBPM<0.05)/it
chomBPM=sum(hom_pBPM<0.05)/it
clinBPM=sum(lin_pBPM>0.05)/it
vioBPM=sum(norm_pBPM<0.05|hom_pBPM<0.05|lin_pBPM>0.05)/it
WScnormBPM = c(WScnormBPM, cnormBPM); WSchomBPM = c(WSchomBPM, chomBPM); WSclinBPM = c(WSclinBPM, clinBPM); WSvioBPM = c(WSvioBPM, vioBPM)

# Save as a data frame
dfAss <- data.frame(t(rbind(WScnormBPM, WSchomBPM, WSclinBPM, WSvioBPM)), c(WS, 'Full'))
column_namesAss <- c("Norm", "Hom", "Lin", "Tot", "WS")
colnames(dfAss) <- column_namesAss
save('dfAss',file = 'Exp_Ass.RData')


# =============================================================================
# Real Data Analysis
# =============================================================================

# Load required packages and data
library(crqa)
library(R.matlab)
library(lmtest)
results=list()
setwd('')
data=readMat('Experimental_R.mat')
UE=c(); CC=c();  E=c(); WCC=c()

# MdRQA parameters
# Number of embedding dimensions
emb=15
# Delay
tau=2
# Radius Embedded
re=c(4,5,6)
# Radius Unembedded
ru=c(0.4,0.7,0.9)


#' Windowed zero-lag cross-correlation (Pearson) with overlap, Fisher z, then average
#'
#' x,y Numeric vectors of equal length.
#' window_size Integer >= 2. Window length in data points.
#' overlap Integer in [0, window_size-1]. Number of points shared between consecutive windows.
#' min_complete Integer >= 2. Minimum number of non-NA paired samples required per window.
#' use. Which stat to utilize for the analysis. Passed to stats::cor(). Default "pairwise.complete.obs".
#' method Correlation method. Default "pearson".
#' clamp_r Numeric in (0,1). Clamps |r| to < clamp_r to avoid infinite Fisher z at |r|=1.
#' returns a list with:
#'   - avg_r: average Fisher-z across windows, inverse-transformed back to r
#'   - avg_z: mean Fisher z across windows
#'   - windows: data.frame with per-window indices, n_complete, r, z
#'
wcc_fisher <- function(
    x,
    y,
    window_size,
    overlap = 0L,
    min_complete = 2L,
    use = "pairwise.complete.obs",
    method = "pearson",
    clamp_r = 0.999999
) {
  # ---- validation ----
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors.")
  if (length(x) != length(y)) stop("x and y must have the same length.")
  n <- length(x)
  
  window_size <- as.integer(window_size)
  overlap <- as.integer(overlap)
  min_complete <- as.integer(min_complete)
  
  if (window_size < 2L) stop("window_size must be >= 2.")
  if (overlap < 0L || overlap >= window_size) stop("overlap must be in [0, window_size-1].")
  if (min_complete < 2L) stop("min_complete must be >= 2.")
  if (window_size > n) stop("window_size cannot exceed the length of the vectors.")
  
  step <- window_size - overlap
  starts <- seq.int(1L, n - window_size + 1L, by = step)
  ends <- starts + window_size - 1L
  
  # ---- compute per-window r and Fisher z ----
  r_vals <- rep(NA_real_, length(starts))
  z_vals <- rep(NA_real_, length(starts))
  n_complete <- rep(NA_integer_, length(starts))
  
  for (i in seq_along(starts)) {
    idx <- starts[i]:ends[i]
    xi <- x[idx]
    yi <- y[idx]
    
    ok <- is.finite(xi) & is.finite(yi)
    n_ok <- sum(ok)
    n_complete[i] <- n_ok
    
    if (n_ok >= min_complete) {
      r <- suppressWarnings(stats::cor(xi, yi, use = use, method = method))
      # Handle edge cases: constant window or cor() returning NA
      if (is.finite(r)) {
        # Clamp to avoid atanh(±1) = ±Inf
        r <- max(min(r, clamp_r), -clamp_r)
        r_vals[i] <- r
        z_vals[i] <- atanh(r)
      }
    }
  }
  
  valid <- is.finite(z_vals)
  if (!any(valid)) {
    return(list(
      avg_r = NA_real_,
      avg_z = NA_real_,
      windows = data.frame(
        start = starts,
        end = ends,
        n_complete = n_complete,
        r = r_vals,
        z = z_vals
      )
    ))
  }
  
  avg_z <- mean(z_vals[valid])
  avg_r <- tanh(avg_z)
  
  list(
    avg_r = avg_r,
    avg_z = avg_z,
    windows = data.frame(
      start = starts,
      end = ends,
      n_complete = n_complete,
      r = r_vals,
      z = z_vals
    )
  )
}


numTS=length(data$dyadBPMmin)/4
real_vec = 2:(numTS+1)
it=1
for (i in c(1:numTS)){
  reaID = real_vec[i]
  if (is.numeric(data$dyadBPMmin[[i+2*numTS]][[1]][1,])&is.numeric(data$dyadBPMmin[[i+3*numTS]][[1]][1,])){
    ts1BPMmin=scale(t(data$dyadBPMmin[[i+2*numTS]][[1]]))[,1]
    ts2BPMmin=scale(t(data$dyadBPMmin[[i+3*numTS]][[1]]))[,1]
    
    # CCF
    CCFBPMmin=ccf(ts1BPMmin,ts2BPMmin,2,pl=FALSE)
    
    currentCC=c(i,i,CCFBPMmin$acf[3])
    CC=rbind(CC,currentCC)
    
    # WCC
    currentWCC = c(i, wcc_fisher(ts1BPMmin, ts2BPMmin, 20, 10)$avg_r)
    WCC=rbind(WCC,currentWCC)
    
    for(j in c(1:3)){
      # MdRQA
      # Binding
      tsBPMmin=cbind(ts1BPMmin,ts2BPMmin)
      # Unembedded
      MdRQABPMmin=crqa(tsBPMmin,tsBPMmin,rescale=0,radius=ru[j],tw = 1, delay = 1,embed=1,normalize = 0,method='mdcrqa')
      # Embedded
      MdRQABPMmin_emb=crqa(tsBPMmin,tsBPMmin,rescale=0,radius=re[j], tw = 1, delay = tau,embed=emb,normalize = 0,method='mdcrqa')
      
      currentUE=c(reaID,reaID,MdRQABPMmin$RR,MdRQABPMmin$LAM,ru[j])
      UE=rbind(UE,currentUE)
      currentE=c(MdRQABPMmin_emb$RR,MdRQABPMmin_emb$LAM,re[j])
      E=rbind(E,currentE)}
    print(it)
    it=it+1
  }}
save(emb,tau,ru,re,CC,WCC,UE,E,file = 'Experimental_real_R.RData')

condRQA=c(rep('Real',48*3))
dfRQA <- data.frame(condRQA,UE,E)
column_namesRQA <- c("Condition", "Subject","Sur",  "%REC(U)_BPMmin","%LAM(U)_BPMmin",'r(U)', "%REC(E)_BPMmin", "%LAM(E)_BPMmin",'r(E)')
colnames(dfRQA) <- column_namesRQA

condCC=c(rep('Real',48))
dfCC <- data.frame(condCC,CC)
column_namesCC <- c("Condition", "Subject","Sur", "CC_BPMmin")
colnames(dfCC) <- column_namesCC

condWCC=c(rep('Real',48))
dfWCC <- data.frame(condWCC,WCC)
column_namesWCC <- c("Condition", "Subject", "WCC_BPMmin")
colnames(dfWCC) <- column_namesWCC

# Duplicate to match 600 cells of mdRQA
indices <- rep(1:nrow(dfCC), each = 3)

dfCC3 <- dfCC[indices, ]
dfWCC3 <- dfWCC[indices, ]


save('dfRQA',file = 'Exp_RealMdRQA_R.RData')
save('dfCC',file = 'Exp_RealCC_R.RData')
save('dfWCC',file = 'Exp_RealWCC_R.RData')

save('dfCC3',file = 'Exp_RealCC3_R.RData')
save('dfWCC3',file = 'Exp_RealWCC3_R.RData')

dfCombined_real=cbind(dfCC3,dfWCC3,dfRQA)
save('dfCombined_real',file = 'Exp_RealCombined_R.RData')


# =============================================================================
# Surrogate Data Analysis
# =============================================================================


# Load required packages and data
library(crqa)
library(R.matlab)
library(lmtest)
results=list()
setwd('C:/Users/alont/Desktop/Work/Leuphana/mtdComp/RevPlos/Experimental')
data=readMat('Experimental_R.mat')
UE=c(); CC=c(); E=c(); WCC=c()
rm()
emb=15
tau=2

# Radius Embedded
re=c(4,5,6)
# Radius Unembedded
ru=c(0.4,0.7,0.9)

numTS=length(data$dyadBPMmin)/4
orig_real_vec = 2:(numTS+1)
real_vec <- orig_real_vec[real_vec != 29]
# Creating a new surrogate vector
# sur_vec = sample(real_vec)
# while (any(real_vec == sur_vec)){
#   sur_vec=sample(real_vec)
# }
# # Add 28 into the right position
# sur_vec <- append(sur_vec,29,after = which(orig_real_vec == 29)-1)
# 
# save(sur_vec,file = 'FS2sur_vector.RData')
load('FS2sur_vector.RData')
real_vec <- append(real_vec,29,after = which(orig_real_vec == 29)-1)


it=1
for (i in c(1:length(real_vec))){
  surID = sur_vec[i]
  sur=which(sur_vec[i] == real_vec)
  if (is.numeric(data$dyadBPMmin[[sur+2*numTS]][[1]][1,])&is.numeric(data$dyadBPMmin[[sur+3*numTS]][[1]][1,])&is.numeric(data$dyadBPMmin[[i+2*numTS]][[1]][1,])&is.numeric(data$dyadBPMmin[[i+3*numTS]][[1]][1,])){
    reaID= real_vec[i] 
    # The analysis utilizes the BPM with the smallest window (minimum RR)
    ts1BPMmin=t(data$dyadBPMmin[[i+2*numTS]][[1]])
    ts2BPMmin=t(data$dyadBPMmin[[sur+3*numTS]][[1]])
    minLengthBPMmin=c(1:min(c(length(ts1BPMmin),length(ts2BPMmin))))
    ts1BPMmin=scale(ts1BPMmin[minLengthBPMmin])[,1]; ts2BPMmin=scale(ts2BPMmin[minLengthBPMmin])[,1]
    
    # CCF
    CCFBPMmin=ccf(ts1BPMmin,ts2BPMmin,2,pl=FALSE)
    
    currentCC=c(reaID,surID,CCFBPMmin$acf[3])
    CC=rbind(CC,currentCC)
    
    # WCC
    currentWCC = c(wcc_fisher(ts1BPMmin, ts2BPMmin, 20, 10)$avg_r)
    WCC=rbind(WCC,currentWCC)
    
    for(j in c(1:3)){
      # MdRQA
      # Binding
      tsBPMmin=cbind(ts1BPMmin,ts2BPMmin)
      # Unembedded
      MdRQABPMmin=crqa(tsBPMmin,tsBPMmin,rescale=0,radius=ru[j],tw = 1, delay = 1,embed=1,normalize = 0,method='mdcrqa')
      # Embedded
      MdRQABPMmin_emb=crqa(tsBPMmin,tsBPMmin,rescale=0,radius=re[j], tw = 1, delay = tau,embed=emb,normalize = 0,method='mdcrqa')
      
      
      currentUE=c(reaID,surID,MdRQABPMmin$RR,MdRQABPMmin$LAM,ru[j])
      UE=rbind(UE,currentUE)
      currentE=c(MdRQABPMmin_emb$RR,MdRQABPMmin_emb$LAM,re[j])
      E=rbind(E,currentE)}
    print(it)
    it=it+1
  }}
save(emb,tau,ru,re,CC,WCC,UE,E,file = 'Experimental_sur_R.RData')
# Create data frame
condRQA=c(rep('Sur',48*3))
dfRQA <- data.frame(condRQA,UE,E)
column_namesRQA <- c("Condition", "Subject","Sur","%REC(U)_BPMmin", "%LAM(U)_BPMmin",'r(U)', "%REC(E)_BPMmin", "%LAM(E)_BPMmin",'r(E)')
colnames(dfRQA) <- column_namesRQA

condCC=c(rep('Sur',48))
dfCC <- data.frame(condCC,CC)
column_namesCC <- c("Condition", "Subject","Sur",  "CC_BPMmin")
colnames(dfCC) <- column_namesCC

condWCC=c(rep('Sur',48))
dfWCC <- data.frame(condWCC,WCC)
column_namesWCC <- c("Condition","WCC_BPMmin")
colnames(dfWCC) <- column_namesWCC

# Duplicate to match 600 cells of mdRQA
indices <- rep(1:nrow(dfCC), each = 3)

dfCC3 <- dfCC[indices, ]
dfWCC3 <- dfWCC[indices, ]

save('dfRQA',file = 'Exp_SurMdRQA_R.RData')
save('dfCC',file = 'Exp_SurCC_R.RData')
save('dfWCC',file = 'Exp_SurWCC_R.RData')
save('dfCC3',file = 'Exp_SurCC3_R.RData')
save('dfWCC3',file = 'Exp_SurWCC3_R.RData')

dfCombined_sur=cbind(dfCC3,dfWCC3,dfRQA)
save('dfCombined_sur',file = 'Exp_SurCombined_R.RData')

dfCombined=cbind(dfCombined_sur,dfCombined_real)
save('dfCombined',file = 'Exp_CombinedPaired_R.RData')


