# =============================================================================
# Assumption Tests
# =============================================================================


# Linear correlation assumption tests - Linearity, Homoscedasticity, Norm of res
library(crqa)
library(R.matlab)
library(lmtest)
setwd('')
data=readMat('Synthetic_Data.mat')


# Windowed versions of:
# - Pearson correlation p-value
# - Breusch–Pagan p-value (Homoscedasticity)
# - Shapiro p-value (normality of residuals)
# and returns BOTH:
# (a) per-window results (one row per dyad × window)
# (b) aggregated violation rates

library(lmtest)

windowed_tests <- function(x, y, window_size, overlap = 0L) {
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
    
    # Protect from series with sd=0 
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

# Non windowed linear assumption tests

library(crqa)
library(R.matlab)

setwd('')
data <- readMat('Synthetic_Data.mat')
WScnormBPM <- list(); WSchomBPM <- list(); WSclinBPM <- list(); WSvioBPM <- list()
WS = c(5, 10, 20, 40, 80, 160, 320)

for (window in WS){
  all_BPM <- list()
  all_RR  <- list()
  
  it <- 0L
  for (i in 1:300) {
    ts1BPM <- as.numeric(t(data$dyadBPMmin[[i]][[1]]))
    ts2BPM <- as.numeric(t(data$dyadBPMmin[[i + 300]][[1]]))
    
    # Ensure BPM lengths match too 
    mb <- min(length(ts1BPM), length(ts2BPM))
    ts1BPM <- ts1BPM[1:mb]
    ts2BPM <- ts2BPM[1:mb]
    
    wBPM <- windowed_tests(ts1BPM, ts2BPM, window_size = window, overlap = window/2)
    
    wBPM$dyad <- i
    all_BPM[[i]] <- wBPM
    
    it <- it + 1L
    print(it)
  }
  
  BPM_df <- do.call(rbind, all_BPM)
  
  cnormBPM <- mean(BPM_df$shapiro_p < 0.05, na.rm = TRUE) # Normality of residuals
  
  chomBPM  <- mean(BPM_df$bp_p < 0.05, na.rm = TRUE) # Homoscedasticity 
  
  clinBPM  <- mean(BPM_df$cor_p > 0.05, na.rm = TRUE)  # not significant correlation
  
  vioBPM <- mean((BPM_df$shapiro_p < 0.05) | (BPM_df$bp_p < 0.05) | (BPM_df$cor_p > 0.05), na.rm = TRUE)
  
  # Results:
  WScnormBPM = c(WScnormBPM, cnormBPM); WSchomBPM = c(WSchomBPM, chomBPM); WSclinBPM = c(WSclinBPM, clinBPM); WSvioBPM = c(WSvioBPM, vioBPM)
  
  head(BPM_df)}

hom_pBPM=c();lin_pBPM=c();norm_pBPM=c();hom_pRR=c();lin_pRR=c();norm_pRR=c()
it=0
for (i in c(1:300)){
  ts1BPM=t(data$dyadBPMmin[[i]][[1]])
  ts2BPM=t(data$dyadBPMmin[[i+300]][[1]])
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
}
cnormBPM=sum(norm_pBPM<0.05)/it
chomBPM=sum(hom_pBPM<0.05)/it
clinBPM=sum(lin_pBPM>0.05)/it
vioBPM=sum(norm_pBPM<0.05|hom_pRR<0.05|lin_pBPM>0.05)/it
cnormBPM; chomBPM; clinBPM; vioBPM
WScnormBPM = c(WScnormBPM, cnormBPM); WSchomBPM = c(WSchomBPM, chomBPM); WSclinBPM = c(WSclinBPM, clinBPM); WSvioBPM = c(WSvioBPM, vioBPM)

dfAss <- data.frame(t(rbind(WScnormBPM, WSchomBPM, WSclinBPM, WSvioBPM)), c(WS, 'Full'))
column_namesAss <- c("Norm", "Hom", "Lin", "Tot", "WS")
colnames(dfAss) <- column_namesAss
save('dfAss',file = 'Syn_Ass.RData')


# =============================================================================
# Synthetic data analysis
# =============================================================================


# Load required packages and data
library(crqa)
library(R.matlab)
library(lmtest)
results=list()
setwd('')
data=readMat('Synthetic_Data.mat')
UE=c(); CC=c();  E=c(); WCC=c()

# MdRQA parameters

# Radius Embedded
re=c(8)
# Radius Unembedded
ru=c(0.7)
# Number of embedding dimensions
emb=25
# Delay
tau=2

#' Windowed zero-lag cross-correlation (Pearson) with overlap, Fisher z, then average
#'
#' x,y Numeric vectors of equal length.
#' window_size Integer >= 2. Window length in data points.
#' overlap Integer in [0, window_size-1]. Number of points shared between consecutive windows.
#' min_complete Integer >= 2. Minimum number of non-NA paired samples required per window.
#' use. Which stat to utilize for the analysis. Passed to stats::cor(). Default "pairwise.complete.obs".
#' method Correlation method. Default "pearson".
#' clamp_r Numeric in (0,1). Clamps |r| to < clamp_r to avoid infinite Fisher z at |r|=1.
#' return a list with:
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

it=1
for (i in c(1:300)){
  # The analysis utilizes the BPM with the smallest window (minimum RR)
  ts1BPMmin=scale(t(data$dyadBPMmin[[i]][[1]]))[,1] 
  ts2BPMmin=scale(t(data$dyadBPMmin[[i+300]][[1]]))[,1]
  
  
  # CCF
  CCFBPMmin=ccf(ts1BPMmin,ts2BPMmin,2,pl=FALSE)
  currentCC=c(i,CCFBPMmin$acf[3])
  CC=rbind(CC,currentCC)
  
  # WCC
  currentWCC = c(i, wcc_fisher(ts1BPMmin, ts2BPMmin, 40, 20)$avg_r)
  WCC=rbind(WCC,currentWCC)
  
  # MdRQA
  # Binding
  tsBPMmin=cbind(ts1BPMmin,ts2BPMmin)
  
  # Unembedded
  MdRQABPMmin=crqa(tsBPMmin,tsBPMmin,rescale=0,radius=ru,tw = 1, delay = 1,embed=1,normalize = 0,method='mdcrqa')
  # Embedded
  MdRQABPMmin_emb=crqa(tsBPMmin,tsBPMmin,rescale=0,radius=re, tw = 1, delay = tau,embed=emb,normalize = 0,method='mdcrqa')
  
  
  currentUE=c(i,MdRQABPMmin$RR,MdRQABPMmin$LAM,ru)
  UE=rbind(UE,currentUE)
  currentE=c(i,MdRQABPMmin_emb$RR,MdRQABPMmin_emb$LAM,re)
  E=rbind(E,currentE)
  print(it)
  it=it+1}

save(emb,tau,ru,re,CC,WCC,UE,E,file = 'Synthetic_flex_R.RData')

# Create a data frame

condRQA=c(rep('Coupled',100),rep('Uncoupled',100),rep('Flex',100))
dfRQA <- data.frame(condRQA,UE,E[,2:4])
column_namesRQA <- c("Condition", "Subject",  "%REC(U)_BPMmin", "%LAM(U)_BPMmin",'r(U)',"%REC(E)_BPMmin", "%LAM(E)_BPMmin",'r(E)')
colnames(dfRQA) <- column_namesRQA

condCC=c(rep('Coupled',100),rep('Uncoupled',100),rep('Flex',100))
dfCC <- data.frame(condCC,CC)
column_namesCC <- c("Condition", "Subject", "CC_BPMmin")
colnames(dfCC) <- column_namesCC

condWCC=c(rep('Coupled',100),rep('Uncoupled',100),rep('Flex',100))
dfWCC <- data.frame(condWCC,WCC)
column_namesWCC <- c("Condition", "Subject", "WCC_BPMmin")
colnames(dfWCC) <- column_namesWCC

save('dfRQA',file = 'SyntheticMdRQA_R.RData')
save('dfCC',file = 'SyntheticCC_R.RData')
save('dfWCC',file = 'SyntheticWCC_R.RData')


dfCombined=cbind(dfCC,dfWCC,dfRQA)
save('dfCombined',file = 'Comb_R.RData')
