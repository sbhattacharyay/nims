function (PA, label, flag, demo = NA, method = "zipln", 
          time.range = c("09:00", "20:59"), K = 3, D = 5, 
          mark.missing = 0, thresh = 10000, graph.diagnostic = TRUE, 
          seed = 1234, m = 5, maxit = 6) 
{
  PA = as.matrix(PA)
  label = as.data.frame(label)
  flag = as.matrix(flag)
  if (is.null(nrow(demo))) {
    demo.include = FALSE
    print("no demo...")
  }
  if (!is.null(nrow(demo))) {
    demo.include = TRUE
    demo = as.data.frame(demo)
  }
  nw = mark.missing
  w = abs(1 - nw)
  start.time = as.POSIXct("2015-01-01")
  seq.time = seq.POSIXt(start.time, length.out = 1440, by = "min")
  min.seq = format(seq.time, "%H:%M")
  min.range = which(min.seq %in% time.range)
  daytime = min.range[1]:min.range[2]
  PA[PA >= thresh] = thresh
  if (demo.include) {
    key1 = colnames(label)[1]
    key2 = colnames(demo)[1]
    labelDemo = merge(label, demo, by.x = key1, by.y = key2, 
                      all.x = TRUE)
    demo.daily = labelDemo[, colnames(demo)]
    if (nrow(PA) != nrow(demo.daily)) {
      stop("Dimension does not match! Please check the demographic data.")
    }
    if (sum(is.na(demo)) != 0) {
      stop("There are missing values (NA) in demo. Imputation is needed, otherwise set 'demo.include=FALSE'.")
    }
  }
  daylabel = label[, 2]
  wk = ifelse((daylabel%%7) >= 2, 1, 0)
  wk = as.factor(wk)
  if (demo.include) {
    xmat = data.frame(demo.daily[, -1], wk)
  }
  else {
    xmat = data.frame(wk)
  }
  print("First iteration starts with zip + pmm... ")
  intimp = PA
  for (s in daytime) {
    cat(paste(min.seq[s], "..."))
    t1 = s
    t0 = s - 1
    if (s == daytime[1]) {
      y0 = PA[, t0]
      y0[flag[, t0] == nw] = NA
      icd.y0 = data.frame(y0, xmat)
      imp.y0 = mice(icd.y0, seed = seed, method = "pmm", 
                    maxit = maxit, m = 1, printFlag = FALSE)
      cd.y0 = complete(imp.y0, 1)$y0
    }
    else {
      cd.y0 = cd.y1
    }
    y1 = PA[, t1]
    r1 = (flag[, t1] == w)
    y1[!r1] = NA
    icd.y1 = data.frame(y1, ln.y0 = log(cd.y0 + 1), xmat)
    imp.y1 = mice(icd.y1, seed = seed, method = "2l.zip.pmm", 
                  maxit = maxit, m = 1, printFlag = FALSE)
    cd.y1 = complete(imp.y1)$y1
    if (graph.diagnostic == TRUE) {
      plot(log(cd.y1 + 1), log(cd.y0 + 1), col = ifelse(r1, 
                                                        3, 2), pch = ifelse(r1, 1, 8), xlab = min.seq[t1], 
           ylab = min.seq[t0], main = "log(count+1)")
      legend("topleft", legend = c("observed", 
                                   "imputed"), col = c(3, 2), pch = c(1, 8))
    }
    intimp[, t1] = complete(imp.y1)$y1
  }
  print("done.")
  print(paste("Preparing the zipln imputation with K=", 
              K, "lag and lead"))
  zmat = matrix(NA, nrow(intimp), ncol(intimp))
  for (s in (daytime[1] - K):(tail(daytime, 1) + K)) {
    yt = round(intimp[, s])
    data.t = data.frame(yt, xmat)
    zip.t = zeroinfl(yt ~ . | ., data = data.t)
    lam.t = predict(zip.t, type = "count")
    zmat[, s] = log(yt + 1) - log(lam.t + 1)
    cat(".", s)
  }
  listimp = vector("list", m)
  names(listimp) = paste("imp", 1:m, sep = "")
  for (k in 1:m) listimp[[k]] = intimp
  print("done")
  print(paste("Second iteration starts with", method, 
              " ... "))
  for (s in daytime) {
    cat(paste(min.seq[s], "..."))
    if (s == daytime[1]) {
      cd.yt_1 = intimp[, s - 1]
    }
    else {
      cd.yt_1 = cd.yt
    }
    tvec = c(s + (-K:K))
    names(tvec) = c(paste("t_", K:1, sep = ""), 
                    "t0", paste("t.", 1:K, sep = ""))
    ytmat = intimp[, tvec]
    ytmat1 = ytmat
    ytmat1[flag[, tvec] == nw] = NA
    colnames(ytmat) = c(paste("yt_", K:1, sep = ""), 
                        "yt", paste("yt.", 1:K, sep = ""))
    colnames(ytmat1) = colnames(ytmat)
    ry = (flag[, s] == w)
    zs = zmat[, tvec]
    colnames(zs) = c(paste("zt_", K:1, sep = ""), 
                     "zt", paste("zt.", 1:K, sep = ""))
    icd.yt = data.frame(yt = ytmat1[, (K + 1)], xmat, ln.yt_1 = log(cd.yt_1 + 
                                                                      1))
    ini = mice(icd.yt, seed = seed, maxit = 0)
    predmat = ini$predictorMatrix
    predmat[1, "ln.yt_1"] = 3
    if (method == "zipln") {
      imp.yt = mice(icd.yt, seed = seed, method = "2l.zipln", 
                    predictorMatrix = predmat, maxit = maxit, m = m, 
                    K = K, zs = zs, printFlag = FALSE)
    }
    if (method == "zipln.pmm") {
      imp.yt = mice(icd.yt, seed = seed, method = "2l.zipln.pmm", 
                    predictorMatrix = predmat, maxit = maxit, m = m, 
                    K = K, zs = zs, D = D, printFlag = FALSE)
    }
    cd.yt = complete(imp.yt, m)$yt
    if (graph.diagnostic == TRUE) {
      plot(log(cd.yt + 1), log(cd.yt_1 + 1), col = ifelse(r1, 
                                                          3, 2), pch = ifelse(r1, 1, 8), xlab = min.seq[s], 
           ylab = min.seq[s - 1], main = "log(count+1)")
      legend("topleft", legend = c("observed", 
                                   "imputed"), col = c(3, 2), pch = c(1, 8))
    }
    for (j in 1:m) listimp[[j]][, s] = complete(imp.yt, j)$yt
  }
  print(paste("Imputation complete!", m, "datasets are created."))
  return(listimp)
}