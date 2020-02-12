#
# yo_utils.R
#
# Utility functions for analysis of the young/old dataset
#
# 2019-11-07  WTR


# gating functions
gate_clean = function(ff, params = c("SSC-A", "KI67FITC", "CD127BV421", "CCR4AF647", "PD1PE"), show = FALSE, show.fn = NULL) {
  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 800, height = 800)
    }
    res = clean.fp(ff = ff, parameters = params, show = show)
    if (!is.null(show.fn)) {
      dev.off()
    }
  } else {
    res = clean.fp(ff = ff, parameters = params, show = show)
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 1))
  }
  ff = Subset(res, rectangleGate("clean" = c(0.5, Inf)))
  exprs(ff) = exprs(ff)[, -which(colnames(ff) == "clean")]

  ff
}

gate_singlet = function(ff, show = FALSE, show.fn = NULL, ...) {
  n_orig = nrow(ff)
  p_singlet = c("FSC-W", "SSC-W")
  bb_singlet = blob.boundary(ff, parameters = p_singlet, location = c(1, 1), height = .05)
  ssc_thresh = max(bb_singlet[, 2])
  fsc_thresh = max(bb_singlet[, 1])
  singlet_gate = rectangleGate("SSC-W" = c(-Inf, ssc_thresh), "FSC-W" = c(-Inf, fsc_thresh))
  ff_singlet = Subset(ff, singlet_gate)
  n_singlet = nrow(ff_singlet)

  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 600, height = 600)
    }
    pplot(ff, p_singlet, tx = 'linear', ty = 'linear', xlim =c(0, 5), ylim = c(0, 5), ...)
    lines(bb_singlet, col = 'red', lwd = 2)
    xline(fsc_thresh, lty = 'dotdash')
    yline(ssc_thresh, lty = 'dotdash')
    text(0.5, 4.5, labels = sprintf("%.2f%% Singlets", 100 * n_singlet / n_orig), pos = 4)
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 1))

    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_singlet
}

# look for the LIVEDEAD negative, CD3 positive blob
gate_live_cd3 = function(ff, show = FALSE, show.fn = NULL, ...) {
  params = c("CD3Q605", "LIVEDEAD")

  pre_cd3 = bx(1000)
  pre_live = bx(2000)
  pre = Subset(ff, rectangleGate("CD3Q605" = c(pre_cd3, Inf), "LIVEDEAD" = c(-Inf, pre_live)))
  bb = blob.boundary(pre, parameters = params, location = bx(c(5000, 500)), height = 0.2)
  # inflate generously
  idist = 0.25
  bb_infl = inflate.contour(get.hull(bb), dist = idist)
  gate_live = polygonGate(.gate = bb_infl)
  ff_live = Subset(ff, gate_live)

  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 600, height = 600)
    }
    pplot(ff, plist = params, xlim = c(0, 5.4), ylim = c(0, 5.4), ...)
    lines(bb)
    lines(bb_infl, lwd = 3, col = 'red')
    xline(pre_cd3, lty = 'dotdash')
    yline(pre_live, lty = 'dotdash')
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 1))

    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_live
}

# so as not to exclude CD3-, just try to exclude LIVEDEAD+ events
#
gate_live = function(ff, no.debris = TRUE, pre_x = 0.5, pre_y = 0.5,
                     show = FALSE, show.fn = NULL) {

  if (no.debris) {
    pre = Subset(ff, rectangleGate("FSC-A" = c(pre_x, Inf), "SSC-A" = c(pre_y, Inf)))
  } else {
    pre = ff
  }
  kde_3   = normalize.kde(bkde(exprs(pre)[,"CD3Q605"], band = 0.1, grid = 1001))
  kde_11b = normalize.kde(bkde(exprs(pre)[,"CD11BAPCCY7"], band = 0.1, grid = 1001))
  kde_14  = normalize.kde(bkde(exprs(pre)[,"CD14Q800"], band = 0.1, grid = 1001))

  pk_3 = max(find.local.maxima(kde_3, thresh = .01)$x)
  pk_11b = max(find.local.maxima(kde_11b, thresh = .01)$x)
  pk_14 = max(find.local.maxima(kde_14, thresh = .01)$x)

  # get thresholds
  lmin_3 = find.local.minima(kde_3, thresh = .005)$x
  lmin_11b = find.local.minima(kde_11b, thresh = .005)$x
  lmin_14 = find.local.minima(kde_14, thresh = .005)$x
  if(length(lmin_3) == 0) {lmin_3 = bx(2000)}
  if(length(lmin_11b) == 0) {lmin_11b = bx(2000)}
  if(length(lmin_14) == 0) {lmin_14 = bx(2000)}
  thresh_3   = max(lmin_3)
  thresh_11b = max(lmin_11b)
  thresh_14  = max(lmin_14)

  # make the OR gate, then see what's what with LIVEDEAD
  g_3   = rectangleGate("CD3Q605" = c(thresh_3, Inf))
  g_11b = rectangleGate("CD11BAPCCY7" = c(thresh_11b, Inf))
  g_14  = rectangleGate("CD14Q800" = c(thresh_14, Inf))
  g_or = g_14 | g_11b | g_3
  res = Subset(ff, g_or)

  # find the two live-cell blob centers
  ysep = 2.5
  xsep = bx(4000)
  top = Subset(res, rectangleGate("SSC-A" = c(ysep, Inf)), "LIVEDEAD" = c(-Inf, xsep))
  bot = Subset(res, rectangleGate("SSC-A" = c(-Inf, ysep)), "LIVEDEAD" = c(-Inf, xsep))
  bb1 = blob.boundary(ff = bot, parameters = c("LIVEDEAD", "SSC-A"), location = c(bx(600), 1.0), gridsize = c(501, 501), height = .2)
  bb2 = blob.boundary(ff = top, parameters = c("LIVEDEAD", "SSC-A"), location = c(bx(2000), 3.2), gridsize = c(501, 501), height = .2)
  cen1 = centroid(bb1)
  cen2 = centroid(bb2)

  # connect the two centers with a line, then displace to the right by an amount
  # depending on the "radius of the lower blob
  rad = (max(bb1[,1]) - min(bb1[,1])) / 2
  infl = 0.4
  rad = rad + infl
  fitdat = data.frame(x = c(cen1[1], cen2[1]), y = c(cen1[2], cen2[2]))
  mod = lm(x ~ y, data = fitdat)
  corners = data.frame(y = c(0, 5.5))
  pred = predict(mod, corners)
  corners = cbind(x = pred, corners)
  corners$x = corners$x + rad

  # complete the gating polygon
  corners[3, ] = c(-0.2, corners[2, 2])
  corners[4, ] = c(-0.2, corners[1, 2])
  corners[5, ] = corners[1, ]
  corners = as.matrix(corners)
  colnames(corners) = c("LIVEDEAD", "SSC-A")
  pgate = polygonGate(.gate = corners)
  live = Subset(pre, pgate)

  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 600, height = 800)
    }
    par(mfrow = c(3, 2))

    pplot(ff, c("FSC-A","SSC-A"), tx = 'linear', ty = 'linear')
    xline(pre_x, col = 'red')
    yline(pre_y, col = 'red')

    plot(kde_3, type = 'l', col = 'red', xlim = c(0, 5.4), xaxt = 'n', xlab = '', ylab = '', main = "CD3")
    ax(axis = 1, type = 'biexp')
    xline(pk_3, lty = 'dotdash')
    xline(thresh_3, lty = 'dotdash', col = 'blue')

    plot(kde_11b, type = 'l', col = 'red', xlim = c(0, 5.4), xaxt = 'n', xlab = '', ylab = '', main = "CD11b")
    ax(axis = 1, type = 'biexp')
    xline(pk_11b, lty = 'dotdash')
    xline(thresh_11b, lty = 'dotdash', col = 'blue')

    plot(kde_14, type = 'l', col = 'red', xlim = c(0, 5.4), xaxt = 'n', xlab = '', ylab = '', main = "CD14")
    ax(axis = 1, type = 'biexp')
    xline(pk_14, lty = 'dotdash')
    xline(thresh_14, lty = 'dotdash', col = 'blue')

    pplot(ff, c("LIVEDEAD", "SSC-A"), ty = 'linear', main = "original", xlim = c(-1, 5.4))
    pplot(res, c("LIVEDEAD", "SSC-A"), ty = 'linear', main = "live gate", xlim = c(-1, 5.4))

    lines(bb1)
    lines(bb2)
    lines(corners, lwd = 3)

    par(mfrow = c(1, 1), mar = c(5, 4, 4, 1))

    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  live
}

# since we included CD3 in gate_live(), this is just a bit of cleanup
gate_scat = function(ff, show = FALSE, show.fn = NULL) {
  params = c("FSC-A", "SSC-A")

  pre = Subset(ff, rectangleGate("FSC-A" = c(0.4, Inf)))
  bb = blob.boundary(pre, parameters = params, location = c(1, 1), height = 0.2)
  # inflate generously
  idist = 0.25
  bb_infl = inflate.contour(get.hull(bb), dist = idist)
  gate_scat = polygonGate(.gate = bb_infl)
  ff_scat = Subset(ff, gate_scat)

  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 600, height = 600)
    }
    pplot(ff, c("FSC-A","SSC-A"), tx = 'linear', ty = 'linear', xlim = c(0, 5.4), ylim = c(0, 5.4))
    lines(bb)
    lines(bb_infl, lwd = 3, col = 'red')
    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_scat
}

################################################################################
################################################################################
# pre-gating stimulation data
################################################################################
################################################################################
gate_lymph_stim = function(ff, is.day8 = FALSE, show = FALSE, show.fn = NULL, ...) {
  params = c("FSC-A", "SSC-A")

  # find a FSC-A threshold to exclude debris, etc
  kde = bkde(exprs(ff)[, params[1]], band = .1)
  kde$y = kde$y / max(kde$y)
  res = find.local.minima(kde = kde)$x
  default_xthresh = 1.0
  if (length(res) > 0) {
    xthresh = max(max(res[res < 2.0]), default_xthresh)  # largest minimum below 2.0, but at least 1.0
  } else {
    xthresh = default_xthresh
  }
  pre = Subset(ff, rectangleGate("FSC-A" = c(xthresh, Inf)))

  # estimate the angle of the major axis of the lymphs
  x = exprs(pre)[, params[1]]
  y = exprs(pre)[, params[2]]
  mod = lm(y ~ x)
  degrees = atan2(mod$coefficients[2], 1) * 180 / pi

  # estimate the ratio of the blob variances
  rat_blob = 2 * ((IQR(x) / IQR(y)) ^ 2)

  if (is.day8) {
    band = 0.02
    height = 0.1
  } else {
    band = 0.01
    height = 0.66
  }
  bb = blob.boundary(pre, parameters = params, rotate = -degrees, bandwidth = band * c(rat_blob, 1),
                     location = c(2, 0.5), gridsize = c(1001, 1001), height = height)
  bb = get.hull(bb)
  bb = smooth.contour(bb)
  # inflate generously
  idist = 0.15
  bb_infl = inflate.contour(get.hull(bb), dist = idist)
  gate_scat = polygonGate(.gate = bb_infl)
  ff_scat = Subset(ff, gate_scat)

  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 600, height = 600)
    }
    pplot(ff, c("FSC-A","SSC-A"), tx = 'linear', ty = 'linear', xlim = c(0, 5.4), ylim = c(0, 5.4), ...)
    xline(xthresh, lty = 'dotdash', col = 'gray')
    lines(bb)
    lines(bb_infl, lwd = 3, col = 'red')
    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_scat
}

gate_live_stim = function(ff, show = FALSE, show.fn = NULL, ...) {
  param_ld = "DRAQ7 AF700"
  param_x = "CD4 BUV805"

  # find the top of the live blob
  bb = blob.boundary(ff, parameters = c(param_x, param_ld), location = bx(c(5000, 0)))
  thresh_ld = max(bb[, 2])

  ff_live = Subset(ff, rectangleGate("DRAQ7 AF700" = c(-Inf, thresh_ld)))

  if (show) {
    if (!is.null(show.fn)) {
      png(filename = show.fn, width = 600, height = 600)
    }
    pplot(ff, c(param_x, param_ld), xlim = c(-1, 5.4), ylim = c(-1, 5.4), ...)
    lines(bb, col = 'gray')
    yline(thresh_ld, lwd = 2, col = 'red')
    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_live
}

# this function may be specific to Wade's Macbook
init_reticulate = function() {
  # NOTRUN:
  # py_config()
  use_python("/usr/local/bin/python", required = TRUE)
}

do_optsne_reduction = function(centers, perplexity = 30, learning_rate = 341, show = FALSE) {
  require(ReductionWrappers)
  set.seed(137)   # so we'll get the same map for the same data
  res = optSNE(centers, n_components = 2, learning_rate = 341, verbose = 0)  # learning rate advised by running with default

  colnames(res) = c("t_sne_1", "t_sne_2")
  if (show) {
    plot(res, pch = 20, col = 'red')
    points(res)
  }

  res
}


################################################################################
################################################################################
# Following functions copied from EVICT/analysis_utils_tcell.R
################################################################################
################################################################################
# input is a flowFP (fp) and the corresponding flowSet (fs)
#  method = "median" returns median +- quartiles.
#  method = "mean" returns mean +- standard deviation
calculate_bin_phenotypes = function(fp, fs, method=c("median", "mean")) {
  parameters = parameters(fp)
  n_bins = 2 ^ nRecursions(fp)
  n_parameters = length(parameters)
  center = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(center) = parameters
  range = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(range) = parameters

  # for convenience, lump the frames in fs together and recalculate the fingerprint
  ff = as(fs, "flowFrame")
  fp = flowFP(fcs = ff, model = as(fp, "flowFPModel"))
  for (i in 1:n_bins) {
    idx = which(tags(fp)[[1]] == i)
    if (length(idx) == 0) {
      next
    }

    for (j in 1:n_parameters) {
      p = parameters[j]
      vals = exprs(ff)[idx, p]
      if (method == "median") {
        center[j, i] = median(vals, na.rm = TRUE)
        range[j, i] = (quantile(x = vals, probs = 0.75, na.rm = TRUE) -
                         quantile(x = vals, probs = 0.25, na.rm = TRUE)) / 2
      } else {

      }
    }
  }
  return(list(center = center, range = range))
}

# given a collection of bin centers, perform T-SNE dimensionality reduction
do_tsne_reduction = function(centers, perplexity = 30, show=FALSE) {
  require(Rtsne)
  set.seed(137)   # so we'll get the same map for the same data
  res = Rtsne(dist(centers), perplexity = perplexity)$Y
  colnames(res) = c("t_sne_1", "t_sne_2")
  if (show) {
    plot(res, pch = 20, col = 'red')
    points(res)
  }

  res

}

decorate_sample_panoply = function(tubes, mod, map, mfi, clst=NULL, superclus=NULL,
                                   colorscale=FALSE, superscale=FALSE, isLabeled=FALSE,
                                   cex=0.9, ...) {
  laymat = make_laymat(k = length(mfi), double = FALSE, allow_wedge = FALSE)
  layout(laymat)
  par(mar = c(0, 0, 0, 0) + 0.1)

  if (is.null(clst)) {
    # display color-coded bin populations
    count = map_sample_panoply(tubes, mod, map, cex = 1.5)
  } else {
    # make a cluster map
    count = NULL
    draw_cluster_map(map = map, clst = clst, superclus = superclus)
  }

  kde = bkde2D(map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))
  min_value = 0
  max_value = 5
  par(mar = c(0, 0, 2, 0) + 0.1)
  mfi_lut = 1:length(mfi)
  for (i in 1:length(mfi)) {
    pname = names(mfi)[mfi_lut[i]]
    contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE,
            xaxt = 'n', yaxt = 'n',
            col = 'darkgray', main = pname, cex.main = 2)

    mod_min_value = min_value
    cols = pcolor(mfi[[mfi_lut[i]]], min_value = mod_min_value, max_value = max_value)
    points(map, cex = 0.9 * cex)
    points(map, col = cols, pch = 20, cex = cex)
  }
  if(colorscale){draw_color_scale(min_col_value = min_value, max_col_value = max_value)}
  if(superscale){
    if(length(superclus)<=12){
      RColorBrewer::display.brewer.pal(length(superclus),"Set3")
      if(isLabeled==TRUE){
        text(labels=names(superclus),x=1:length(superclus),y=1,cex=1.5,srt=90)
      }else{
        text(labels=1:length(superclus),x=1:length(superclus),y=1,cex=2.5)
      }
    }else{
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      cl<-col_vector[1:length(superclus)]
      pie(rep(1,length(superclus)), col=cl)
    }


  }
  invisible(count)
}

cluster_map = function(map, h = NULL, k = NULL) {
  if (is.null(h) & is.null(k)) {
    stop("Must provide EITHER n or k\n")
  }
  require(cluster)
  ag = agnes(map)
  clst = cutree(as.hclust(ag), h = h, k = k)
  n_clust = max(clst)

  centers = matrix(NA, nrow = n_clust, ncol = ncol(map))
  boundaries = list()
  cindex = list()
  colnames(centers) = c("t_sne_1", "t_sne_2")
  for (i in 1:n_clust) {
    idx = which(clst == i)
    cindex[[i]] = idx
    centers[i, ] = c(median(map[idx, 1]), median(map[idx, 2]))
    boundaries[[i]] = get.hull(map[idx, ])
  }

  return(list(clst = clst, c_index = cindex, centers = centers, boundaries = boundaries))
}

make_laymat = function(k, double = FALSE, allow_wedge = FALSE) {
  #cat(paste("asking for k of ",k,"\n"))
  if (k == 4) {    # special case
    laymat = matrix(0, nrow = 3, ncol = 3)
    laymat[2, 2] = 1
    laymat[1, 1] = 2
    laymat[1, 3] = 3
    laymat[3, 1] = 4
    laymat[3, 3] = 5
    laymat[3, 2] = 6

    return(laymat)
  }
  if (double) {
    n = ceiling(k/8 + 2)
    if (n <= 5) {n = 6}
  } else {
    n = ceiling(k/4 + 1)
  }

  laymat = matrix(0, nrow = n, ncol = n)

  # make the central figure
  if (double) {
    for (i in 3:(n - 2)) {
      for (j in 3:(n - 2)) {
        laymat[i, j] = 1
      }
    }
  } else {
    for (i in 2:(n - 1)) {
      for (j in 2:(n - 1)) {
        laymat[i, j] = 1
      }
    }
  }

  # top
  if (double) {
    laymat[1, ] = (1:n)  + 1
    laymat[2, ] = ((n + 1):(2 * n)) + 1
  } else {
    laymat[1, ] = (1:n)  + 1
  }

  # middle
  if (double) {
    m = (2 * n + 1) + 1
    for (i in 3:(n - 2)) {
      for (j in c(1, 2, n - 1, n)) {
        laymat[i, j] = m
        m = m + 1
      }
    }
  } else {
    if(k==4){
      laymat[2, 2]<-4
      return(laymat)
    }else{
      m = (n + 1) + 1
      for (i in 2:(n - 1)) {
        for (j in c(1, n)) {
          laymat[i, j] = m
          m = m + 1
        }
      }
    }
  }

  # bottom
  if (double) {
    laymat[(n - 1), ] = m:(m + n - 1)
    laymat[n, ]       = (m + n):(m + (2 * n) - 1)
  } else {
    laymat[n, ] = m:(m + n - 1)
  }

  if (allow_wedge) {
    if (max(laymat) == k + 1) {
      cat("warning, not enough room for wedge\n")
    }
    laymat[laymat > k + 1] = 0
    laymat[n, n] = k + 2
  }

  laymat
}

# assumes biexp vert scale
draw_color_scale = function(min_col_value = 0, max_col_value = 5, ...) {
  ll = -0.5
  ul = bx(262143)

  vec = seq(ll, ul, length.out = 500)
  cols = pcolor(pvalue = vec, min_value = min_col_value, max_value = max_col_value)

  opar = par(mar = c(0, 15, 0, 0) + .1)
  plot(0, 0, pch = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "Fluorescence Intensity",
       xlim = c(0, 5), ylim = c(ll, ul), ...
  )
  for (i in 1:length(vec)) {
    y = vec[i]
    segments(x0 = 0, y0 = y, x1 = 1, y1 = y, col = cols[i], lwd = 3)
  }
  ax(axis = 2, instrument = 'diva', type = 'biexp', ...)
  par(opar)
}

draw_cluster_map = function(map, clst, dot_col='gray', superclus=NULL) {
  kde = bkde2D(map, bandwidth = c(2.5, 2.5), gridsize = c(501, 501))
  contour(kde$x1, kde$x2, kde$fhat, drawlabels = FALSE, col = 'darkgray', add = FALSE, xaxt = 'n', yaxt = 'n')
  if (!is.null(superclus)) {
    # color-code the superclusters
    #cl = hsv(h = seq(0, .6667, length.out = length(superclus)), s = 1, v = .75)
    if(length(superclus)<=12){
      cl <- RColorBrewer::brewer.pal(length(superclus),"Set3")
    }else{
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      cl<-col_vector[1:length(superclus)]
    }
    dot_col = rep('gray', length.out = nrow(map))
    for (i in 1:length(superclus)) {
      idx = bins_in_supercluster(clst = clst, sc = superclus[[i]])
      dot_col[idx] = cl[i]
    }
  }
  points(map, pch = 20, col = dot_col, cex = 2)
  for (i in 1:nrow(clst$centers)) {
    lines(clst$boundaries[[i]], col = 'darkgreen')
  }
  text(x = clst$centers[, 1], y = clst$centers[, 2], labels = 1:nrow(clst$centers), vfont = c("serif", "bold"), cex = 1.5)

}

# define a color function to transform a parameter value
pcolor = function(pvalue, min_value = 0, max_value = 4) {
  require(fields)
  top_col = 'red'
  bot_col = 'darkgreen'
  mid_col = 'yellow'
  zero_col = 'darkgray'
  zero_col = bot_col

  len = length(pvalue)
  if (length(which(pvalue < min_value)) == len) {
    col_values = rep(zero_col, len)
  } else if (length(which(pvalue > max_value)) == len) {
    col_values = rep(top_col, len)
  } else {
    pvalue[pvalue < min_value] = min_value
    pvalue[pvalue > max_value] = max_value

    col_values = color.scale(
      z = pvalue, col = two.colors(
        n = 100,
        start = bot_col,
        end = top_col,
        middle = mid_col),
      zlim = c(min_value, max_value)
    )
    col_values[which(pvalue <= min_value)] = zero_col
  }


  col_values
}

# do some parallel coordinate plots
parallel_pheno = function(mfi, idx_bin = 1:length(mfi[[1]]), parameters = names(mfi), col = 'red', mfi_colors = FALSE, show_contours = TRUE, bars = TRUE, ...) {

  if (is.numeric(parameters)) {
    pnames = colnames(ff)[parameters]
  } else {
    pnames = parameters
  }
  if (show_contours) {
    cont_col = 'darkgray'
  } else {
    cont_col = 'white'
  }
  mfi = data.frame(mfi)
  # hyphens in parameter names get turned into dots here.  Change back
  colnames(mfi) = sub(pattern = ".", replacement = "-", x = colnames(mfi), fixed = TRUE)

  parallel_contours(mfi, parameters = parameters, col = cont_col, ...)
  mfi = mfi[idx_bin, parameters]

  med_vec = vector(mode = 'numeric')
  q1_vec = vector(mode = 'numeric')
  q3_vec = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    tmp = fivenum(mfi[, i])
    med_vec[i] = tmp[3]
    q1_vec[i] = tmp[2]
    q3_vec[i] = tmp[4]
  }
  # draw the median
  if (bars) {
    if (mfi_colors) {
      col = pcolor(med_vec, min_value = 0, max_value = 5)
    }
    add_bars(vals = med_vec, yvals = 1:length(parameters), col = col)
  } else {
    x = med_vec
    y = 1:length(parameters)
    lines(x, y, col = col, lwd = 3)
  }

  # draw the flags
  for (i in 1:length(parameters)) {
    if (bars) {
      draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = NA, cex = 2, lwd = 2)
    } else {
      draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = med_vec[i], cex = 2, lwd = 2)
    }
  }
}

parallel_contours = function(mfi, parameters = names(mfi), col = 'blue', ...) {
  if (is.numeric(parameters)) {
    pnames = colnames(ff)[parameters]
  } else {
    pnames = parameters
  }
  plot(0, 0, pch = '', xlim = c(0, bx(262143)), ylim = c(1 - .3, length(parameters) + .3),
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '')
  ax(1, instrument = 'diva', type = 'biexp')
  axis(side = 2, labels = pnames, at = 1:length(pnames), las = 1, ...)

  x = matrix(NA, nrow = nrow(mfi) * length(parameters), ncol = 2)
  k = 1
  for (i in 1:nrow(mfi)) {
    for (p in 1:length(parameters)) {
      x[k, 1] = mfi[i, p]
      x[k, 2] = p
      k = k + 1
    }
  }
  kde = bkde2D(x = x, bandwidth = c(.1, 1), gridsize = c(501, 501))
  kde$fhat = kde$fhat / max(kde$fhat)
  contour(x = kde$x1, y = kde$x2, z = kde$fhat, col = col,
          xaxt = 'n', yaxt = 'n',
          drawlabels = FALSE,
          # levels = seq(.01, .2, length.out = 20),
          add = TRUE)

}

draw_y_grid = function(lo = -4, hi = 0) {
  require(fields)
  minor = seq(2, 9, by = 1)
  major = 10 ^ seq(lo, hi, by = 1)
  for (maj in major) {
    yline(maj, col = "darkgray", lwd = 2)
    for (mn in minor) {
      yline(maj * mn, col = 'lightgray', lwd = 1)
    }
  }
}

add_bars = function(vals, yvals, col) {
  hw = 0.4
  if (length(col) == 1) {
    col = rep(col, length(vals))
  }
  for (i in 1:length(vals)) {
    rect(xleft = 0, ybottom = yvals[i] - hw, xright = vals[i], ytop = yvals[i] + hw, col = col[i], border = 'black')
  }
}

draw_flag = function(y, q1, q3, med = NA, ...) {
  segments(x0 = q1, y0 = y, x1 = q3, y1 = y, ...)
  if (!is.na(med)) {points(med, y, pch = 20, ...)}
}

################################################################################
################################################################################
# END: functions copied from EVICT/analysis_utils_tcell.R
################################################################################
################################################################################

# ff is the flowFrame from which the fingerprint fp was computed.
# clst contains the assignment of fp bins to clusters
tag_event_cluster = function(ff, fp, clst) {
  ctag = vector('numeric', length = nrow(ff))   # per-event assigment to clusters
  tg = tags(fp)[[1]]
  for (i in 1:length(clst$c_index)) {
    cbins = clst$c_index[[i]]
    idx = which(tg %in% cbins)
    ctag[idx] = i
  }

  ctag
}

################################################################################
################################################################################
# here starts igraph stuff
################################################################################
################################################################################
build_graph = function(mfi) {
  require(igraph)
  # make mfi into a distance matrix
  mat = matrix(NA, nrow = length(mfi[[1]]), ncol = 0)

  for (i in 1:length(mfi)) {
    mat = cbind(mat, mfi[[i]])
  }
  colnames(mat) = names(mfi)
  dst = as.matrix(stats::dist(mat))
  g = graph_from_adjacency_matrix(adjmatrix = dst, mode = 'undirected', weighted = TRUE)

  # some algorithms interpret weight as distance, and others as strength.  Let's
  # preserve both interpretations so they can be conveniently swapped.  Initial weight
  # will be distance.  This is appropriate for MST.
  # modularity treats weights as strengths
  # betweenness treats weights as distances
  dstnce = E(g)$weight
  strngt = 1 / dstnce
  edge_attr(g, 'distance') <- dstnce
  edge_attr(g, 'strength') <- strngt

  g
}

# create a new graph from a communities object.  g is the MST graph from which the
# communities object 'comm' was induced.
# note:  marker expressions should already have been added to g.
make_graph_from_community = function(comm, g) {
  markers = vertex_attr_names(g)
  markers = markers[-which(markers == "name")]

  n_vertices = max(comm$membership)
  gcomm = make_empty_graph(n = n_vertices, directed = FALSE)
  sizes = vector('numeric')
  indices = list()
  for (i in 1:n_vertices) {
    indices[[i]] = idx = which(comm$membership == i)
    V(gcomm)[i]$size = length(indices[[i]])
    for (j in 1:length(markers)) {
      vertex_attr(gcomm, markers[j])[i] <- mean(vertex_attr(g, markers[j])[idx])
    }
  }

  # add edges.  Edges are the cummulative strengths between communities.
  k = 1
  for (i in 1:(n_vertices - 1)) {
    for (j in (i + 1):n_vertices) {
      res = edges_between_communities(comm, g, indices[[i]], indices[[j]])
      if(res$strength != 0) {
        # add an edge
        gcomm = add_edges(gcomm, c(i, j))
        E(gcomm)$strength[k] = res$strength
        k = k + 1
        cat(i, j, "\n")
      }
    }
  }
  E(gcomm)$distance = 1/E(gcomm)$strength
  E(gcomm)$weight = E(gcomm)$strength
  mst(gcomm)
}

# c1 is the vector of indices of community 1, c2 of community 2
edges_between_communities = function(comm, g, c1, c2) {
  edges = E(g)[c1 %--% c2]
  strength = sum(E(g)$strength[edges])

  return(list(edges = edges, strength = strength))
}

# if the graph is induced from a distance matrix, then the edge weights are
# equal to distances.  However, the strongest edges should be the ones closest
# in metric distance, so invert the weights in this case.
set_weight_as = function(g, weight = c("distance", "strength")) {
  weight = match.arg(weight)
  E(g)$weight = edge_attr(g, weight)

  g
}

# best so far
attach_layout_fr = function(g, overwrite = TRUE) {
  set.seed(137)
  g = add_layout_(g, with_fr(), overwrite = overwrite)

  g
}

# also nice
attach_layout_lgl = function(g, overwrite = TRUE) {
  set.seed(137)
  g = add_layout_(g, with_lgl(), overwrite = overwrite)

  g
}

attach_layout = function(g, layfun, overwrite = TRUE) {
  set.seed(137)
  g = add_layout_(g, layfun, overwrite = overwrite)

  g
}


plot_community_graph = function(g, marker, vs = 3, ms = 0, log.size = TRUE, vertex.frame = TRUE, cex.main) {
  g = set_weight_as(g, "distance")

  if (log.size) {
    vsize = vs * log10(V(g)$size)
    vsize[vsize < ms] = ms
  } else {
    vsize = V(g)$size
    mx = max(vsize)
    mn = min(vsize)
    vsize = vs * vsize / mx
    if (ms != 0) {
      vsize[vsize < ms] = ms
    }
  }
  if (vertex.frame) {
    vfc = 'black'
  } else {
    vfc = NA
  }
  plot(g, vertex.size = vsize, vertex.color = pcolor(vertex_attr(gcomm, marker)),
       vertex.frame.color = vfc, vertex.label = NA)
  title(main = marker, cex.main = cex.main)
}


plot_comm_spread = function(g, markers = NULL, vs = 3, ms = 1, log.size = TRUE, vertex.frame = FALSE, cex.main) {
  if (is.null(markers)) {
    markers = vertex_attr_names(g)
    markers = markers[-which(markers == "size")]
  }

  # calculate plot layout
  n = length(markers)
  sq = sqrt(n)
  frac = sq - floor(sq)
  if(frac == 0) {
    ac = dn = floor(sq)
  } else {
    ac = floor(sq) + 1
    dn = ceiling(n / ac)
  }

  par(mfrow = c(dn, ac), mar = c(0, 0, 2, 0))
  for (i in 1:length(markers)) {
    marker = markers[i]
    plot_community_graph(g, marker, vs = vs, ms = ms, log.size = log.size, vertex.frame = vertex.frame, cex.main)
  }
}

################################################################################
################################################################################
# conventional clustering analyses
################################################################################
################################################################################

# ag is an agnes object
# look for a change in slope of nclust vs cut height
# breaks is the number of height breaks
advise_n_clust = function(ag, breaks = 500, show = TRUE, ...) {
  obj = as.hclust(ag)
  ht = seq(min(ag$height), max(ag$height), length.out = breaks)
  nclust = vector('numeric')
  for (i in 1:breaks) {
    nclust[i] = max(cutree(obj, h = ht[i]))
  }
  df = data.frame(nclust = nclust, height = ht)
  df = df[breaks:1, ]

  # eliminate redundant heights
  nclust = unique(df$nclust)
  ht = vector('numeric')
  for (i in 1:length(nclust)) {
    ht[i] = min(df$height[df$nclust == nclust[i]])
  }
  df2 = data.frame(nclust = nclust, height = ht)


  # fit the asymptotic behavior, and look for where it starts to diverge
  ignore_last = 3
  npts = 20

  idx2 = (nrow(df2) - npts - ignore_last):(nrow(df2) - ignore_last)
  ln2 = predict(lm(height ~ nclust, data = df2, subset = idx2), newdata = df2)
  dfh = data.frame(x = nclust, y = ln2)

  resid = df2$height - dfh$y

  # where does residual exceed more than x% of max?
  crit1 = 0.10
  crit2 = 0.05
  idx1 = min(which(resid < crit1 * max(resid)))
  nclust1 = df2$nclust[idx1]
  idx2 = min(which(resid < crit2 * max(resid)))
  nclust2 = df2$nclust[idx2]
  if (show) {
    plot(df2, ...)
    lines(dfh$x, dfh$y, col = 'blue')

    segments(x0 = df2$nclust[idx1], y0 = 0, x1 = df2$nclust[idx1], y1 = df2$height[idx1])
    segments(x0 = df2$nclust[idx2], y0 = 0, x1 = df2$nclust[idx2], y1 = df2$height[idx2])

    points(df2$nclust[idx1], df2$height[idx1], pch = 20, col = 'red', cex = 2)
    points(df2$nclust[idx2], df2$height[idx2], pch = 20, col = 'red', cex = 2)
    points(df2$nclust[idx1], df2$height[idx1], cex = 2)
    points(df2$nclust[idx2], df2$height[idx2], cex = 2)

    x = 10
    y = 0.75 * max(df2$height)
    text(x = x, y = y, labels = sprintf("Recommend between %d and %d clusters", nclust1, nclust2), pos = 4, cex = 1.5)
  }
  return(list(n1 = nclust1, n2 = nclust2))
}

# calculate the within-cluster sum of squares
# https://discuss.analyticsvidhya.com/t/what-is-within-cluster-sum-of-squares-by-cluster-in-k-means/2706
calculate_wss = function(clustering, max_clust = 100) {
  mat = clustering$data
  ndim = ncol(mat)              # number of dimensions

  wss = vector('numeric')
  for (k in 1:max_clust) {      # for all clusterings
    idx = cutree(as.hclust(clustering), k = k)
    d4 = 0
    for (i in 1:k) {       # for each cluster
      pts = which(idx == i)
      npts = length(pts)
      d3 = 0
      for (j in 1:ndim) {         # for each dimension of data
        xbar = mean(mat[pts, j])
        d2 = sum((mat[pts, j] - xbar) ^ 2)
        d3 = d3 + d2
      }
      d4 = d4 + d3
    }
    wss[k] = d4
  }

  wss
}

# estimate the first derivative using the symmetric difference quotient
lslope = function (kde, normalize=TRUE) {
  x = kde[,1]
  y = kde[,2]
  npts = length(y)
  yp = vector('numeric')
  for (i in 2:(npts - 1)) {
    yp[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
  }
  yp[1] = yp[2]
  yp[npts] = yp[npts - 1]
  yp[yp == -Inf] = 0.0

  if (normalize) {
    yp = yp / max (abs(yp))
  }
  res = list(x=x, y=yp)
  res
}

#construct a community object based on hierarchical (e.g. agnes) clustering
agnes_to_community = function(ag, nclust) {
  comm = list()
  vcount = nrow(ag$data)
  comm$vcount = vcount
  comm$names = as.character(1:vcount)

  membership = cutree(as.hclust(ag), k = nclust)
  comm$membership = membership
  class(comm) = "communties"

  comm
}

# we noted that write.FCS mucks up the range, minRange and maxRange values
# let's see if we can fix that.
# This function only works for YO data, in that
#    1 - it assumes 4 scattering parameters that were linearly transformed
#    2 - 16 FL parameters that were biexponentially transformed
#    3 - one Time parameter (haven't investigated this)
recover_FCS = function(fn) {
  rng = c(rep(262144.0, length = 22), 13304.6)

  minr = rep(-0.17156, length = 23)
  minr[1:6] = 0.000
  minr[c(1, 4)] = -0.0023
  minr[23] = 0.000

  maxr = rep(5.4185, length = 23)
  maxr[1:6] = 5.4000
  maxr[23] = rng[23]

  pmat = matrix(c(rng, minr, maxr), ncol = 3)
  colnames(pmat) = c("range", "minRange", "maxRange")

  ff = read.FCS(fn)
  pData(parameters(ff))[, 3:5] = pmat

  ff
}


