#
# yo_utils.R
#
# Utility functions for analysis of the young/old dataset
#
# 2019-11-07  WTR

# retrieve and prepare files as specified by a list of files and a dbase
get_sample = function(fn, compensate=TRUE, transform=TRUE, derail=TRUE, nice.names = TRUE, verbose=FALSE) {

  ff = read.FCS(fn)
  fl_params = which(colnames(ff) %in% colnames(keyword(ff)$SPILL))
  sc_params = 1:(fl_params[1] - 1)

  if (compensate) {ff = autocomp(ff)}

  if (derail) {
    ff = Subset(ff, rectangleGate("FSC-A"=c(-Inf, 262142), "SSC-A"=c(-Inf, 262142)))
  }
  if (transform) {
    ff = doTransform(ff, cols = sc_params, method = 'linear')
    ff = doTransform(ff, cols = fl_params, method = 'biexp')
  }
  if (nice.names) {
    dnames = colnames(ff)    # detector names

    names = parameters(ff)$desc
    names[sc_params] = c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")
    scat_names = names[sc_params]
    fl_names = names[fl_params]
    parameters(ff)$desc[sc_params] = scat_names

    colnames(ff) = c(scat_names, fl_names, "Time")
    parameters(ff)$desc = c(scat_names, paste(dnames[fl_params], fl_names), "Time")
  }

  ff
}

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
  }
  ff = Subset(res, rectangleGate("clean" = c(0.5, Inf)))
  exprs(ff) = exprs(ff)[, -which(colnames(ff) == "clean")]

  ff
}

gate_singlet = function(ff, show = FALSE, show.fn = NULL) {
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
    pplot(ff, p_singlet, tx = 'linear', ty = 'linear', xlim =c(0, 5), ylim = c(0, 5))
    lines(bb_singlet, col = 'red', lwd = 2)
    xline(fsc_thresh, lty = 'dotdash')
    yline(ssc_thresh, lty = 'dotdash')
    text(0.5, 4.5, labels = sprintf("%.2f%% Singlets", 100 * n_singlet / n_orig), pos = 4)
    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_singlet
}

# look for the LIVEDEAD negative, CD3 positive blob
gate_live = function(ff, show = FALSE, show.fn = NULL) {
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
    pplot(ff, plist = params, xlim = c(0, 5.4), ylim = c(0, 5.4))
    lines(bb)
    lines(bb_infl, lwd = 3, col = 'red')
    xline(pre_cd3, lty = 'dotdash')
    yline(pre_live, lty = 'dotdash')
    if (!is.null(show.fn)) {
      dev.off()
    }
  }

  ff_live
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



