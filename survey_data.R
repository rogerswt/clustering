#
# survey_data.R
#
# survey the dataset.
#
# 2019-11-07  WTR

library(wadeTools)

source("~/git/R/tailor/utils/yo_utils.R")

# proj_base will depend on the particular analysis platform
proj_base = "~/Data/Independent_Consulting/Penn/Matei/"
data_base = tight(proj_base, "data/young_old/FR-FCM-ZZGS/")
pic_base = tight(proj_base, "results/young_old/")

files = grep(pattern = "fcs", x = dir(data_base), fixed = TRUE, value = TRUE)


################################################################################
################################################################################
# This first section makes figures for perusing, and makes a sampled flowSet
################################################################################
################################################################################

ff_list = list()
for (i in 1:length(files)) {
  cat("working on", i, "...")
  fn = files[i]
  fbase = sub(pattern = ".fcs", replacement = "", x = fn, fixed = TRUE)
  ff = get_sample(tight(data_base, fn))

  # create a flowSet of subsampled data to look for anomalous distributions
  ff_list[[i]] = Subset(ff, sampleFilter(10000))

  # create some bivariate figures to visually look for outliers
  # NOTE: we are writing directly to the png file, since rendering in the RStudio
  # plot window can be quite slow
  fn = sprintf("%s%03d%s", tight(pic_base, "raw_"), i, tight("_", fbase, ".png"))
  png(filename = fn, width = 1000, height = 600)
  par(mfrow = c(1, 2))
  pplot(ff, c("FSC-A", "SSC-A"))
  pplot(ff, c("CD3Q605", "SSC-A"))
  dev.off()
  cat("done.\n")
}
fs_raw = flowSet(ff_list)

################################################################################
################################################################################
# This next section uses flowFP to look for outliers
################################################################################
################################################################################
# Use flowFP to look for cases where the primary gating markers differ
# significantly from the average
fp_params = c("SSC-A","CD3Q605","LIVEDEAD")
mod = flowFPModel(fcs = fs_raw, parameters = fp_params, nRecursions = 6)
fp  = flowFP(fcs = fs_raw, model = mod)
flowFP::plot(fp, type = "qc", red_limit = 3)
dev.print(png, tight(pic_base, "qc_raw.png"), width = 800, height = 800)

# based on this we see that instance #15 is wacky.  We will censor this instance.


################################################################################
################################################################################
# This next section is a test of gating, and also makes a gated sampled flowSet
################################################################################
################################################################################

# now let's do a trial gating, after which we will repeat the QC step to look
# for any other staining irregularities
ff_list = list()
for (i in 1:length(files)) {
  cat("working on", i, "...")
  fn = files[i]
  fbase = sub(pattern = ".fcs", replacement = "", x = fn, fixed = TRUE)
  ff = get_sample(tight(data_base, fn))

  fn_clean = sprintf("%s%03d%s", tight(pic_base, "clean_"), i, tight("_", fbase, ".png"))
  fn_singlet = sprintf("%s%03d%s", tight(pic_base, "singlet_"), i, tight("_", fbase, ".png"))
  fn_live = sprintf("%s%03d%s", tight(pic_base, "live_"), i, tight("_", fbase, ".png"))
  ff = gate_clean(ff, show = TRUE, show.fn = fn_clean)
  ff = gate_singlet(ff, show = TRUE, show.fn = fn_singlet)
  ff = gate_live(ff, show = TRUE, show.fn = fn_live)

  # make some bivariates
  fn_biv = sprintf("%s%03d%s", tight(pic_base, "gated_"), i, tight("_", fbase, ".png"))
  png(filename = fn_biv, width = 1000, height = 1000)
    par(mfrow = c(2, 2))
    pplot(ff, c("CD4PETR", "CD8Q705"), xlim = c(-1, 5.4), ylim = c(-1, 5.4))
    pplot(ff, c("CD45RAQ655", "CCR7PECY7"), xlim = c(-1, 5.4), ylim = c(-1, 5.4))
    pplot(ff, c("PD1PE", "CD28PECY5"), xlim = c(-1, 5.4), ylim = c(-1, 5.4))
    pplot(ff, c("CD95PERCPCY55", "CD11BAPCCY7"), xlim = c(-1, 5.4), ylim = c(-1, 5.4))
    par(mfrow = c(1, 1))
  dev.off()

  ff_list[[i]] = Subset(ff, sampleFilter(10000))
  cat("done.\n")
}

fs_gated = flowSet(ff_list)
fp_params = colnames(ff)[c(7:9, 11:22)]  # all except LIVEDEAD
mod = flowFPModel(fcs = fs_gated, parameters = fp_params, nRecursions = 6)
fp  = flowFP(fcs = fs_gated, model = mod)
flowFP::plot(fp, type = "qc", red_limit = 3)
dev.print(png, tight(pic_base, "qc_gated.png"), width = 800, height = 800)

# persist the sampled flowSets for possible further examination
save(fs_raw, fs_gated, file = tight(pic_base, "sampled_flowsets.rda"))

################################################################################
################################################################################
# This bit calculates and sorts the qc value so we can poke at it visually
################################################################################
################################################################################

# calculate qcval as in the qc figure
fpMatrix = counts(fp, transformation='log2norm')
fpMatrix[which(is.infinite(fpMatrix))] = NA
qcval = apply(fpMatrix, 1, na.rm=TRUE, sd)

idx = sort(qcval, decreasing = TRUE, index.return = TRUE)$ix

################################################################################
################################################################################
# This next section makes univariate distributions
################################################################################
################################################################################

# take a look at univariate distributions
agg = suppressWarnings(as(fs_gated[-15], "flowFrame"))   # remove instance 15
exprs(agg) = exprs(agg)[,which(colnames(agg) != "Original")]
par(mfrow = c(4,4), mar = c(2,2,2,0))
for (p in fp_params) {
  kde = bkde(exprs(agg)[,p], bandwidth = .01, gridsize = 1001)
  kde$y = kde$y / max(kde$y)
  plot(kde, type = 'l', xlim = c(-1, 5.4), ylim = c(0, 1),
       xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
       main = p)
  fac = 10.0
  lines(kde$x, fac * kde$y, col = 'red')
  xline(0, lty = 'dotdash')
  ax(axis = 1, instrument = 'diva', type = 'biexp')
}
dev.print(png, tight(pic_base, "univariate_gated.png"), width = 600, height = 600)

# add a comment here
