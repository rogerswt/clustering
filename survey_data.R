#
# survey_data.R
#
# survey the dataset.
#
# 2019-11-07  WTR

source("~/git/R/tools/sourceTools.R")
library(flowFP)

source("~/git/R/clustering/yo_utils.R")

#  Change the following locations to correspond with where you've downloaded
#  the Young/Old dataset from FlowRepository.  I recommend that you
#     a) create a directory somewhere, called "clustering"
#     b) download the data from FlowRepository into that directory,
#     c) create a subdirectory called "results" to collect the
#        figures that you make.
#
proj_base = "~/Data/clustering/"
data_base = tight(proj_base, "FR-FCM-ZZGS/")
pic_base = tight(proj_base, "results/")

files = grep(pattern = "fcs", x = dir(data_base), fixed = TRUE, value = TRUE)

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
  png(filename = tight(pic_base, "raw_", i, "_", fbase, ".png"), width = 1000, height = 600)
  par(mfrow = c(1, 2))
  pplot(ff, c("FSC-A", "SSC-A"))
  pplot(ff, c("CD3Q605", "SSC-A"))
  dev.off()
  cat("done.\n")
}
fs = flowSet(ff_list)

# Use flowFP to look for cases where the primary gating markers differ
# significantly from the average
fp_params = c("SSC-A","CD3Q605","LIVEDEAD")
mod = flowFPModel(fcs = fs, parameters = fp_params, nRecursions = 6)
fp  = flowFP(fcs = fs, model = mod)
flowFP::plot(fp, type = "qc", red_limit = 3)
dev.print(png, tight(pic_base, "qc_raw.png"), width = 800, height = 800)

# based on this we see that instance #15 is wacky.  We will censor this instance.

# now let's do a trial gating, after which we will repeat the QC step to look
# for any other staining irregularities
ff_list = list()
for (i in 1:length(files)) {
  cat("working on", i, "...")
  fn = files[i]
  fbase = sub(pattern = ".fcs", replacement = "", x = fn, fixed = TRUE)
  ff = get_sample(tight(data_base, fn))

  fn_clean = tight(pic_base, "clean_", i, "_", fbase, ".png")
  fn_singlet = tight(pic_base, "singlet_", i, "_", fbase, ".png")
  fn_live = tight(pic_base, "live_", i, "_", fbase, ".png")
  fn_scat = tight(pic_base, "scat_", i, "_", fbase, ".png")
  ff = gate_clean(ff, show = TRUE, show.fn = fn_clean)
  ff = gate_singlet(ff, show = TRUE, show.fn = fn_singlet)
  ff = gate_live(ff, show = TRUE, show.fn = fn_live)
  ff = gate_scat(ff, show = TRUE, show.fn = fn_scat)

  ff_list[[i]] = Subset(ff, sampleFilter(10000))
  cat("done.\n")
}

# we are re-using fs, fp, mod from before
fs = flowSet(ff_list)
fp_params = colnames(ff)[c(7:9, 11, 13:22)]  # all parameters not used for gating
mod = flowFPModel(fcs = fs, parameters = fp_params, nRecursions = 6)
fp  = flowFP(fcs = fs, model = mod)
flowFP::plot(fp, type = "qc", red_limit = 3)
dev.print(png, tight(pic_base, "qc_gated.png"), width = 800, height = 800)


