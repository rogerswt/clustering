#
# gate_data.R
#
# Apply the previously optimized gating strategy, and write files out for downstream
# processing.
#
#  2019-12-07  WTR
#

source("~/git/R/tools/sourceTools.R")
library(flowFP)

source("~/git/R/clustering/yo_utils.R")

#  As in suvey_data.R, change the following locations to correspond with where you've downloaded
#  the Young/Old dataset from FlowRepository.  I recommend that you
#     a) create a directory somewhere, called "clustering"
#     b) download the data from FlowRepository into that directory,
#     c) create a subdirectory called "gated_fcs" to collect the
#        gated data.
#
proj_base = "~/Data/clustering/"
data_base = tight(proj_base, "FR-FCM-ZZGS/")
gated_base = tight(proj_base, "gated_fcs/")

files = grep(pattern = "fcs", x = dir(data_base), fixed = TRUE, value = TRUE)

for (i in 1:length(files)) {
  cat("working on", i, "...")
  fn = files[i]
  fbase = sub(pattern = ".fcs", replacement = "", x = fn, fixed = TRUE)
  ff = get_sample(tight(data_base, fn))

  # do the gating EXACTLY like was done in survey_data.R
  ff = gate_clean(ff, show = FALSE)
  ff = gate_singlet(ff, show = FALSE)
  ff = gate_live(ff, show = FALSE)
  ff = gate_scat(ff, show = FALSE)

  # write the gated data
  outfile = tight(gated_base, "gated_", i, "_", fbase, ".fcs")
  write.FCS(ff, filename = outfile)
  cat("done.\n")
}


