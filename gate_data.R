#
# gate_data.R
#
# Apply the previously optimized gating strategy, and write files out for downstream
# processing.
#
#  2019-12-07  WTR
#

library(wadeTools)
library(flowFP)
library(fields)
source("~/git/R/clustering/yo_utils.R")

# proj_base will depend on the particular analysis platform
proj_base = "~/Data/clustering/"
data_base = tight(proj_base, "FR-FCM-ZZGS/")
gated_base = tight(proj_base, "gated_fcs/")

files = grep(pattern = "fcs", x = dir(data_base), fixed = TRUE, value = TRUE)

# for (i in 1:length(files)) {
for (i in 1:5) {
    cat("working on", i, "...")
  fn = files[i]
  fbase = sub(pattern = ".fcs", replacement = "", x = fn, fixed = TRUE)
  ff = get_sample(tight(data_base, fn))

  # do the gating EXACTLY like was done in survey_data.R
  ff = gate_clean(ff, show = FALSE)
  ff = gate_singlet(ff, show = FALSE)
  ff = gate_live(ff, show = FALSE)

  # write the gated data
  # omit index numbers so that wi=2
  # e can use the CSV manifest file for reading downstream
  outfile = tight(gated_base, fbase, ".fcs")
  write.FCS(ff, filename = outfile)
  cat("done.\n")
}


