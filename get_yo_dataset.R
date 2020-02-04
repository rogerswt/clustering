#
# get_yo_dataset.R
#
# Download the Young/Old dataset from FlowRepository
#
# 2020-02-04  WTR
#

library(FlowRepositoryR)

# choose a directory to contain the downloaded data.  For example:
data_dir = "~/Data/clustering/"

ds_id = "FR-FCM-ZZGS"
ds = flowRep.get(ds_id)
ds = download(ds)

# print a summary of the download
summary(ds)

