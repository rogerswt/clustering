#
# figs_for_clustering_ppt.R
#
# Generate figures for the clustering lecture
#
# 2020-01-17  WTR
#

library(wadeTools)
library(cluster)
library(RColorBrewer)
library(dbscan)
library(mclust)

cols = brewer.pal(9, "Set1")
# cols = brewer.pal(20, name = "Spectral")

pic_dir = "~/Data/clustering/figs_clustering/"

# retrieve and repair a gated file from the YO dataset
ff = read.FCS("~/Data/clustering/gated_fcs/Tphe09943-005-00_F4_R.fcs")
ff = replace_phedata_ranges(ff)

pplot(ff, c("CD3Q605", "CD4PETR"), xlim = c(-1, 5.4), ylim = c(-1, 5.4), showZero = TRUE)
dev.print(png, tight(pic_dir, "pplot_cd3_cd4.png"), width = 600, height = 600)

dat1 = exprs(ff)[, c("CD3Q605", "CD4PETR")]

set.seed(137)
dat2 = dat1[sample(1:nrow(dat1), size = 100000), ]    # reduced to 100k events
dat3 = dat1[sample(1:nrow(dat1), size = 10000), ]     # reduced to 10k events
dat4 = dat1[sample(1:nrow(dat1), size = 5000), ]      # reduced to 5k events
dat5 = dat1[sample(1:nrow(dat1), size = 1000), ]      # reduced to 1k events

# kmeans
opar = par(mfrow = c(2,2))
cl_km = kmeans(dat1, centers = 5)
plot(dat1, pch = '.', col = cols[cl_km$cluster], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

cl_km = kmeans(dat1, centers = 5)
plot(dat1, pch = '.', col = cols[cl_km$cluster], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

cl_km = kmeans(dat1, centers = 5)
plot(dat1, pch = '.', col = cols[cl_km$cluster], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

cl_km = kmeans(dat1, centers = 5)
plot(dat1, pch = '.', col = cols[cl_km$cluster], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

dev.print(png, tight(pic_dir, "kmeans_1.png"), width = 1000, height = 1000)
par(opar)

library(cluster)
# pam

opar = par(mfrow = c(2,2))

pam = pam(dat4, k = 5)
plot(dat4, pch = 20, cex = .3, col = cols[pam$clustering], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

pam = pam(dat4, k = 5)
plot(dat4, pch = 20, cex = .3, col = cols[pam$clustering], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

pam = pam(dat4, k = 5)
plot(dat4, pch = 20, cex = .3, col = cols[pam$clustering], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

pam = pam(dat4, k = 5)
plot(dat4, pch = 20, cex = .3, col = cols[pam$clustering], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

dev.print(png, tight(pic_dir, "pam_1.png"), width = 1000, height = 1000)
par(opar)


# agnes
#
opar = par(mfrow = c(1, 2))
ag = agnes(dat5)
res = cutree(as.hclust(ag), k = 5)
pltree(ag, labels = rep("", length = nrow(dat5)))
yline(1.125, lty = 'dotdash', col = 'red')

plot(dat5, pch = 20, cex = .5, col = cols[res], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

dev.print(png, tight(pic_dir, "agnes_1.png"), width = 1000, height = 600)
par(opar)

opar = par(mfrow = c(1, 2))
ag = agnes(dat5)
res = cutree(as.hclust(ag), h = 0.8)
pltree(ag, labels = rep("", length = nrow(dat5)))
yline(0.8, lty = 'dotdash', col = 'red')

plot(dat5, pch = 20, cex = .5, col = cols[res], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

dev.print(png, tight(pic_dir, "agnes_2.png"), width = 1000, height = 600)
par(opar)

# fanny
opar = par(mfrow = c(1, 1))
fan = fanny(dat4, k = 5)  # 5k events, 5 clusters
plot(dat4, pch = 20, cex = .5, col = cols[fan$clustering], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')
dev.print(png, tight(pic_dir, "fanny_1.png"), width = 600, height = 600)

par(opar)

# density-based
opar = par(mfrow = c(1, 1))
res = dbscan(dat3, eps = .05, minPts = 20)

plot(dat3, pch = 20, cex = .2, col = cols[res$cluster + 1], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')
dev.print(png, tight(pic_dir, "dbscan_1.png"), width = 600, height = 600)

par(opar)

# mclust
opar = par(mfrow = c(1, 1))

mod1 = Mclust(dat3)
plot(mod1, what = "classification", xlim = c(-1,5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n')
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

dev.print(png, tight(pic_dir, "mclust_1.png"), width = 600, height = 600)

par(opar)


################################################################################
################################################################################
# construct a synthetic dataset with lots of gaussian clusters,
# in multiple dimensions
################################################################################
################################################################################
library(MASS)
# create a multidimensional uniformly random covariance matrix
make_covar = function(ndim, minv = .075, maxv = .4) {
  lsig = runif(n = ndim, min = minv, max = maxv)
  p = qr.Q(qr(matrix(rnorm(ndim ^ 2), ndim)))
  Sigma = crossprod(p, p * lsig)

  Sigma
}

set.seed(137)
ndim = 10
nclust = 20
max_per_cluster = 500
min_per_cluster = 50
centers = matrix(NA, nrow = nclust, ncol = ndim)
for (i in 1:nclust) {
  centers[i, ] = runif(n = ndim, min = 0, max = 5)
}
dat = matrix(NA, ncol = ndim, nrow = 0)
for (i in 1:nclust) {
  npts = sample(min_per_cluster:max_per_cluster, size = 1)
  tmp = mvrnorm(n = npts, mu = rep(centers[i, ]), Sigma = make_covar(ndim = ndim))
  dat = rbind(dat, tmp)
}

res = Mclust(dat, G = 5:30)

# illustrate different distance measures
library(factoextra)

opar = par(mfrow = c(2, 2))
ag1 = agnes(dat4)
res1 = cutree(as.hclust(ag1), k = 10)

plot(dat4, pch = 20, cex = .5, col = cols[res1], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n', main = "Euclidean")
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

ag2 = agnes(get_dist(dat4, method = "manhattan"), diss = TRUE)
res2 = cutree(as.hclust(ag2), k = 10)

plot(dat4, pch = 20, cex = .5, col = cols[res2], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n', main = "Manhattan")
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

ag3 = agnes(get_dist(dat4, method = "pearson"), diss = TRUE)
res3 = cutree(as.hclust(ag3), k = 10)

plot(dat4, pch = 20, cex = .5, col = cols[res3], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n', main = "Pearson")
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

ag4 = agnes(get_dist(dat4, method = "minkowski"), diss = TRUE)
res4 = cutree(as.hclust(ag4), k = 10)

plot(dat4, pch = 20, cex = .5, col = cols[res4], xlim = c(-1, 5.4), ylim = c(-1, 5.4), xaxt = 'n', yaxt = 'n', main = "Minkowski")
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'biexp')

dev.print(png, tight(pic_dir, "distance_metrics.png"), width = 1000, height = 1000)

par(opar)



