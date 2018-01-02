# Create a downsampled data set of the terated Henon map.
data <- henon(2e6, params = c(-1.4, 0.3, 0.1 * runif(2)))
data <- data[sample.int(n = nrow(data), size = 2e5), ]

# Here we create the ATRIA nearest neighbor searcher object.
searcher = create_searcher(data, metric="manhattan", cluster_max_points = 64)

# We need to get an idea about typical small and large distances
# in the data set in order to generate distance bins.
dist.limits <- distlimits(searcher, data)
dist.breaks <- logspace(dist.limits[1], dist.limits[2], 32)

# Invoke the correlation sum computation.
res <- corrsum(searcher=searcher,
               data=data,
               dist.breaks = dist.breaks,
               min.actual.pairs = 2000,
               min.nr.samples.at.scale = 256,
               max.nr.samples.at.scale = 1024,
               batch.size = 128)

# Print and visualize results.
print(res)
x <- log2(res$dists)
y <-log2(res$correlation.sum)
plot(x, y)

# Fit a linear model to the data in log-log scale. The slope should give
# us an estimate of the correlation dimension.
print(lm(y ~ x, data = data.frame(x = x, y = y)[2:18,]))

# Delete the searcher object.
release_searcher(searcher)
