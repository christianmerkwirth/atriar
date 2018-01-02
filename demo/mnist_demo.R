library(h5)

this.file <- c('/Users/cmerk/Downloads/mnist-784-euclidean.hdf5')

file <- h5file(this.file, mode = 'r')
train <- file['train'][]
test <- file['test'][]
expected.distances <- file['distances'][]
expected.neighbors <- file['neighbors'][]

data = rbind(train, test)

searcher = create_searcher(data, metric="manhattan", cluster_max_points = 64)

dist.breaks <- logspace(2000, 24000, 32)

res <- corrsum(searcher=searcher,
               data=data,
               dist.breaks = dist.breaks,
               min.actual.pairs = 2000,
               min.nr.samples.at.scale = 128,
               max.nr.samples.at.scale = 512,
               batch.size = 64)

# Print and visualize results.
print(res)
x <- log2(res$dists)
y <-  log2(res$correlation.sum)
plot(x, y)

# Fit a linear model to the data in log-log scale. The slope should give
# us an estimate of the correlation dimension.
print(lm(y ~ x, data = data.frame(x = x, y = y)[2:18,]))

# Delete the searcher object.
release_searcher(searcher)



