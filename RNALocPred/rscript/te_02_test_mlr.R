library(mlr)
library(tidyverse)
library(ggthemes)
root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

hsa_mirna_task <- read_rds(path = file.path(human_path_mirna, "te_rds_01_hsa_mirna_task.rds.gz"))

# hsa_mirna_task
# 
# calculate the performance measures
n = getTaskSize(bh.task)
lrn = makeLearner("regr.gbm", n.tress = 1000)
mod = train(lrn, bh.task, subset = seq(1, n, 2))
pred = predict(mod, task = bh.task, subset = seq(2, n, 2))
pred
performance(pred = pred, measures = list(medse, mse, mae, timetrain), model = mod)


# cluster
lrn = makeLearner("cluster.kmeans", centers = 3)
mod = train(lrn, mtcars.task)
pred = predict(mod, task = mtcars.task)
performance(pred, measures = dunn, task = mtcars.task)

# ORC AUC needs posterior probabilities are predicted.
# # measures only suitable for binary problems
lrn = makeLearner("classif.rpart", predict.type = "prob")
mod = train(lrn, task = sonar.task)
pred = predict(mod, task = sonar.task)
performance(pred, measures = auc)

lrn = makeLearner("classif.lda", predict.type = "prob")
n = getTaskSize(sonar.task)
mod = train(lrn, task = sonar.task, subset = seq(1,n,2))
pred = predict(mod, task = sonar.task, subset = seq(2, n, 2))

pred
performance(pred, measures = list(fpr, fnr, mmce))
d = generateThreshVsPerfData(pred, measures = list(fpr, fnr, mmce))
plotThreshVsPerf(d)
plotThreshVsPerfGGVIS(d)

r = calculateROCMeasures(pred)

################
#Resampling
################
r = resample(learner = "regr.lm",task = bh.task, resampling = cv10)

rdesc = makeResampleDesc(method = "Subsample", iters = 5)
lrn = makeLearner("classif.rpart", parms = list(split = "information"))
r = resample(
  learner = lrn, 
  task = sonar.task, 
  resampling = rdesc, 
  measures = list(mmce, fpr, fnr, timetrain),
  models = TRUE,
  show.info = FALSE)

r = resample(
  learner = "classif.lda", 
  task = iris.task, 
  resampling = makeResampleDesc("CV", iters = 3, stratify = TRUE),
  show.info = FALSE)


# resampleInstance
rdesc = makeResampleDesc("CV", iters = 3)
rin = makeResampleInstance(rdesc, task = iris.task)

r.lda = resample(learner = "classif.lda", task = iris.task, resampling = rdesc, show.info = FALSE)
r.rpart = resample(learner = "classif.rpart", task = iris.task, resampling = rdesc, show.info = FALSE)



############################
# Tuning Hyperparameters
#########################
# search space for C
ps <- makeParamSet(makeNumericParam("C", lower = 0.01, upper = 0.1))

# random search with 100 iterations
ctrl <- makeTuneControlRandom(maxit = 100L)

rdesc <- makeResampleDesc("CV", iters = 3L)

# Specifying the seach space
discrete_ps <- makeParamSet(
  makeDiscreteParam("C", values = c(0.5, 1.0, 1.5, 2.0)),
  makeDiscreteParam("sigma", values = c(0.5, 1.0, 1.5, 2.0))
)

num_ps <- makeParamSet(
  makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10^x),
  makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10^x)
)


ctrl <- makeTuneControlGrid()
ctrl <- makeTuneControlGrid(resolution = 15L)
rdesc <- makeResampleDesc("CV", iters = 3L)



discrete_ps <- makeParamSet(
  makeDiscreteParam("C", values = c(0.5, 1.0, 1.5, 2.0)),
  makeDiscreteParam("sigma", values = c(0.5, 1.0, 1.5, 2.0))
)



ctrl <- makeTuneControlGrid(resolution = 1000L)
rdesc <- makeResampleDesc("CV", iters = 5L)

res <- tuneParams("classif.ksvm", task = iris.task, resampling = rdesc, par.set = discrete_ps, control = ctrl)


num_ps <- makeParamSet(
  makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10^x),
  makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10^x)
)
ctrl = makeTuneControlRandom(maxit = 100L)
res = tuneParams("classif.ksvm", task = iris.task, resampling = rdesc, par.set = num_ps, control = ctrl, measures = list(acc, setAggregation(acc, test.sd)))
res = tuneParams("classif.ksvm", task = iris.task, resampling = rdesc, par.set = num_ps,
                 control = ctrl, measures = list(acc, mmce), show.info = FALSE)
data = generateHyperParsEffectData(res)
plotHyperParsEffect(data, x = "iteration", y = "acc.test.mean",
                    plot.type = "line")
lrn <- setHyperPars(makeLearner("classif.ksvm"), par.vals = res$x)
lrn
m <- train(lrn, task = iris.task)
predict(m, task = iris.task)

#################
#benchmark
#################

lrns <- list(makeLearner("classif.lda"), makeLearner("classif.rpart"))

rdesc <- makeResampleDesc("Holdout")

bmr <- benchmark(learners = lrns, tasks = sonar.task, resamplings = rdesc)


lrns2 <- list(makeLearner("classif.randomForest"), makeLearner("classif.qda"))
bmr2 <- benchmark(lrns2, sonar.task, rdesc)

mergeBenchmarkResults(list(bmr, bmr2))
rin = getBMRPredictions(bmr)[[1]][[1]]$instance
rin

# big example
lrns <- list(
  makeLearner("classif.lda", id = "lda"),
  makeLearner("classif.rpart", id = "rpart"),
  makeLearner("classif.randomForest", id = "randomForest")
)

#get task from mlbench
ring.task <- convertMLBenchObjToTask("mlbench.ringnorm", n = 600)
wave.task <- convertMLBenchObjToTask("mlbench.waveform", n = 600)

tasks <- list(iris.task, sonar.task, pid.task, ring.task, wave.task)

rdesc <- makeResampleDesc("CV", iters = 10)
meas <- list(mmce, ber, timetrain)

bmr <- benchmark(learners = lrns, tasks = tasks, resamplings = rdesc, measures = meas)

perf <- getBMRPerformances(bmr, as.df = T)

plt <- plotBMRBoxplots(bmr, measure = mmce)
plotBMRBoxplots(bmr, measure = ber,  pretty.names = FALSE) 

levels(plt$data$task.id) = c("Iris", "Ringnorm", "Waveform", "Diabetes", "Sonar")
levels(plt$data$learner.id) = c("LDA", "CART", "RF")

plt + ylab("Error rate")


plotBMRSummary(bmr)


m <- convertBMRToRankMatrix(bmr, mmce)
plotBMRRanksAsBarChart(bmr, pos = "tile")
plotBMRSummary(bmr, trafo = "rank", jitter = 0)

plotBMRRanksAsBarChart(bmr, mmce, pos = "dodge")
friedmanTestBMR(bmr)
friedmanPostHocTestBMR(bmr, p.value = 0.1)

g = generateCritDifferencesData(bmr, p.value = 0.1, test = "nemenyi")
plotCritDifferences(g) + coord_cartesian(xlim = c(-1,5), ylim = c(0,2))


g = generateCritDifferencesData(bmr, p.value = 0.1, test = "bd", baseline = "randomForest")
plotCritDifferences(g) + coord_cartesian(xlim = c(-1,5), ylim = c(0,2))

#custom plots
perf <- getBMRPerformances(bmr, as.df = TRUE)
qplot(mmce, colour = learner.id, facets = . ~ task.id,
      data = perf[perf$task.id %in% c("iris-example", "Sonar-example"),], geom = "density") +
  theme(strip.text.x = element_text(size = 8))

df = reshape2::melt(perf, id.vars = c("task.id", "learner.id", "iter"))
df = df[df$variable != "ber",]
head(df)

qplot(variable, value, data = df, colour = learner.id, geom = "boxplot",
      xlab = "measure", ylab = "performance") +
  facet_wrap(~ task.id, nrow = 2)

perf = getBMRPerformances(bmr, task.id = "Sonar-example", as.df = TRUE)
df = reshape2::melt(perf, id.vars = c("task.id", "learner.id", "iter"))
df = df[df$variable == "mmce",]
df = reshape2::dcast(df, task.id + iter ~ variable + learner.id)
head(df)
GGally::ggpairs(df, 3:5)

start <- proc.time()
rdesc = makeResampleDesc("CV", iters = 50)
r = resample("classif.lda", iris.task, rdesc)
proc.time() - start

parallelStartMulticore(10)
start <- proc.time()
rdesc = makeResampleDesc("CV", iters = 50)
r = resample("classif.lda", iris.task, rdesc)
proc.time() - start
parallelStop()


# Advanced tuning
ps <- makeParamSet(
  makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeDiscreteParam("kernel", values = c("vanilladot", "polydot", "rbfdot")),
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x, requires = quote(kernel == "rbfdot")),
  makeIntegerParam("degree", lower = 2L, upper = 5L, requires = quote(kernel == "polydot"))
)

ctrl <- makeTuneControlIrace(maxExperiments = 200L)
rdesc <- makeResampleDesc("Holdout")

parallelStartSocket(10, show.info = FALSE)
start <- proc.time()
res <- tuneParams("classif.ksvm", iris.task, rdesc, par.set = ps, control = ctrl, show.info = FALSE)
proc.time() - start
parallelStop()

names(res)
res$x
res$y
res$learner
res$control
res$threshold
head(as.data.frame(res$opt.path))

# Tuning across whole model spaces
base.learners <-
  list(
    makeLearner("classif.ksvm"),
    makeLearner("classif.randomForest")
  )
lrn <- makeModelMultiplexer(base.learners = base.learners)
ps <- makeModelMultiplexerParamSet(
  lrn,
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeIntegerParam("ntree", lower = 1L, upper = 500L)
)
print(ps)
rdesc <- makeResampleDesc("CV", iters = 2L)
ctrl <- makeTuneControlIrace(maxExperiments = 200L)

parallelStartSocket(20, show.info = FALSE)
start <- proc.time()
res <- tuneParams(lrn, iris.task, resampling = rdesc, par.set = ps, control = ctrl, show.info = FALSE)
proc.time() - start
parallelStop()
print(head(as.data.frame(res$opt.path)))

# multi criteria
ps = makeParamSet(
  makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x)
)
ctrl = makeTuneMultiCritControlRandom(maxit = 30L)
rdesc = makeResampleDesc("Holdout")
res = tuneParamsMultiCrit("classif.ksvm", task = sonar.task, resampling = rdesc, par.set = ps, measures = list(fpr, fnr), control = ctrl, show.info = FALSE)
res

plotTuneMultiCritResult(res)

# Feature selection
# calculating the feature importance
fv <- generateFilterValuesData(iris.task, method = "information.gain")
fv2 <- generateFilterValuesData(iris.task, method = c("information.gain", "chi.squared"))

plotFilterValues(fv2)

filtered.task <- filterFeatures(iris.task, method = "information.gain", abs = 2)
filtered.task <- filterFeatures(iris.task, fval = fv, perc = 0.25)
filtered.task <- filterFeatures(iris.task, fval = fv, threshold = 0.5)


lrn = makeFilterWrapper(learner = "classif.fnn", fw.method = "information.gain", fw.abs = 2)
rdesc = makeResampleDesc("CV", iters = 10)
r = resample(learner = lrn, task = iris.task, resampling = rdesc, show.info = FALSE, models = TRUE)
r$aggr

sfeats <- sapply(r$models, getFilteredFeatures)
sfeats


################
################
################

