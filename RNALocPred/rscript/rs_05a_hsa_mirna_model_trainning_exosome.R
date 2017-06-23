library(mlr)
library(tidyverse)
library(ggthemes)
library(parallelMap)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

hsa_mirna_feature <- read_rds(path = file.path(human_path_mirna, "rds_06_hsa_mirna_features_exosome.rds.gz"))


task <- makeClassifTask(id = "hsa_mirna_exosome", data = hsa_mirna_feature %>% as.data.frame(), target = "subloc")

generateFilterValuesData(task = task, method = "information.gain") -> task_feature
plotFilterValues(fvalues = task_feature) +
  labs(title = "Human miRNA Features, filter = information.gain") +
  theme_gdocs() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9))

lrns <- makeLearners(cls = c("rpart", "lda", "ksvm", "svm", "randomForest", "naiveBayes"), type = "classif", predict.type = "prob")

cv10 <- makeResampleDesc(method = "CV", iters = 10, stratify = TRUE)
meas <- list(auc, mmce)
parallelStartSocket(cpus = 20, show.info = FALSE)
bmr <- benchmark(learners = lrns, tasks = task, resamplings = cv10, measures = meas, show.info = FALSE)
parallelStop()
ggsave(filename = "fig_11a_hsa_seed_learner_comparison_1_exosome_mmce.png", path = human_path_mirna, device = "png")
plotBMRBoxplots(bmr = bmr, measure = mmce) +
  theme_gdocs() +
  labs(x = "Algorithm") +
  aes(color = learner.id) +
  scale_color_gdocs()

bmr_df <- generateThreshVsPerfData(bmr, measures = list(fpr, tpr, mmce, auc)) 
plt <- plotROCCurves(obj = bmr_df)
plt + geom_text(aes(0.8, 0.2, label = paste("AUC = ", round(auc, 2))), check_overlap = T) + theme_gdocs()
ggsave(filename = "fig_11a_hsa_seed_learner_comparison_2_exosome_auc.png", path = human_path_mirna, device = "png")
plotThreshVsPerf(bmr_df)

# Hyperparameter tuning ksvm and randomForest
ps <- makeParamSet(
  makeNumericParam(id = "C", lower = -10, upper = 10, trafo = function(x) 10 ^ x),
  makeNumericParam(id = "sigma", lower = -10, upper = 10, trafo = function(x) 10 ^ x)
)
ctrl <- makeTuneControlGrid(resolution = 200L)
rdesc <- makeResampleDesc(method = "CV", iters = 10L, stratify = T)
parallelStartSocket(cpus = 40)
res <- tuneParams(learner = "classif.ksvm", task = task, resampling = rdesc, par.set = ps, control = ctrl, measures = list(mmce, acc))
parallelStop()
# Ksvm
ps = makeParamSet(
  makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeDiscreteParam("kernel", values = c("vanilladot", "polydot", "rbfdot")),
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x,
                   requires = quote(kernel == "rbfdot")),
  makeIntegerParam("degree", lower = 2L, upper = 5L,
                   requires = quote(kernel == "polydot"))
)

ctrl = makeTuneControlIrace(maxExperiments = 200L)
rdesc = makeResampleDesc("Holdout")
parallelStartSocket(cpus = 20, show.info = FALSE)
res = tuneParams("classif.ksvm", task, rdesc, par.set = ps, control = ctrl, show.info = FALSE)
parallelStop()
print(head(as.data.frame(res$opt.path)))



