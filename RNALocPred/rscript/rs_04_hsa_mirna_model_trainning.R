library(mlr)
library(tidyverse)
library(ggthemes)
library(parallelMap)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

hsa_mirna_feature <- read_rds(path = file.path(human_path_mirna, "rds_04_hsa_mirna_feature.rds.gz"))

#####################################
# two calsses
# in cell or out cell use all dataset
# ##################################
hsa_mirna_feature_cell <- 
  hsa_mirna_feature %>% 
  mutate(subloc = ifelse(subloc %in% c("Exosome", "Microvesicle", "Extracellular vesicle"), "Extracell", "Innercell")) %>% 
  as.data.frame()


hsa_mirna_cell_task <- 
  makeClassifTask(
    id = "hsa_mirna", 
    data = hsa_mirna_feature_cell, 
    target = "subloc"
    )
# for test.
write_rds(hsa_mirna_cell_task, path = file.path(human_path_mirna, "te_rds_01_hsa_mirna_task.rds.gz"), compress = "gz")

# feature selection
cell_task_feature_value <- generateFilterValuesData(task = hsa_mirna_cell_task, method = "information.gain")
# selected feature will lose some information.
plotFilterValues(cell_task_feature_value) +
  labs(
    title = "Human miRNA 29 Features, filter = information.gain"
  ) + 
  theme_gdocs() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9))
ggsave(filename = "fig_08_hsa_seed_feature_selection.png", path = human_path_mirna, device = "png")

# ----------------------basic benchmark experiment with seach hyperparameter--------------------------
# learners pool
lrns <- makeLearners(cls = c("rpart", "lda", "ksvm", "svm", "randomForest", "naiveBayes"), type = "classif", predict.type = "prob")
# resampling stratigies
cv10 <- makeResampleDesc(method = "CV", iters = 10, stratify = TRUE)
# measures
meas <- list(auc, mmce)
# benchmark
parallelStartSocket(cpus = 20, show.info = FALSE)
hsa_mirna_bmr <- benchmark(learners = lrns, tasks = hsa_mirna_cell_task, resamplings = cv10, measures = meas, show.info = FALSE)
parallelStop()

# algorithm comparison mmce
plotBMRBoxplots(bmr = hsa_mirna_bmr, measure = mmce) + theme_gdocs() + labs(x = "Algorithm") + aes(color = learner.id) + scale_color_gdocs()
ggsave(filename = "fig_09_hsa_seed_learner_comparison_1_mmce.png", path = human_path_mirna, device = "png")

# threshold vs performance
bmr_df <- generateThreshVsPerfData(hsa_mirna_bmr, measures = list(fpr, tpr, mmce, auc))
plt <- plotROCCurves(obj = bmr_df) 
plt + geom_text(aes(0.8,0.2, label = paste("AUC = ", round(auc,2))), check_overlap = T) + theme_gdocs()
ggsave(filename = "fig_09_hsa_seed_learner_comparison_2_auc.png", path = human_path_mirna, device = "png")

plotThreshVsPerf(bmr_df)


# based on the fig, the threshold should be set to 0.5

# -------------------tuning hyperparameter to promote performance---------------------
# Tuning paramters 
# Hyperparameter tuning is very slow, skip the hypertuning.
#-----------------------------------------------------------------------
# ksvm_ps = makeParamSet(
#   makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
#   makeDiscreteParam("kernel", values = c("vanilladot", "polydot", "rbfdot")),
#   makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x, requires = quote(kernel == "rbfdot"))
# )
# ksvm_rdesc <-makeResampleDesc("Holdout", stratify = TRUE)
# ksvm_ctrl <- makeTuneControlRandom(maxit = 100)
# ksvm_ctrl = makeTuneControlIrace(maxExperiments = 200L)
# ksvm_lrn <- makeTuneWrapper(learner = makeLearner("classif.ksvm", predict.type = "prob"), resampling = cv10, measures = meas, par.set = ksvm_ps, control = ksvm_ctrl, show.info = FALSE)
# parallelStartSocket(cpus = 10)
# res <- tuneParams(learner = ksvm_lrn, task = hsa_mirna_cell_task, resampling = ksvm_rdesc, measures = meas, par.set = ksvm_ps, control = ksvm_ctrl)
# parallelStop()
# lrns <- list(ksvm_lrn, makeLearner("classif.rpart", predict.type = "prob"))
# parallelStartMulticore(cpus = 50)
# bmr <- benchmark(learners = lrns, tasks = hsa_mirna_cell_task, resamplings = cv10, measures = meas)
# parallelStop()




# parameter tuning