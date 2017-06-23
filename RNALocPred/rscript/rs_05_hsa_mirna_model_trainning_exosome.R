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
hsa_mirna_feature_exosome <- 
  hsa_mirna_feature %>% 
  mutate(subloc = ifelse(subloc %in% c("Exosome"), "Exosome", "ExtraExosome")) %>% 
  as.data.frame()

hsa_mirna_cell_task <- 
  makeClassifTask(
    id = "hsa_mirna", 
    data = hsa_mirna_feature_exosome, 
    target = "subloc"
  )

# feature selection
cell_task_feature_value <- generateFilterValuesData(task = hsa_mirna_cell_task, method = "information.gain")
# selected feature will lose some information.
plotFilterValues(cell_task_feature_value) +
  labs(
    title = "Human miRNA 29 Features, filter = information.gain"
  ) + 
  theme_gdocs() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9))
# ggsave(filename = "fig_08_hsa_seed_feature_selection.png", path = human_path_mirna, device = "png")

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
ggsave(filename = "fig_11_hsa_seed_learner_comparison_1_exosome_mmce.png", path = human_path_mirna, device = "png")

# threshold vs performance
bmr_df <- generateThreshVsPerfData(hsa_mirna_bmr, measures = list(fpr, tpr, mmce, auc))
plt <- plotROCCurves(obj = bmr_df) 
plt + geom_text(aes(0.8,0.2, label = paste("AUC = ", round(auc,2))), check_overlap = T) + theme_gdocs()
ggsave(filename = "fig_11_hsa_seed_learner_comparison_2_exosome_auc.png", path = human_path_mirna, device = "png")

plotThreshVsPerf(bmr_df)
