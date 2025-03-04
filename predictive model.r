library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3verse)

# feature selection
# fetch the samples carrying top20/bottom20 TSG in t cell therapy screen
mlr_data <- as.data.frame(mouse_data)
mlr_data$response <- factor(t.the[rownames(mouse_data),'sample'])

task <- as_task_classif(mlr_data,target = 'response',positive = '1')
learner_rf <- lrn("classif.ranger",importance='impurity',
                  predict_type = "response")
afs = auto_fselector(
  fselector = fs("random_search"),
  learner = learner_rf,
  resampling = rsmp ("cv",folds=5),
  measure = msr('classif.acc'),
  term_evals = 4)
afs$train(task)
features=afs$fselect_result$features
# obtain a set of 28 biomarkers
# subset their expression and combine with TME composition in human data
human_data=cbind(human_data[,features[[1]]],human_tme_cibersort)
mlr_test <- as.data.frame(human_data)
mlr_test$response <- factor(ifelse(meta[rownames(mlr_test),]$response=='1','1','-1'))
mlr_task1 <- as_task_classif(mlr_test,target = 'response',positive = '1')
# split train data 
split <- partition(mlr_task1, ratio = 0.8)  
train_task <- mlr_task1$clone()
test_task <- mlr_task1$clone()
train_task$filter(split$train)
test_task$filter(split$test)

# set the ensemble model
learner_gausspr <- lrn('classif.gausspr',predict_type='prob',
                       kernel=to_tune(c('rbfdot', 'polydot', 'vanilladot', 'tanhdot', 'laplacedot', 'besseldot', 'anovadot', 'splinedot')),
                       sigma=to_tune(0.1,10),
                       tol=to_tune(0.000001,0.1))
learner_bart <- lrn('classif.bart',predict_type='prob',
                    ntree=to_tune(100,800),
                    k=to_tune(0,10))
learner_rf <- lrn("classif.ranger",
                  num.trees = to_tune(100,1000),
                  predict_type = "prob",
                  mtry=to_tune(2,4),
                  min.node.size=to_tune(8,20),
                  max.depth=to_tune(3,6))
learner_svm <- lrn("classif.svm",predict_type = "prob",
                   type = "C-classification",
                   cost=to_tune(0.1,10),
                   gamma=to_tune(0.1,10),
                   kernel=to_tune(c("polynomial","radial","sigmoid")),
                   degree=to_tune(1,3))
learner_kknn <- lrn('classif.kknn',predict_type='prob',
                    k= to_tune(3,50),
                    distance=to_tune(0,3),
                    kernel=to_tune(c("rectangular", "gaussian", "rank", "optimal")))
learner_bayes <- lrn('classif.naive_bayes',predict_type = "prob",
                     laplace=to_tune(0,1))
learner_xgboost <- lrn("classif.xgboost",predict_type = "prob",
                       eta=to_tune(0,1),
                       gamma=to_tune(0,5),
                       max_depth=to_tune(1,8),
                       min_child_weight=to_tune(5,20),
                       subsample=to_tune(0,1),
                       colsample_bytree=to_tune(0.5,1),
                       nrounds=to_tune(20,30),
                       eval_metric=to_tune(c("auc","error",'logloss')))
learner_catboost <- lrn('classif.catboost',predict_type = "prob",
                        depth=to_tune(1,16),
                        bagging_temperature=to_tune(0,5),
                        learning_rate=to_tune(0.001,1),
                        l2_leaf_reg=to_tune(0,10))
learner_ksvm <- lrn('classif.ksvm',predict_type = "prob",
                    type=to_tune(c('C-svc', 'C-bsvc', 'spoc-svc', 'kbb-svc')),
                    kernel=to_tune(c('rbfdot', 'polydot', 'vanilladot', 'laplacedot', 'besseldot', 'anovadot')),
                    C=to_tune(-0.5,10))
learners = list(learner_rf,learner_gausspr,learner_bart,
                learner_bayes, learner_svm,learner_ksvm,
                learner_kknn,learner_xgboost,learner_catboost)

granp <- ppl('branch',learners)
#granp$plot()
glearner<-as_learner(granp)
#glearner$param_set

at<-auto_tuner(tuner = tnr("random_search"),
               learner = glearner,
               resampling = rsmp("cv",folds=10),
               measure = msr("classif.auc"),
               term_evals = 50)
at$train(train_task)

# fetch the auc in train dataset
best_learner <- at$learner
resampling = rsmp("cv",folds=5)
rr = resample(task = train_task,learner = best_learner,
              resampling = resampling , store_models = TRUE)
rr$score(msr('classif.auc'))$classif.auc
# fetch the auc in test dataset
tmpp=at$predict(test_task)
roc_curve <- roc(factor(tmpp$truth), tmpp$prob[,1])
roc_curve$auc

