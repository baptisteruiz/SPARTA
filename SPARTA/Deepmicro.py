# importing numpy, pandas
from unittest import TextTestRunner
import numpy as np
import random

import pandas as pd

# importing sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score

import shap

# importing util libraries
import datetime
import time
import os
import joblib


class DeepMicrobiome(object):
    def __init__(self, data, seed, data_dir, rand_state):
        self.t_start = time.time()
        self.filename = str(data)
        self.data = str(data)
        self.seed = seed
        self.rand_state = rand_state # To seed the gridSearch
        self.data_dir = data_dir
        self.prefix = ''
        self.representation_only = False


    def loadCustomDataWithLabels(self, Y_data, label_data, Ytest_ext, dtype=None, seed_valid_split = 0):
        # read file

        raw = Y_data
        #raw['random'] = np.random.random(size = len(raw))
        label = label_data

        # label data validity check
        if not len(label.shape) > 1:
            label_flatten = label.reshape((label.shape[0]))
        else:
            print('FileSpecificationError: The label file contains more than 1 column.')
            exit()
        
        indices = np.arange(len(label_flatten))

        # print("seed :" + str(seed_valid_split) )
        # train and test split
        self.X_train, self.X_test, self.y_train, self.y_test, self.indices_train, self.indices_test = train_test_split(raw.values.astype(dtype),
                                                                                label_flatten.astype('int'), indices, test_size=len(Ytest_ext.flatten()),
                                                                                random_state=seed_valid_split,
                                                                                stratify=label_flatten)

    # Classification
    def classification(self, Xtest_ext, Ytest_ext, seed, perf_dict, hyper_parameters, best_feature_records, var_ranking_method, method='svm', cv=5, scoring='roc_auc', n_jobs=1, cache_size=10000, best_auc=0, threshold_opt=0, multi_param = False):
        clf_start_time = time.time()

        # print("# Tuning hyper-parameters")
        # print(self.X_train.shape, self.y_train.shape)

        ###Test to adjust for multiclass classification
        if len(np.unique(np.concatenate((self.y_train, self.y_test)))) > 2:
             scoring = "roc_auc_ovr"
             multi_param = TextTestRunner

        # Support Vector Machine
        if method == 'svm':
            clf = GridSearchCV(SVC(probability=True, cache_size=cache_size, random_state=self.rand_state[0]), hyper_parameters, cv=StratifiedKFold(cv, shuffle=True, random_state=self.rand_state[1]), scoring=scoring, n_jobs=n_jobs, verbose=0)
            clf.fit(self.X_train, self.y_train)

        # Random Forest
        if method == 'rf':
            clf = GridSearchCV(RandomForestClassifier(n_jobs=-1, random_state=self.rand_state[0], class_weight='balanced'), hyper_parameters, cv=StratifiedKFold(cv, shuffle=True, random_state=self.rand_state[1]), scoring=scoring, n_jobs=n_jobs, verbose=0)
            clf.fit(self.X_train, self.y_train)

            if var_ranking_method == 'gini':
                best_features = pd.DataFrame(clf.best_estimator_.feature_importances_)
                best_feature_records.append(best_features)
                # with open(self.data_dir + "results/" + self.data + "_best_features_random_rf.txt", 'a') as f:
                #     best_features.to_csv(f, header=None)

            elif var_ranking_method == 'shap':
                explainer = shap.TreeExplainer(clf.best_estimator_)
                shap_values = explainer.shap_values(Xtest_ext)
                shap_df = pd.DataFrame(shap_values[:,:,1])
                shaps_summed = pd.DataFrame(shap_df.sum(axis=0))
                best_feature_records.append(shaps_summed)

        # Evaluate performance of the best model on test set
        y_prob = clf.predict_proba(self.X_test)

        # Get best threshold for predictions (Youden J method)
        fpr, tpr, thresholds = roc_curve(self.y_test, y_prob[:, 1], pos_label=1)
        youdenJ = tpr - fpr
        index = np.argmax(youdenJ)
        thresholdOpt = round(thresholds[index], ndigits = 4)

        #Class predictions on each set
        y_control = (clf.predict_proba(self.X_train)[:,1] >= thresholdOpt).astype('int')
        y_true, y_pred = self.y_test, (clf.predict_proba(self.X_test)[:,1] >= thresholdOpt).astype('int')
        y_pred_ext = (clf.predict_proba(Xtest_ext)[:,1] >= thresholdOpt).astype('int')

        #print("y_prob: ", y_prob)
        # Performance Metrics: AUC, ACC, Recall, Precision, F1_score
        if multi_param:
            metrics = [round(roc_auc_score(y_true, y_prob, multi_class = "ovr"), 4),
                    round(accuracy_score(y_true, y_pred), 4),
                    round(recall_score(y_true, y_pred, average = "micro"), 4),
                    round(precision_score(y_true, y_pred, average = "micro"), 4),
                    round(f1_score(y_true, y_pred, average = "micro"), 4), ]

        else:
            metrics = [round(roc_auc_score(y_true, y_prob[:, 1]), 4),
                   round(accuracy_score(y_true, y_pred), 4),
                   round(recall_score(y_true, y_pred), 4),
                   round(precision_score(y_true, y_pred), 4),
                   round(f1_score(y_true, y_pred), 4), ]

        #print("metric : ", metrics[0], " vs ", best_auc)

        ###ONLY UNCOMMENT IF YOU WANT TO CHECK THE PERFORMANCES OF ALL TRAINED MODELS.
        #joblib.dump(clf, self.data_dir + self.data + "_saved_classifier_run_"+str(seed)+".joblib")

        if metrics[0] >= best_auc:
            best_auc = metrics[0]
            threshold_opt = thresholdOpt
            saved_classifier_filepath = os.path.join(self.data_dir, self.data + "_saved_classifier.joblib")
            joblib.dump(clf, saved_classifier_filepath)
            #print("Model saved!")

        # time stamp
        metrics.append(str(datetime.datetime.now()))

        # running time
        metrics.append(round( (time.time() - self.t_start), 2))

        # classification time
        metrics.append(round( (time.time() - clf_start_time), 2))

        # best hyper-parameter append
        metrics.append(str(clf.best_params_))

        ####PUT BEST FEATURES BACK HERE IF ERROR

        if multi_param:
            perf_dict[str(seed)] = [str(clf.best_params_), self.indices_test, thresholdOpt, round(roc_auc_score(self.y_train, clf.predict_proba(self.X_train), multi_class = "ovr"), 4), round(roc_auc_score(y_true, y_prob, multi_class = "ovr"), 4), round(roc_auc_score(Ytest_ext, clf.predict_proba(Xtest_ext), multi_class = "ovr"), 4)]
        else:
            perf_dict[str(seed)] = [str(clf.best_params_), self.indices_test, thresholdOpt, round(roc_auc_score(self.y_train, clf.predict_proba(self.X_train)[:,1]), 4), round(roc_auc_score(y_true, y_prob[:, 1]), 4), round(roc_auc_score(Ytest_ext, clf.predict_proba(Xtest_ext)[:,1]), 4)]


        return (best_auc, threshold_opt, perf_dict, best_feature_records)


# run exp function
def run_exp(seed, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, Y_data, label_data, dataset_name,
            best_feature_records, data_dir, method, var_ranking_method, real_seed, seed_valid):

    random.seed(real_seed)
    vec_rand = random.sample(range(1000), 2)
    dm = DeepMicrobiome(data=dataset_name, seed=seed, data_dir=data_dir, rand_state = vec_rand)
    random.seed(seed_valid)
    seed_valid_split = random.sample(range(1000), 1)[0]


    dm.loadCustomDataWithLabels(Y_data=Y_data, label_data=label_data, Ytest_ext=Ytest_ext, dtype=None, seed_valid_split = seed_valid_split )

    # time check after data has been loaded
    dm.t_start = time.time()

    rf_hyper_parameters = [{'n_estimators': [s for s in range(100, 1001, 200)],
                            'max_features': ['sqrt', 'log2'],
                            'min_samples_leaf': [1, 2, 3, 4, 5],
                            'criterion': ['gini']
                            }, ]
    
    svm_hyper_parameters = [{'C': [2 ** s for s in range(-5, 6, 2)], 'kernel': ['linear']},
                            {'C': [2 ** s for s in range(-5, 6, 2)], 'gamma': [2 ** s for s in range(3, -15, -2)],'kernel': ['rbf']}]
    
    if method == 'rf':
        best_auc, threshold_opt, perf_dict, best_feature_records = dm.classification(Xtest_ext, Ytest_ext, seed, perf_dict, rf_hyper_parameters, best_feature_records, var_ranking_method, method='rf', cv=5,
                            n_jobs=-2, scoring='roc_auc', best_auc=best_auc, threshold_opt=threshold_opt)
        
    elif method == 'svm':
        best_auc, threshold_opt, perf_dict, best_feature_records = dm.classification(Xtest_ext, Ytest_ext, seed, perf_dict, svm_hyper_parameters, best_feature_records, var_ranking_method, method='svm', cv=5,
                            n_jobs=-2, scoring='roc_auc', best_auc=best_auc, threshold_opt=threshold_opt)

    # Create dict containing training and validation sets.
    labels_sets = {}
    labels_sets['training_set'] = ','.join(Y_data.iloc[dm.indices_train].index.tolist())
    labels_sets['validation_set'] = ','.join(Y_data.iloc[dm.indices_test].index.tolist())

    labels_sets['training_set_labels'] = ','.join(map(str, label_data[dm.indices_train]))
    labels_sets['validation_set_labels'] = ','.join(map(str, label_data[dm.indices_test]))

    return (best_auc, threshold_opt, perf_dict, best_feature_records, labels_sets)
