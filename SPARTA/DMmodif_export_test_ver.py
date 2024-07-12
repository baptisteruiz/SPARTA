# importing numpy, pandas, and matplotlib
from unittest import TextTestRunner
import numpy as np
from numpy.core.arrayprint import _none_or_positive_arg
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
# import tensorflow as tf



import sys
import tqdm

# importing sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.random_projection import GaussianRandomProjection
from sklearn import cluster
from sklearn import preprocessing
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score

import shap


# from sklearn.manifold import TSNE
# import umap

# importing keras
# import keras
# import keras.backend as K
# from keras.wrappers.scikit_learn import KerasClassifier
# from keras.callbacks import EarlyStopping, ModelCheckpoint, LambdaCallback
# from keras.models import Model, load_model

# importing util libraries
import datetime
import time
import math
import os
import importlib
import joblib

# importing custom library
#import DNN_models
import exception_handle

#fix np.random.seed for reproducibility in numpy processing
#np.random.seed(7)


class DeepMicrobiome(object):
    def __init__(self, data, seed, data_dir):
        self.t_start = time.time()
        self.filename = str(data)
        self.data = str(data)
        self.seed = seed
        self.data_dir = data_dir
        self.prefix = ''
        self.representation_only = False

        

    def loadData(self, feature_string, label_string, label_dict, dtype=None):
        # read file
        filename = self.data_dir + "data/" + self.filename
        if os.path.isfile(filename):
            raw = pd.read_csv(filename, sep='\t', index_col=0, header=None)
        else:
            print("FileNotFoundError: File {} does not exist".format(filename))
            exit()

        # select rows having feature index identifier string
        X = raw.loc[raw.index.str.contains(feature_string, regex=False)].T
        #X['random'] = np.random.random(size = len(X))
        
        # get class labels
        Y = raw.loc[label_string] #'disease'
        Y = Y.replace(label_dict)

        # train and test split
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X.values.astype(dtype), Y.values.astype('int'), test_size=0.2, random_state=self.seed, stratify=Y.values)
        self.printDataShapes()

    def loadCustomData(self, dtype=None):
        # read file
        filename = self.data_dir + "data/" + self.filename
        if os.path.isfile(filename):
            raw = pd.read_csv(filename, sep=',', index_col=False, header=None)
            #raw['random'] = np.random.random(size = len(raw))
        else:
            print("FileNotFoundError: File {} does not exist".format(filename))
            exit()


        # load data
        self.X_train = raw.values.astype(dtype)

        # put nothing or zeros for y_train, y_test, and X_test
        self.y_train = np.zeros(shape=(self.X_train.shape[0])).astype(dtype)
        self.X_test = np.zeros(shape=(1,self.X_train.shape[1])).astype(dtype)
        self.y_test = np.zeros(shape=(1,)).astype(dtype)
        self.printDataShapes(train_only=True)

    def loadCustomDataWithLabels(self, Y_data, label_data, Ytest_ext, dtype=None):
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

        print('DeepMicro input Y data:', raw)
        print('DeepMicro input labels:', label_flatten)

        indices = np.arange(len(label_flatten))
        # train and test split
        self.X_train, self.X_test, self.y_train, self.y_test, self.indices_train, self.indices_test = train_test_split(raw.values.astype(dtype),
                                                                                label_flatten.astype('int'), indices, test_size=len(Ytest_ext.flatten()),
                                                                                random_state=self.seed,
                                                                                stratify=label_flatten)
        
        
        print('Validation indices:', self.indices_test)
        #self.printDataShapes()




    # #ADDED BY US: t-SNE projection
    # def tsne(self):

    #     self.prefix = self.prefix + 't-SNE_'

    #     X_train = TSNE().fit_transform(self.X_train)
    #     X_test = TSNE().fit_transform(self.X_test)

    #     # applying the eigenvectors to the whole training and the test set.
    #     self.X_train = X_train
    #     self.X_test = X_test
    #     self.printDataShapes()

    # #ADDED BY US: UMAP projection
    # def umap(self):
    #     self.prefix = self.prefix + 'UMAP_'
    #     reducer = umap.UMAP()
        
    #     X_train = reducer.fit_transform(self.X_train)
    #     X_test = reducer.fit_transform(self.X_test)

    #     # applying the eigenvectors to the whole training and the test set.
    #     self.X_train = X_train
    #     self.X_test = X_test
    #     self.printDataShapes()

    # #Principal Component Analysis
    # def pca(self, ratio=0.99):
    #     # manipulating an experiment identifier in the output file
    #     self.prefix = self.prefix + 'PCA_'

    #     # PCA
    #     pca = PCA()
    #     pca.fit(self.X_train)
    #     n_comp = 0
    #     ratio_sum = 0.0

    #     for comp in pca.explained_variance_ratio_:
    #         ratio_sum += comp
    #         n_comp += 1
    #         if ratio_sum >= ratio:  # Selecting components explaining 99% of variance
    #             break

    #     pca = PCA(n_components=n_comp)
    #     pca.fit(self.X_train)

    #     X_train = pca.transform(self.X_train)
    #     X_test = pca.transform(self.X_test)

    #     # applying the eigenvectors to the whole training and the test set.
    #     self.X_train = X_train
    #     self.X_test = X_test
    #     self.printDataShapes()

    # #Gausian Random Projection
    # def rp(self):
    #     # manipulating an experiment identifier in the output file
    #     self.prefix = self.prefix + 'RandP_'
    #     # GRP
    #     rf = GaussianRandomProjection(eps=0.5)
    #     rf.fit(self.X_train)

    #     # applying GRP to the whole training and the test set.
    #     self.X_train = rf.transform(self.X_train)
    #     self.X_test = rf.transform(self.X_test)
    #     self.printDataShapes()

    # #Shallow Autoencoder & Deep Autoencoder
    # def ae(self, dims = [50], epochs= 40000, batch_size=100, verbose=2, loss='mean_squared_error', latent_act=False, output_act=False, act='relu', patience=20, val_rate=0.2, no_trn=False, opt="adam"):

    #     # manipulating an experiment identifier in the output file
    #     if patience != 20:
    #         self.prefix += 'p' + str(patience) + '_'
    #     if len(dims) == 1:
    #         self.prefix += 'AE'
    #     else:
    #         self.prefix += 'DAE'
    #     if loss == 'binary_crossentropy':
    #         self.prefix += 'b'
    #     if latent_act:
    #         self.prefix += 't'
    #     if output_act:
    #         self.prefix += 'T'
    #     self.prefix += str(dims).replace(", ", "-") + '_'
    #     if act == 'sigmoid':
    #         self.prefix = self.prefix + 's'

    #     # filename for temporary model checkpoint
    #     modelName = self.prefix + self.data + '.h5'

    #     # clean up model checkpoint before use
    #     if os.path.isfile(modelName):
    #         os.remove(modelName)

    #     # callbacks for each epoch
    #     callbacks = [EarlyStopping(monitor='val_loss', patience=patience, mode='min', verbose=1),
    #                  ModelCheckpoint(modelName, monitor='val_loss', mode='min', verbose=1, save_best_only=True)]

    #     # spliting the training set into the inner-train and the inner-test set (validation set)
    #     X_inner_train, X_inner_test, y_inner_train, y_inner_test = train_test_split(self.X_train, self.y_train, test_size=val_rate, random_state=self.seed, stratify=self.y_train)

    #     # insert input shape into dimension list
    #     dims.insert(0, X_inner_train.shape[1])

    #     # create autoencoder model
    #     self.autoencoder, self.encoder = DNN_models.autoencoder(dims, act=act, latent_act=latent_act, output_act=output_act)
    #     self.autoencoder.summary()

    #     if no_trn:
    #         return

    #     # compile model
    #     self.autoencoder.compile(optimizer=opt, loss=loss)

    #     # fit model
    #     self.history = self.autoencoder.fit(X_inner_train, X_inner_train, epochs=epochs, batch_size=batch_size, callbacks=callbacks,
    #                          verbose=verbose, validation_data=(X_inner_test, X_inner_test))
    #     # save loss progress
    #     self.saveLossProgress()

    #     # load best model
    #     self.autoencoder = load_model(modelName)
    #     layer_idx = int((len(self.autoencoder.layers) - 1) / 2)
    #     self.encoder = Model(self.autoencoder.layers[0].input, self.autoencoder.layers[layer_idx].output)

    #     # applying the learned encoder into the whole training and the test set.
    #     self.X_train = self.encoder.predict(self.X_train)
    #     self.X_test = self.encoder.predict(self.X_test)

    # # Variational Autoencoder
    # def vae(self, dims = [10], epochs=40000, batch_size=100, verbose=2, loss='mse', output_act=False, act='relu', patience=25, beta=1.0, warmup=True, warmup_rate=0.01, val_rate=0.2, no_trn=False):

    #     # manipulating an experiment identifier in the output file
    #     if patience != 25:
    #         self.prefix += 'p' + str(patience) + '_'
    #     if warmup:
    #         self.prefix += 'w' + str(warmup_rate) + '_'
    #     self.prefix += 'VAE'
    #     if loss == 'binary_crossentropy':
    #         self.prefix += 'b'
    #     if output_act:
    #         self.prefix += 'T'
    #     if beta != 1:
    #         self.prefix += 'B' + str(beta)
    #     self.prefix += str(dims).replace(", ", "-") + '_'
    #     if act == 'sigmoid':
    #         self.prefix += 'sig_'

    #     # filename for temporary model checkpoint
    #     modelName = self.prefix + self.data + '.h5'

    #     # clean up model checkpoint before use
    #     if os.path.isfile(modelName):
    #         os.remove(modelName)

    #     # callbacks for each epoch
    #     callbacks = [EarlyStopping(monitor='val_loss', patience=patience, mode='min', verbose=1),
    #                  ModelCheckpoint(modelName, monitor='val_loss', mode='min', verbose=1, save_best_only=True,save_weights_only=True)]

    #     # warm-up callback
    #     warm_up_cb = LambdaCallback(on_epoch_end=lambda epoch, logs: [warm_up(epoch)])  # , print(epoch), print(K.get_value(beta))])

    #     # warm-up implementation
    #     def warm_up(epoch):
    #         val = epoch * warmup_rate
    #         if val <= 1.0:
    #             K.set_value(beta, val)
    #     # add warm-up callback if requested
    #     if warmup:
    #         beta = K.variable(value=0.0)
    #         callbacks.append(warm_up_cb)

    #     # spliting the training set into the inner-train and the inner-test set (validation set)
    #     X_inner_train, X_inner_test, y_inner_train, y_inner_test = train_test_split(self.X_train, self.y_train,
    #                                                                                 test_size=val_rate,
    #                                                                                 random_state=self.seed,
    #                                                                                 stratify=self.y_train)

    #     # insert input shape into dimension list
    #     dims.insert(0, X_inner_train.shape[1])

    #     # create vae model
    #     self.vae, self.encoder, self.decoder = DNN_models.variational_AE(dims, act=act, recon_loss=loss, output_act=output_act, beta=beta)
    #     self.vae.summary()

    #     if no_trn:
    #         return

    #     # fit
    #     self.history = self.vae.fit(X_inner_train, epochs=epochs, batch_size=batch_size, callbacks=callbacks, verbose=verbose, validation_data=(X_inner_test, None))

    #     # save loss progress
    #     self.saveLossProgress()

    #     # load best model
    #     self.vae.load_weights(modelName)
    #     self.encoder = self.vae.layers[1]

    #     # applying the learned encoder into the whole training and the test set.
    #     _, _, self.X_train = self.encoder.predict(self.X_train)
    #     _, _, self.X_test = self.encoder.predict(self.X_test)



    # # Convolutional Autoencoder
    # def cae(self, dims = [32], epochs=40000, batch_size=100, verbose=2, loss='mse', output_act=False, act='relu', patience=25, val_rate=0.2, rf_rate = 0.1, st_rate = 0.25, no_trn=False, opt="adam"):

    #     # manipulating an experiment identifier in the output file
    #     self.prefix += 'CAE'
    #     if loss == 'binary_crossentropy':
    #         self.prefix += 'b'
    #     if output_act:
    #         self.prefix += 'T'
    #     self.prefix += str(dims).replace(", ", "-") + '_'
    #     if act == 'sigmoid':
    #         self.prefix += 'sig_'

    #     # filename for temporary model checkpoint
    #     modelName = self.prefix + self.data + '.h5'

    #     # clean up model checkpoint before use
    #     if os.path.isfile(modelName):
    #         os.remove(modelName)

    #     # callbacks for each epoch
    #     callbacks = [EarlyStopping(monitor='val_loss', patience=patience, mode='min', verbose=1),
    #                  ModelCheckpoint(modelName, monitor='val_loss', mode='min', verbose=1, save_best_only=True,save_weights_only=True)]


    #     # fill out blank
    #     onesideDim = int(math.sqrt(self.X_train.shape[1])) + 1
    #     enlargedDim = onesideDim ** 2
    #     self.X_train = np.column_stack((self.X_train, np.zeros((self.X_train.shape[0], enlargedDim - self.X_train.shape[1]))))
    #     self.X_test = np.column_stack((self.X_test, np.zeros((self.X_test.shape[0], enlargedDim - self.X_test.shape[1]))))

    #     # reshape
    #     self.X_train = np.reshape(self.X_train, (len(self.X_train), onesideDim, onesideDim, 1))
    #     self.X_test = np.reshape(self.X_test, (len(self.X_test), onesideDim, onesideDim, 1))
    #     self.printDataShapes()

    #     # spliting the training set into the inner-train and the inner-test set (validation set)
    #     X_inner_train, X_inner_test, y_inner_train, y_inner_test = train_test_split(self.X_train, self.y_train,
    #                                                                                 test_size=val_rate,
    #                                                                                 random_state=self.seed,
    #                                                                                 stratify=self.y_train)

    #     # insert input shape into dimension list
    #     dims.insert(0, (onesideDim, onesideDim, 1))

    #     # create cae model
    #     self.cae, self.encoder = DNN_models.conv_autoencoder(dims, act=act, output_act=output_act, rf_rate = rf_rate, st_rate = st_rate)
    #     self.cae.summary()
    #     if no_trn:
    #         return

    #     # compile
    #     self.cae.compile(optimizer=opt, loss=loss)

    #     # fit
    #     self.history = self.cae.fit(X_inner_train, X_inner_train, epochs=epochs, batch_size=batch_size, callbacks=callbacks, verbose=verbose, validation_data=(X_inner_test, X_inner_test, None))

    #     # save loss progress
    #     self.saveLossProgress()

    #     # load best model
    #     self.cae.load_weights(modelName)
    #     if len(self.cae.layers) % 2 == 0:
    #         layer_idx = int((len(self.cae.layers) - 2) / 2)
    #     else:
    #         layer_idx = int((len(self.cae.layers) - 1) / 2)
    #     self.encoder = Model(self.cae.layers[0].input, self.cae.layers[layer_idx].output)

    #     # applying the learned encoder into the whole training and the test set.
    #     self.X_train = self.encoder.predict(self.X_train)
    #     self.X_test = self.encoder.predict(self.X_test)
    #     self.printDataShapes()

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
            clf = GridSearchCV(SVC(probability=True, cache_size=cache_size, random_state=0), hyper_parameters, cv=StratifiedKFold(cv, shuffle=True, random_state=0), scoring=scoring, n_jobs=n_jobs, verbose=0)
            clf.fit(self.X_train, self.y_train)



        # Random Forest
        if method == 'rf':
            clf = GridSearchCV(RandomForestClassifier(n_jobs=-1, random_state=0, class_weight='balanced'), hyper_parameters, cv=StratifiedKFold(cv, shuffle=True, random_state=0), scoring=scoring, n_jobs=n_jobs, verbose=0)
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


        # # Multi-layer Perceptron
        # if method == 'mlp':
        #     model = KerasClassifier(build_fn=DNN_models.mlp_model, input_dim=self.X_train.shape[1], verbose=0, )
        #     clf = GridSearchCV(estimator=model, param_grid=hyper_parameters, cv=StratifiedKFold(cv, shuffle=True), scoring=scoring, n_jobs=n_jobs, verbose=100, multiclass = 'ovr')
        #     clf.fit(self.X_train, self.y_train, batch_size=32)

        # print("Best parameters set found on development set:")
        # print()
        # print(clf.best_params_)

        

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

        # print("pred:",y_pred)
        # print("real:",y_true)
        
        ##export details of predictions vs real into text file

        # with open(self.data_dir + "results/" + self.data + "_compare_pred_vs_test.txt", 'a') as f:
        #     f.writelines('TRAINING: \n')
        #     f.writelines('pred:' + str(y_control) +'\n')
        #     f.writelines('real:' + str(self.y_train) +'\n')
        #     f.writelines('VALIDATION: \n')
        #     f.writelines('pred:' + str(y_pred) +'\n')
        #     f.writelines('real:' + str(y_true) +'\n')
        #     f.writelines('Test set indices: ' + str(self.indices_test) +'\n')
        #     f.writelines('EXTERNAL TEST SET: \n')
        #     f.writelines('pred:' + str(y_pred_ext) +'\n')
        #     f.writelines('real:' + str(Ytest_ext.reshape((Ytest_ext.shape[0]))) +'\n')
        

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
            joblib.dump(clf, self.data_dir + '/' + self.data + "_saved_classifier.joblib")
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
            

        # Write performance metrics as a file (old output)
        # res = pd.DataFrame([metrics], index=[self.prefix + method])

        
        # with open(self.data_dir + "results/" + self.data + "_result.txt", 'a') as f:
        #     res.to_csv(f, header=None)
        
        # with open(self.data_dir + "results/" + self.data + "_test_indices.txt", 'a') as f:
        #     self.indices_test

        # print('Accuracy metrics')
        # print('AUC, ACC, Recall, Precision, F1_score, time-end, runtime(sec), classfication time(sec), best hyper-parameter')
        # print(metrics)

        if multi_param:
            perf_dict[str(seed)] = [str(clf.best_params_), self.indices_test, thresholdOpt, round(roc_auc_score(self.y_train, clf.predict_proba(self.X_train), multi_class = "ovr"), 4), round(roc_auc_score(y_true, y_prob, multi_class = "ovr"), 4), round(roc_auc_score(Ytest_ext, clf.predict_proba(Xtest_ext), multi_class = "ovr"), 4)]
        else:
            perf_dict[str(seed)] = [str(clf.best_params_), self.indices_test, thresholdOpt, round(roc_auc_score(self.y_train, clf.predict_proba(self.X_train)[:,1]), 4), round(roc_auc_score(y_true, y_prob[:, 1]), 4), round(roc_auc_score(Ytest_ext, clf.predict_proba(Xtest_ext)[:,1]), 4)]

        return(best_auc, threshold_opt, perf_dict, best_feature_records)

    def printDataShapes(self, train_only=False):
        print("X_train.shape: ", self.X_train.shape)
        if not train_only:
            print("y_train.shape: ", self.y_train.shape)
            print("X_test.shape: ", self.X_test.shape)
            print("y_test.shape: ", self.y_test.shape)

    # ploting loss progress over epochs
    def saveLossProgress(self):
        #print(self.history.history.keys())
        #print(type(self.history.history['loss']))
        #print(min(self.history.history['loss']))

        loss_collector, loss_max_atTheEnd = self.saveLossProgress_ylim()

        # save loss progress - train and val loss only
        figureName = self.prefix + self.data + '_' + str(self.seed)
        plt.ylim(min(loss_collector)*0.9, loss_max_atTheEnd * 2.0)
        plt.plot(self.history.history['loss'])
        plt.plot(self.history.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train loss', 'val loss'],
                   loc='upper right')
        plt.savefig(self.data_dir + "results/" + figureName + '.png')
        plt.close()

        if 'recon_loss' in self.history.history:
            figureName = self.prefix + self.data + '_' + str(self.seed) + '_detailed'
            plt.ylim(min(loss_collector) * 0.9, loss_max_atTheEnd * 2.0)
            plt.plot(self.history.history['loss'])
            plt.plot(self.history.history['val_loss'])
            plt.plot(self.history.history['recon_loss'])
            plt.plot(self.history.history['val_recon_loss'])
            plt.plot(self.history.history['kl_loss'])
            plt.plot(self.history.history['val_kl_loss'])
            plt.title('model loss')
            plt.ylabel('loss')
            plt.xlabel('epoch')
            plt.legend(['train loss', 'val loss', 'recon_loss', 'val recon_loss', 'kl_loss', 'val kl_loss'], loc='upper right')
            plt.savefig(self.data_dir + "results/" + figureName + '.png')
            plt.close()

    # supporting loss plot
    def saveLossProgress_ylim(self):
        loss_collector = []
        loss_max_atTheEnd = 0.0
        for hist in self.history.history:
            current = self.history.history[hist]
            loss_collector += current
            if current[-1] >= loss_max_atTheEnd:
                loss_max_atTheEnd = current[-1]
        return loss_collector, loss_max_atTheEnd




# run exp function
def run_exp(seed, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, Y_data, label_data, dataset_name, best_feature_records, data_dir, method, var_ranking_method):


    dm = DeepMicrobiome(data=dataset_name, seed=seed, data_dir=data_dir)
    dm.loadCustomDataWithLabels(Y_data=Y_data, label_data=label_data, Ytest_ext=Ytest_ext, dtype=None)

    # time check after data has been loaded
    dm.t_start = time.time()

    # CHANGED BY US: parameters for autoencoder representation learning models
    # opt=tf.keras.optimizers.Adam(learning_rate=1e-5)

    # # Representation learning (Dimensionality reduction)
    # if args.pca:
    #     dm.pca()
    # if args.ae:
    #     dm.ae(dims=[int(i) for i in args.dims.split(',')], act=args.act, epochs=args.max_epochs, loss=args.aeloss,
    #           latent_act=args.ae_lact, output_act=args.ae_oact, patience=args.patience, no_trn=args.no_trn, opt=opt)
    # if args.vae:
    #     dm.vae(dims=[int(i) for i in args.dims.split(',')], act=args.act, epochs=args.max_epochs, loss=args.aeloss, output_act=args.ae_oact,
    #            patience= 25 if args.patience==20 else args.patience, beta=args.vae_beta, warmup=args.vae_warmup, warmup_rate=args.vae_warmup_rate, no_trn=args.no_trn)
    # if args.cae:
    #     dm.cae(dims=[int(i) for i in args.dims.split(',')], act=args.act, epochs=args.max_epochs, loss=args.aeloss, output_act=args.ae_oact,
    #            patience=args.patience, rf_rate = args.rf_rate, st_rate = args.st_rate, no_trn=args.no_trn, opt=opt)
    # if args.rp:
    #     dm.rp()

    # #ADDED BY US
    # if args.tsne:
    #     dm.tsne()

    # if args.umap:
    #     dm.umap()

    # write the learned representation of the training set as a file
    # if args.save_rep:
    #     if numRLrequired == 1:
    #         rep_file = dm.data_dir + "results/" + dm.prefix + dm.data + "_rep.csv"
    #         pd.DataFrame(dm.X_train).to_csv(rep_file, header=None, index=None)
    #         print("The learned representation of the training set has been saved in '{}'".format(rep_file))
    #     else:
    #         print("Warning: Command option '--save_rep' is not applied as no representation learning or dimensionality reduction has been conducted.")

    # Classification
    # if args.no_clf or (args.data == None and args.custom_data_labels == None):
    #     print("Classification task has been skipped.")
    # else:
    #     # turn off GPU
        # os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
        #importlib.reload(keras)

        # best model performances


    # training classification models
    
    # elif args.method == "rf":

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
    # elif args.method == "mlp":
    #     best_auc, threshold_opt = dm.classification(hyper_parameters=mlp_hyper_parameters, method='mlp', cv=args.numFolds,
    #                       n_jobs=args.numJobs, scoring=args.scoring, best_auc=best_auc, threshold_opt=threshold_opt)
    # elif args.method == "svm_rf":
    #     best_auc, threshold_opt = dm.classification(hyper_parameters=svm_hyper_parameters, method='svm', cv=args.numFolds,
    #                       n_jobs=args.numJobs, scoring=args.scoring, cache_size=args.svm_cache, best_auc=best_auc, threshold_opt=threshold_opt)
    #     best_auc, threshold_opt = dm.classification(hyper_parameters=rf_hyper_parameters, method='rf', cv=args.numFolds,
    #                       n_jobs=args.numJobs, scoring=args.scoring, best_auc=best_auc, threshold_opt=threshold_opt)
    # else:
    #     best_auc, threshold_opt = dm.classification(hyper_parameters=svm_hyper_parameters, method='svm', cv=args.numFolds,
    #                       n_jobs=args.numJobs, scoring=args.scoring, cache_size=args.svm_cache, best_auc=best_auc, threshold_opt=threshold_opt)
    #     best_auc, threshold_opt = dm.classification(hyper_parameters=rf_hyper_parameters, method='rf', cv=args.numFolds,
    #                       n_jobs=args.numJobs, scoring=args.scoring, best_auc=best_auc, threshold_opt=threshold_opt)
    #     best_auc, threshold_opt = dm.classification(hyper_parameters=mlp_hyper_parameters, method='mlp', cv=args.numFolds,
    #                       n_jobs=args.numJobs, scoring=args.scoring, best_auc=best_auc, threshold_opt=threshold_opt)

    return(best_auc, threshold_opt, perf_dict, best_feature_records)
