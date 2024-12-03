######### training a classifier to annotate tumor/normal cells
######### Chi-Yun Wu

# Import required libraries
import pandas as pd
import numpy as np 
from numpy.random import seed
import matplotlib.pyplot as plt
import sklearn
import scanpy as sc
import os
#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri

# Import necessary modules
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from keras.models import model_from_json
from math import sqrt


# Keras specific
import keras
from keras.models import Sequential
from keras.layers import Dense
#from keras.utils import to_categorical 
from tensorflow.keras.utils import to_categorical
from keras.callbacks import EarlyStopping, ModelCheckpoint
import tensorflow as tf

import pickle as pkl
import json

from numpy.random import seed
seed(123)
import tensorflow



df=pd.read_csv("~/matrix_all_normalN3_integrated_train.txt", sep='\t') ##manually annotated cells for training
tensorflow.random.set_seed(123)

filepath='~/integrated_label3N_best_model_1layer_node1024_relu_p10_b10000_all_shffled.h5'
modelpath='~/model_integrated_label3N_1024nodes_shuffled_all'
historypath='~/history_label3N_1024nodes_shuffled.json'

#df=df.rename(columns = {'Unnamed: 3379':'label'})
print(df.shape)
label=df['label'].tolist()
label2=[]
for ii in label:
	if ii == 'T':
		label2.append(1)
	elif ii == 'N':
		label2.append(0)
	else:
		label2.append(NA)

df['label']=label2


df2=df.sample(frac=1) ## for shuffling
#df.describe()


target_column = ['label'] 
predictors = list(set(list(df.columns))-set(target_column))
predictors.sort()

#df[predictors] = df[predictors]/df[predictors].max()
#df.describe()

###### anndata
adata=df2[predictors]
adata = sc.AnnData(adata)
#adata.obs
adata.raw = adata.copy()
n_counts = adata.X.sum(axis = 1)
#sc.pp.normalize_per_cell(adata, counts_per_cell_after = 10000)
#adata.obs['size_factors'] = n_counts / 10000
adata.obs['label']=df2[target_column]

#sc.pp.log1p(adata)
#sc.pp.scale(adata)

rname=adata.obs.index
datasets=[]
for ii in rname:
	datasets.append(ii.split("___")[0])

dataset_uni=[]
for ii in datasets:
	if ii not in dataset_uni:
		dataset_uni.append(ii)

datasets=np.array(datasets)


#accu_all=[]
dname='all'
ind=[]
nind=[i for i in range(df.shape[0])]


train=adata[nind,:]
#test=adata[ind,:]


X_train=train.X#df[predictors].iloc[nind,].values
y_train=train.obs['label'].values#df[target_column].iloc[nind,].values

print(X_train.shape); #print(X_test.shape)


y_train = to_categorical(y_train)



####

count_classes = y_train.shape[1]
print(count_classes)

input_dim=X_train.shape[1]


model = Sequential()
model.add(Dense(1024, activation='relu', input_dim=input_dim))
model.add(Dense(2, activation='softmax'))

# Compile the model
model.compile(optimizer='adam', 
              loss='categorical_crossentropy', 
              metrics=['accuracy'])


callbacks = [EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True), ModelCheckpoint(filepath=filepath, monitor='val_loss', save_best_only=True)]

history = model.fit(X_train, y_train, validation_split=0.1, epochs=200, batch_size=10000, verbose=1, callbacks=callbacks)


# save model
model.save(modelpath)#, save_format='tf')

# save history
history_dict = history.history
# Save it under the form of a json file
#historypath="/home/mnt/nzh/nzhanglab/project/cywu/CPTCA/NB/classifier/integrated/history_test.json"
json.dump(history_dict, open(historypath, 'w'))

#history_dict2 = json.load(open(historypath, 'r'))


## test
model = keras.models.load_model('~/model_integrated_label3N_1024nodes_shuffled_all/')
test=adata[range(74526, 82807),:]
X_test=test.X#df[predictors].iloc[ind,].values
y_test=test.obs['label'].values#df[target_column].iloc[ind,].values
y_test = to_categorical(y_test)

pred_test = model.predict(X_test)
scores2 = model.evaluate(X_test, y_test, verbose=0)
dname
print('Accuracy on test data: {}% \n Error on test data: {}'.format(scores2[1], 1 - scores2[1]))
X_test.shape    
scores2

historypath='~/'+'history_label3N_1024nodes_shuffled.json'
history_dict=json.load(open(historypath, 'r'))
history_dict['val_loss'][len(history_dict['val_loss'])-11]




