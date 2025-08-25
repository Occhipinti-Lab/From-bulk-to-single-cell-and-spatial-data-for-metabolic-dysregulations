# -*- coding: utf-8 -*-
"""Pre-process Gene Data
"""

#Import Packages used for Dataset
import pandas as pd #Converting Arrays
import numpy as np  #Arithemetic Operations
import matplotlib.pyplot as plt #visualization
from sklearn.model_selection import train_test_split

from google.colab import drive
drive.mount('/content/drive')

#Importing data
dataset = pd.read_csv('/content/drive/MyDrive/Colab Notebooks/TCGA-noiso-ENSG.csv')

nan_value = float("NaN")
dataset.replace(0, nan_value, inplace=True)
dataset.replace("", nan_value, inplace=True)
print(dataset.shape)

dataset.dropna(subset = ['EnsemblgeneID'], inplace=True)
print(dataset.shape)

perc = 70.0
min_count =  int(((100-perc)/100)*dataset.shape[1] + 1)
dataset = dataset.dropna( axis=0, 
                    thresh=min_count)

from sklearn.impute import KNNImputer

imputer = KNNImputer(n_neighbors=2)
dataset.iloc[:,1:] = imputer.fit_transform(dataset.iloc[:,1:])

dataset.iloc[:,1:]

d=dataset.T

#Importing data
sdata = pd.read_csv('/content/drive/MyDrive/Colab Notebooks/Clinical_Original.csv')

newdata=sdata[sdata.type == 'BRCA']
print(newdata.shape)

newdata=newdata[newdata.gender == 'FEMALE']
print(newdata.head())
print(newdata.tail())
print(newdata.shape)

d.reset_index(level=0, inplace=True)
new_header = d.iloc[0] 
d = d[1:] 
d.columns = new_header
d['EnsemblgeneID'] = d['EnsemblgeneID'].str[:-3]
print(d.head())
print(newdata.head())

newdata.reset_index(level=0, inplace=True)
print(newdata.head())
data=pd.DataFrame()
t=pd.DataFrame()
print(len(newdata))
for i in range(len(newdata)):
  data= d[newdata.bcr_patient_barcode[i]==d.EnsemblgeneID]
  t= t.append(data)
print(t)
print(newdata.bcr_patient_barcode[504])

t=t.drop_duplicates(subset='EnsemblgeneID', keep='first', inplace=False)

print(t.shape)

Data1=newdata[['bcr_patient_barcode','OS.time']]
print(Data1.head())
t=t.rename(columns ={'EnsemblgeneID':'bcr_patient_barcode'})
mergedata=pd.merge(t, Data1, on = 'bcr_patient_barcode')

print(mergedata.head())

Finaldata=mergedata.interpolate()
print(Finaldata.head())
Finaldata.isnull().sum()
print(Finaldata.shape)

m = np.mean(Finaldata['OS.time'].values)
print(m)
Finaldata.loc[(Finaldata['OS.time'] <= m), 'OS.time'] =0
Finaldata.loc[(Finaldata['OS.time'] > m), 'OS.time'] =1
print(Finaldata['OS.time'])

# Rescale data (between 0 and 1)
import pandas as pd 
from numpy import set_printoptions
from sklearn.preprocessing import MinMaxScaler
# separate array into input and output components
X = Finaldata.iloc[:,1:-1]
Y = Finaldata.iloc[:,-1]
scaler = MinMaxScaler(feature_range=(0, 1))
rescaledX = scaler.fit_transform(X)
# summarize transformed data
set_printoptions(precision=3)
print(rescaledX[0:5,:])

norm_data = Finaldata
norm_data.reset_index(level=0, inplace=True)
norm_data.head()

norm_data.drop([ 'index'], axis=1, inplace=True)

df = pd.DataFrame(rescaledX)

nrow,ncolumn = df.shape
for i in range(ncolumn):
  norm_data.iloc[:,i+1] = df.iloc[:,i]
print(norm_data.head())

norm_data=norm_data.filter(regex='^(?!NaN).+', axis=1)


"""Part 2: Feature selection - mRMR"""
# Instal mrmr package
# pip install git+https://github.com/smazzanti/mrmr

from mrmr import mrmr_classif


X=norm_data.iloc[:,1:-1]
y=norm_data.iloc[:,-1]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

# use mrmr classification
selected_features = mrmr_classif(X, y, K = 1000)

columns=selected_features +['OS.time']
print(columns)
data2=norm_data[columns]
data2.tail(20)

data2.to_csv (r'/content/drive/MyDrive/Colab Notebooks/GE1500.csv', index= False, header= True)
