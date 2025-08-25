# -*- coding: utf-8 -*-
"""Multi-objective feature selection algorithm
by Le Doan
"""
# Install neccessary package
# !pip install lifelines
# !pip install scikit-optimize
# !pip install mygene

# Check scikit-learn version
import sklearn
print('The scikit-learn version is {}.'.format(sklearn.__version__))

""" Step 1: Preprare the problem """
# 1.1. Load the libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from skopt import BayesSearchCV
from sklearn.metrics import  classification_report, confusion_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score

from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import OrdinalEncoder
from lifelines.statistics import logrank_test, multivariate_logrank_test, pairwise_logrank_test

from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts

# import mygene

sns.set_style("white")
import warnings
warnings.filterwarnings("ignore")

# Load clinical, fluxomic and gene data
file1 = pd.read_csv('Data/clinical_data.csv',  index_col=0)
file2 = pd.read_csv('Data/FS_Q3.csv')
file3 = pd.read_csv('Data/Processed_gene.csv', index_col=0)
# Load kegg results
kegg_gene =  pd.read_csv(r'Data/kegg_7000.txt', sep='\t')


"""<h3> Preprocess clinical"""

file1 = file1.rename(columns={ 'vital_status':'OS.Status'})
file3 = file3.drop( ['OS.time'], axis = 1)
file2.index = file2['bcr_patient_barcode']
file2.drop( ['bcr_patient_barcode'], axis = 1, inplace=True)

file1.loc[(file1['initial_pathologic_dx_year'] >= 1988) & (file1['initial_pathologic_dx_year'] <= 2000) , 'initial_pathologic_dx_year'] = 1
file1.loc[(file1['initial_pathologic_dx_year'] >= 2001) & (file1['initial_pathologic_dx_year'] <= 2005) , 'initial_pathologic_dx_year'] = 2

year = file1[['initial_pathologic_dx_year']]
enc = OrdinalEncoder()
year_enc = enc.fit_transform(year)

file1['initial_pathologic_dx_year'] = year_enc

"""<h3> Preprocess gene"""


# Select KEGG pathway
# Select significant pathways
gene_list = kegg_gene.query('Benjamini  <0.05')
# .nsmallest(5, columns='Benjamini' )

gene_list_kegg =[]
for g in gene_list['Genes']:
    g_list = g.split(", ")
    gene_list_kegg.append(g_list)

select_gene = [number for sublist in gene_list_kegg for number in sublist]
gene_kegg = set(select_gene)
print(f'Total KEGG selected gene: {len(gene_kegg)}')

# Select mrmr gene
mrmr_gene = set(file3.columns)

inter_gene = list(gene_kegg & mrmr_gene)
print(f'Total combined selected gene: {len(list(gene_kegg & mrmr_gene))}')
df_gene = file3[inter_gene]

"""<h3> Preprocess fluxomic data"""

# Convert reaction name

match =pd.read_csv('Data/Match_ID.csv')
reaction_name = pd.DataFrame(file2.columns, columns=['ID'])
reaction_name = pd.merge(reaction_name,match, how="inner", on=['ID'])
recon_name = list(reaction_name['MAR ID'])
file2.columns =recon_name

data = pd.merge(pd.merge(file1,file2, how="inner", left_index=True, right_index=True),df_gene, how="inner", left_index=True, right_index=True)

# Check missing data
print(data.isna().sum().sum())

df = data.query("`OS.time`!=0 | `OS.Status`!=0")
df.shape

"""<h3> ML pipeline"""

"""ML model function"""

import time
from sklearn.feature_selection import SelectKBest,  f_classif, chi2
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer

# Parameter
RANDOM_STATE =12
cv=5


# Modality features
clinical = ['age_at_initial_pathologic_diagnosis',
 'initial_pathologic_dx_year',
 'race',
 'ajcc_pathologic_tumor_stage',
 'histological_type',
 'menopause_status',
 'OS.Status',
 'tumor_status',
 'margin_status']
gene = df_gene.columns.to_list()
flux = file2.columns.to_list()


# ### pipeline each modality
clinical_transformer = Pipeline(steps=[
    ('features', SelectKBest(chi2, k=9)),
    ('scaler', StandardScaler())])

gene_transformer = Pipeline(steps=[
    ('features', SelectKBest(f_classif)),
    ('scaler', StandardScaler())])

flux_transformer = Pipeline(steps=[
    ('features', SelectKBest(f_classif)),
    ('scaler', StandardScaler())])

preprocessor = ColumnTransformer(
    transformers=[
        ("clinical", clinical_transformer, clinical),
        ("gene", gene_transformer, gene),
        ("flux", gene_transformer, flux),
    ], verbose_feature_names_out=False, # added this line
)

# Use BayesSearchCV to find the best parameters:
# Multi-objective feature selection algorithm
def model_best_estimator(x_train, y_train, RANDOM_STATE=12, cv=5):

    # SVC
    t0 = time.time()

    SVC_grid = {
    'preprocessor__clinical__features__k': list(range(3,10,2)),
    'preprocessor__gene__features__k': list(range(20,500,20)),
    'preprocessor__flux__features__k': list(range(20,200,20)),
    'classifier__C': [ 0.01, 0.1, 1, 10, 100], #10, 100
    'classifier__kernel': ["rbf", "linear" ], #"sigmoid", 'poly'
    'classifier__gamma': [0.001, 0.01, 0.1, 1, 100]
    }

    svc_pipe = Pipeline(
    steps=[("preprocessor", preprocessor), ("classifier", SVC(random_state=RANDOM_STATE, max_iter=10000))])

    grid_log_reg = BayesSearchCV(svc_pipe,
                                SVC_grid, cv=cv, n_jobs=-1)
    grid_log_reg.fit(x_train, y_train)

    # get the SVC with the best parameters.
    svc_model = grid_log_reg.best_estimator_
    t1 = time.time()

    print("Best fit parameter for SVC", svc_model)
    print("Elapsed time {:.2f} s".format(t1 - t0))



    # DecisionTree Classifier:
    t4 = time.time()

    tree_params_grid = {
    'preprocessor__clinical__features__k': list(range(3,10,2)),
    'preprocessor__gene__features__k': list(range(20,500,20)),
    'preprocessor__flux__features__k': list(range(20,200,20)),
    'classifier__criterion': [ "gini", "entropy"],
    'classifier__min_samples_leaf': [4 , 8],
    'classifier__max_depth':list(range(3,30,3))
    }

    dt_pipe = Pipeline(
    steps=[("preprocessor", preprocessor), ("classifier", DecisionTreeClassifier(random_state=RANDOM_STATE))])


    grid_tree = BayesSearchCV(dt_pipe,
                             tree_params_grid, cv=cv)
    grid_tree.fit(x_train, y_train)

    # tree best estimator
    tree_clf = grid_tree.best_estimator_
    t5 = time.time()

    print("\nBest fit parameter for Decision Tree:", tree_clf)
    print("Elapsed time {:.2f} s".format(t5 - t4))

    # Random Forest Classifier
    t6 = time.time()

    rf_params_grid = {
    'preprocessor__clinical__features__k': list(range(3,10,2)),
    'preprocessor__gene__features__k': list(range(20,500,20)),
    'preprocessor__flux__features__k': list(range(20,200,20)),
    'classifier__criterion': [ "gini", "entropy"],
    'classifier__min_samples_leaf': [4, 8 ],
    'classifier__max_depth':list(range(3,30,3)),
    'classifier__n_estimators': [200]
    }

    rf_pipe = Pipeline(
    steps=[("preprocessor", preprocessor), ("classifier", RandomForestClassifier(random_state=RANDOM_STATE))])

    grid_rf = BayesSearchCV(rf_pipe,
                           rf_params_grid, cv=cv)
    grid_rf.fit(x_train, y_train)

    # random forest best estimator
    rf = grid_rf.best_estimator_
    t7 = time.time()

    print("\nBest fit parameter for Random Forest:", rf)
    print("Elapsed time {:.2f} s".format(t7 - t6))


    # GBoost Classifier
    t8 = time.time()
    gb_para = {
    'preprocessor__clinical__features__k': list(range(3,10,2)),
    'preprocessor__gene__features__k': list(range(20,500,20)),
    'preprocessor__flux__features__k': list(range(20,200,20)),
    'classifier__learning_rate':[0.01, 0.1, 1],
    'classifier__min_samples_leaf': [4, 8 ],  # 8 looks good
    'classifier__n_estimators': [200]
    }

    gb_pipe = Pipeline(
    steps=[("preprocessor", preprocessor), ("classifier", GradientBoostingClassifier(random_state=RANDOM_STATE))])


    grid_gb = BayesSearchCV(gb_pipe,
                           gb_para, cv=cv)
    grid_gb.fit(x_train, y_train)

    # GB best estimator
    gb = grid_gb.best_estimator_
    t9 = time.time()

    print("\nBest fit parameter for Random Forest:", gb)
    print("Elapsed time {:.2f} s".format(t9 - t8))


    return [svc_model,  tree_clf, rf, gb]   #knn, nb, tree_clf,


# Evaluate model by using cross validation
def evaluate_model(classifier, x_train, y_train, cv=5):
    classifier.fit(x_train, y_train)
    score = cross_val_score(classifier, x_train, y_train, cv=cv)
    return score

# Get training model results
def train_model(classifier, x_train, y_train, cv=5):
    y_train_pre = cross_val_predict(classifier, x_train, y_train, cv=cv)
    print(classification_report(y_train, y_train_pre, labels=[1,0]))


# Get testing model results
def predict_model(classifier, x_test, y_test):
    y_pre = classifier.predict(x_test)
    print(classification_report(y_test, y_pre, labels=[1,0]))

    # Confusion Matrix
    print('Confusion matrix:', classifier)
    cf_matrix = confusion_matrix(y_test, y_pre, labels=[1,0])
    ax =sns.heatmap(cf_matrix, annot=True, fmt="d", cmap="Blues",
                xticklabels=['High Risk', 'Low Risk'],
                yticklabels=['High Risk', 'Low Risk'])
    ax.set(xlabel="Predicted outputs", ylabel = "Actual outputs")
    plt.show()

    return y_pre


# 4.4 Set up function for training and testing flow
# Start with find the best parameter for ML model, train and get result + visualisation results

def train_test(X_train, y_train, X_test, y_test, RANDOM_STATE=RANDOM_STATE, cv=5):
    # Find best parameter for model
    model_select_result = model_best_estimator(X_train, y_train)

    # log_reg, knn, nb, tree_clf, rf, gb = model_select_result
    log_reg, tree_clf, rf, gb = model_select_result

    results =[]
    accuracy_list =[]
    f1_list =[]


    # Train and get result
    for classifier in model_select_result:
        # print("\nPredict model:", classifier)
        # evaluate_model(classifier, X_train, y_train, cv=cv)
        # print("\nTraining result:")
        # train_model(classifier, X_train, y_train, cv=cv)
        print("Testing result:")
        y_pre = predict_model(classifier, X_test, y_test)
        results.append(y_pre)
        accuracy_list.append(accuracy_score(y_test, y_pre))
        f1_list.append(f1_score(y_test, y_pre, average='weighted'))


    #can add feature importance
    return [log_reg, tree_clf,  rf, gb, results, accuracy_list, f1_list] #knn, nb, tree_clf,

""" Step 4: Evaluate Algorithms"""

drop_list = ['OS.time', 'Risk']

n=30
np.random.seed(1)
seeds = np.random.permutation(1500)[:n]

model_list =[]
accuracy_list = []
f1_list = []

for s in seeds:
    # print(s)
    # Split train test
    X_train1, X_test1 = train_test_split(df, test_size=0.2,
                                     stratify=df['OS.Status'],
                                    random_state=s,
                                    shuffle=True)

    # Get median of X_train and assign Risk variable
    med_train= X_train1['OS.time'].median()
    print('Median value of OS.time in training set', med_train)
    X_train1['Risk'] = np.where(X_train1['OS.time'] >= med_train, 0, 1)
    X_test1['Risk'] = np.where(X_test1['OS.time'] >= med_train, 0, 1)


    X_train = X_train1.drop(drop_list, axis = 1)
    X_test = X_test1.drop(drop_list, axis = 1)

    y_train = X_train1['Risk'].values
    y_test = X_test1['Risk'].values

    print(X_train.shape)
    print(X_test.shape)



    print(f"\n************* Evaluate models seed {s}**************")

    svc_model, tree_clf,  rf, gb, results, acc, f1_res = train_test(X_train, y_train, X_test, y_test)
    # results = train_test(X_train, y_train, X_test, y_test)
    accuracy_list.append(acc)
    f1_list.append(f1_res)
    model_list.append([svc_model, tree_clf, rf, gb])

model_name =['SVM',  'DT', 'RF', 'GB']
pd.DataFrame(accuracy_list, columns=model_name, index=seeds).to_csv('Results/acc_multi.csv')
pd.DataFrame(f1_list, columns=model_name, index=seeds).to_csv('Results/f1_multi.csv')

pd.DataFrame(accuracy_list, columns=model_name, index=seeds)
# .to_csv('Results/acc_multi.csv')
pd.DataFrame(f1_list, columns=model_name, index=seeds)
# .to_csv('Results/f1_multi.csv')

import pickle
for i in range(len(seeds)):
    svc_model,  tree_clf, rf, gb = model_list[i]
    pickle.dump(svc_model, open(f"Model/Multi_SVM_{seeds[i]}.pkl", "wb"))
    pickle.dump(tree_clf, open(f"Model/Multi_DT_{seeds[i]}.pkl", "wb"))
    pickle.dump(rf, open(f"Model/Multi_RF_{seeds[i]}.pkl", "wb"))
    pickle.dump(gb, open(f"Model/Multi_Gb_{seeds[i]}.pkl", "wb"))

