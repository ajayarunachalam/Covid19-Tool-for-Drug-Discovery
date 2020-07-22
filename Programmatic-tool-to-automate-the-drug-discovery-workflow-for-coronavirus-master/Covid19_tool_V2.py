# %%
import pandas as pd
import numpy as np
import urllib
import os 
import requests
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import seaborn as seabornInstance 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn import metrics
get_ipython().magic(u'matplotlib inline')
import time

# %%
link1 = "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22protacxn%22:%22notnull%22},{%22cid%22:%22notnull%22},{%22repacxn%22:%22P0C6U8%22}]},%22order%22:[%22activity,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PROTACXN_P0C6U8_bioactivity_protein%22}"

# %%
for i in range(0,2):
    try:
        os.remove('downloaded1.csv')
        #print("Deleted old File")
        break
    except Exception as e:
        #print("No file ")
        break
    else:
        break

# %%
for i in range(2):
    try:
        data = pd.read_csv(link1)
        break
    #except IncompleteRead as I:
     #   print("Server Overloading , Proceeding")
      #  break
    except Exception as a:
        print(str(a)+" is the error , Trying {} time".format(i))
        continue
    else:
        break
else:
    print("something Wrong , Try running Again [refer error code for more]")

# %%
#data = pd.read_csv("downloaded1.csv",error_bad_lines=False)

# %%
data = pd.DataFrame(data)
data

# %%
new_data = data[['cid','acvalue']]
new_data = new_data.dropna()

# %%
except_val=0
cid_value = new_data['cid'].to_list()
PIC50_value = (-np.log10(new_data['acvalue']*10**-6)).to_list()

PIC50_value = pd.DataFrame(PIC50_value,columns = ["y"])
y_data = PIC50_value

# %%
link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
sub_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
sub_link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
#the link with the required properties.
counter = 1

# %%
final_data_frame = pd.DataFrame()
x_data = pd.DataFrame()
cid_main_counter=0 
for x in cid_value:
    for iter_x in range(1000):
        try:    
            link = link_fixed
            link = link.replace("new_link_me",str(x))
            data1 = pd.read_csv(link)
            data1 = pd.DataFrame(data1)
            x_data = x_data.append(data1)
            data1 = 0
        except Exception as e:
            print("Exception Encountered as {} .! , Trying again !Iteration : {}".format(str(e),str(iter_x)))
            time.sleep(5)
            continue
        else:
            break
    
    else:
        print("Something Wrong with the Trials ! Restart The Algorithm!")

# %%
x_data_saved = x_data
x_data = x_data_saved
n_cid=[]

# %%
#20 sec Wait_time
for i in cid_value:
    max_tries = 10
    
    for iter_ in range(max_tries):
        try:
            #print("inside Try")
            #time.sleep(5)


            link = link_fixed
            link = link.replace("new_link_me",str(i))
            f_data = pd.read_csv(link) 
            
            sub_link = sub_link_fixed
            sub_link = sub_link.replace("replaceme",str(i))
            res = urllib.request.urlopen(sub_link)
            data_sub = res.read()
            data_sub = str(data_sub)
            data_sub = data_sub.replace("\\n",",") ; data_sub = data_sub.replace('b'," ")
            data_sub = data_sub.replace("'",""); data_sub = data_sub.replace(" ","")
            n_count_loop = 0 
            #time.sleep(3)
            #print("Going for the loop")
            for j in data_sub.split(","):
                n_count_loop+=1
                if(len(j)>1):
                    if(n_count_loop<=30):
                        n_cid.append(j)
                    else:
                        break
                else:
                    pass
             
        except:
            if(except_val>=15):
                break
            else:
                print("Re-Trying")
                except_val+=1
                #time.sleep(3)
                continue
        else:
            break
    '''else:
        print("Either The Network is Down(So no point in continuing ) , Or Some Uknown Error Spotted!. Refer Error Code , Continuing with Fetched Data")
        time.sleep(5)
        break'''
    
    
    '''final_data_frame = final_data_frame.append(f_data)
    cid_main_counter+=1
    final_data_frame.to_csv("fdf1.csv")
         
   

    print("Main cid no: {}".format(cid_main_counter))
    #time.sleep(5)'''
    



# %%
k_count = 0
t_count = 0
phase_count=0
f_data = pd.DataFrame()

# %%
for k in n_cid:
    k_count+=1
    t_count+=1
    phase_count+=1
    for sub_iter in range(10):
            try:
                if(t_count>=50):
                    #print("taking Rest")
                    time.sleep(15)
                    t_count=0
                    
                #print("Inside Try"+str(k_count))
            
                link = link_fixed
                link = link.replace("new_link_me",str(k))
                #print("replacement done"+str(k_count))
                f_data_df = pd.read_csv(link)
                f_data = f_data.append(f_data_df)
                print("Passed Without Exception"+str(k_count))
                break

            except Exception as e:
                print(str(e) + "Encountered , Please Wait :+ "+str(k_count))
                time.sleep(20)
                continue

    final_data_frame = final_data_frame.append(f_data)
    cid_main_counter+=1
    final_data_frame.to_csv("fdf1.csv")
    print("Data Fetching Continued"+ str(k_count))
    
    final_data_frame.drop_duplicates(inplace=True)


# %%
print(final_data_frame.shape)
final_data_frame.drop_duplicates(inplace=True)
print("The Final Training DataSet: X_data = {} , Y_data = {}".format(x_data.shape,y_data.shape))
print("The Final Data For Preditcion is : {}".format(final_data_frame.shape))

# %%
try:
    final_data_frame.reset_index(inplace=True)
    final_data_frame.drop('index',axis=1,inplace=True)
except:
    pass


# %%
#DON'T RE-RUN AT ANY COST
final_data_frame_reset = final_data_frame
x_data_reset = x_data
y_data_reset = y_data

# %%
#fghjk

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
RESETTING THE VALUES : TESTING PURPOSE
"""

# %%
'''#RUN TO RESET ONLY : TESTING CAUTION
final_data_frame = final_data_frame_reset
x_data = x_data_reset
y_data = y_data_reset'''

# %%
x_data = pd.DataFrame(x_data)

# %%
y_data = pd.DataFrame(y_data)

# %%
print("X :- "+str(x_data.shape)+" Y :- "+str(y_data.shape) )

# %%
train_file = pd.DataFrame()
x_data = x_data.astype("float64")
y_data = y_data.astype("float64")

# %%
x_data.head()

# %%
y_data.head()

# %%
x_data.head()

# %%
new = y_data['y'].to_list()

# %%
train_file  = x_data
train_file['y'] = new

# %%
cid_reg_list=train_file['CID'].to_list()
train_file.drop('CID',axis=1,inplace =True)


# %%
for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')
train_file

# %%
x = list(train_file.columns)
x = x[:-1]

x_ = train_file[x].values
y_ = train_file['y'].values

# %%
train_file.describe()
x_ = x_.astype("float64")
y_ = y_.astype("float64")

# %%
X_train, X_test, y_train, y_test = train_test_split(x_, y_, test_size=0.1, random_state=0)


# %%
"""
 
"""

# %%
"""
REGRESSION STARTS - 1,2,3 
"""

# %%
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing

# %%
#poly regression - degree 1
scaler = preprocessing.StandardScaler()
degree=1
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)

# %%
#R2 Score
y_pred = polyreg_scaled.predict(X_test)
main_r2 = r2_score(y_test, y_pred, multioutput='uniform_average')
max_ = main_r2 
max_

# %%
from itertools import combinations
comb_list = [[]]

# %%
#creating combinations , function
def sub(arr,r):
    global comb_list
    for i in r:
        comb = list(combinations(arr,i))
        comb_list.append(comb)
    return comb_list
newone = 0
newone = sub(x , [2,3,4,5,6])
del newone[0]

# %%
#Trying combinations for degree one and degree 2 only 

# %%
#degree 1 
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]
for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        scaler = preprocessing.StandardScaler()
        degree=1
        polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
        polyreg_scaled.fit(X_train,y_train)
        y_pred = polyreg_scaled.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]

# %%
#Extracting the best value so far using the regressor
sec_index = max_r2.index(max(max_r2))
fir_index = index_r2[sec_index]


# %%
if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
    maxi_r2 = max(max_r2)

# %%
print("R2 - degree1 = " + str(maxi_r2)+str("\n")+str("Features = ")+str(features))
reg_max_r2 = max(max_r2)

# %%
#We have the best degree 1 model
model_dict = {}
model_dict[1]=[maxi_r2,features]

# %%
for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')

# %%
#degree 2
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]
for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        scaler = preprocessing.StandardScaler()
        degree=2
        polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
        polyreg_scaled.fit(X_train,y_train)
        y_pred = polyreg_scaled.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]

# %%
#Extracting the best value so far using the regressor
sec_index = max_r2.index(max(max_r2))
fir_index = index_r2[sec_index]

# %%
if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
    maxi_r2 = max(max_r2)

# %%
print("R2 - degree2 = " + str(maxi_r2)+str("\n")+str("Features = ")+str(features))
reg_max_r2 = max(max_r2)

# %%
model_dict[2]=[maxi_r2,features]

# %%
for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')

# %%
model_dict

# %%
#degree 3
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]
for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        scaler = preprocessing.StandardScaler()
        degree=3
        polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
        polyreg_scaled.fit(X_train,y_train)
        y_pred = polyreg_scaled.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]

# %%
#Extracting the best value so far using the regressor
sec_index = max_r2.index(max(max_r2))
fir_index = index_r2[sec_index]

# %%
if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
    maxi_r2 = max(max_r2)

# %%
print("R2 - degree3 = " + str(maxi_r2)+str("\n")+str("Features = ")+str(features))
reg_max_r2 = max(max_r2)

# %%
model_dict[3]=[maxi_r2,features]

# %%
for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')

# %%
model_dict

# %%
top_cid_mix_df = pd.DataFrame()

# %%
"""
  
"""

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
DEGREE 1 : 
"""

# %%
model_dict[1]

# %%
features = list(model_dict[1][1])
features

# %%
x_trained =train_file[features].values
y_trained = train_file['y'].values
x_trained.shape,y_trained.shape

# %%
X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)

# %%
#poly regression - degree 1
scaler = preprocessing.StandardScaler()
degree=1
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)

# %%
#R2 Score
y_pred = polyreg_scaled.predict(X_test)
r2_score(y_test, y_pred, multioutput='uniform_average')

# %%
final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]

# %%
pred_data_cid = final_data_frameF2["CID"].to_list()

# %%
#features.append('CID')
new_features = features.append("CID")
features

# %%
final_data_frame1 = final_data_frameF2[features]
for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')
final_data_frame1.drop(columns='CID',inplace=True)
final_data_frame1.shape

# %%
final_pred_1 = polyreg_scaled.predict(final_data_frame1)

# %%
final_pred_1 = list(final_pred_1)
final_data_frame1['CID'] = pred_data_cid
final_data_frame1['y_'] = final_pred_1

# %%
saved_final_data_frame1 = final_data_frame1
sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)
final_larg= sorted_final_df.head(30)

# %%
final_larg = final_larg[['CID','y_']]
top_cid = final_larg["CID"].to_list()

# %%
print("The Top 30 Drug Leads Which are identified with PubChem cid's are (USING DEGREE 1 REGRESSION): ")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))

# %%
print("r2 = "+str(r2_score(y_test, y_pred, multioutput='uniform_average')))

# %%
top_cid_mix_df['Degree 1'] = top_cid

# %%
final_cid_1 = pd.DataFrame()
final_cid_1['CID'] = top_cid
final_cid_1.to_csv("final_cid_degree_1.csv")

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
DEGREE :2
"""

# %%
model_dict[2]

# %%
features = list(model_dict[2][1])
features

# %%
x_trained =train_file[features].values
y_trained = train_file['y'].values
x_trained.shape,y_trained.shape

# %%
X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)

# %%
#poly regression - degree 2
scaler = preprocessing.StandardScaler()
degree=2
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)

# %%
#R2 Score
y_pred = polyreg_scaled.predict(X_test)
r2_score(y_test, y_pred, multioutput='uniform_average')

# %%
final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]

# %%
pred_data_cid = final_data_frameF2["CID"].to_list()

# %%
#features.append('CID')
new_features = features.append("CID")
features

# %%
final_data_frame1 = final_data_frameF2[features]
for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')
final_data_frame1.drop(columns='CID',inplace=True)
final_data_frame1.shape

# %%
final_pred_1 = polyreg_scaled.predict(final_data_frame1)
final_pred_1 = list(final_pred_1)
final_data_frame1['CID'] = pred_data_cid
final_data_frame1['y_'] = final_pred_1

# %%
saved_final_data_frame1 = final_data_frame1
sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)
final_larg= sorted_final_df.head(30)

# %%
final_larg = final_larg[['CID','y_']]
top_cid = final_larg["CID"].to_list()

# %%
print("The Top 30 Drug Leads Which are identified with PubChem cid's are (USING DEGREE 2 REGRESSION): ")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))

# %%
top_cid_mix_df['Degree 2'] = top_cid

# %%
print("r2 = "+str(r2_score(y_test, y_pred, multioutput='uniform_average')))
final_cid_2 = pd.DataFrame()
final_cid_2['CID'] = top_cid
final_cid_2.to_csv("final_cid_degree_2.csv")

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
DEGREE :3
"""

# %%
model_dict[3]

# %%
features = list(model_dict[3][1])
features

# %%
x_trained =train_file[features].values
y_trained = train_file['y'].values
x_trained.shape,y_trained.shape

# %%
X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)

# %%
#poly regression - degree 3
scaler = preprocessing.StandardScaler()
degree=3
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)

# %%
#R2 Score
y_pred = polyreg_scaled.predict(X_test)
r2_score(y_test, y_pred, multioutput='uniform_average')

# %%
final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]

# %%
pred_data_cid = final_data_frameF2["CID"].to_list()

# %%
#features.append('CID')
new_features = features.append("CID")
features

# %%
final_data_frame1 = final_data_frameF2[features]
for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')
final_data_frame1.drop(columns='CID',inplace=True)
final_data_frame1.shape

# %%
final_pred_1 = polyreg_scaled.predict(final_data_frame1)
final_pred_1 = list(final_pred_1)
final_data_frame1['CID'] = pred_data_cid
final_data_frame1['y_'] = final_pred_1

# %%
saved_final_data_frame1 = final_data_frame1
sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)
final_larg= sorted_final_df.head(30)

# %%
final_larg = final_larg[['CID','y_']]
top_cid = final_larg["CID"].to_list()

# %%
print("The Top 30 Drug Leads Which are identified with PubChem cid's are (USING DEGREE 3 REGRESSION): ")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))

# %%
top_cid_mix_df['Degree 3'] = top_cid

# %%
print("r2 = "+str(r2_score(y_test, y_pred, multioutput='uniform_average')))
final_cid_2 = pd.DataFrame()
final_cid_2['CID'] = top_cid
final_cid_2.to_csv("final_cid_degree_2.csv")

# %%
"""
 THE END - COMPOUND FETCHING
"""

# %%
top_cid_mix_df

# %%
top_cid_mix_df.to_csv("TOP_CID_123.csv")

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
 
"""

# %%
"""
CODE FOR MOLECULAR DYNAMICS [RUN FOR ALL 3 SETS OF TOP 30 COMPOUNDS]
"""

# %%
import urllib.request
from os.path import expanduser
import os
home = expanduser("~")
assert (os.path.isdir(home+'/MGLTools-1.5.6')), "AutoDockTools not found!"
#set preparation pathways
prepare_protein_path = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -A bonds_hydrogens -U nphs_lps_waters_nonstdres -r'
prepare_ligand_path = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -A bonds_hydrogens -U nphs_lps -l'

# %%
l=[]
for itter in top_cid:
    l.append(itter)
for i in l:


# download SDF


    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/SDF'.format(str(i))
    urllib.request.urlretrieve(url, '{}'.format(str(i)+'.sdf'))
    
    
#     ligand preparation
    print("Preparing ligand {} for docking".format(i))
    os.system('obabel {} -O {} --gen3d --best'.format(str(i)+'.sdf', str(i)+'.mol2'))
    os.system(prepare_ligand_path + '{}'.format(str(i)+'.mol2'))
    
    

    
    print("Running docking procedure for 6w9c")
    os.system(
    'vina --config coordinatesandsetups6w9c.txt --receptor 6w9c.pdbqt --ligand {} --out {} --log {}'.format(
    str(i)+'.pdbqt', str(i)+'.6w9c.pdbqt', str(i)+'.6w9c.log')
)
    print("Docking to 6w9c protein completed. Files saved as {} as docked complex and {} as a logfile".format(
        str(i)+'.6x9c.pdbqt',str(i)+'.6w9c.log'))
    
    
    print("Running docking procedure for 1p9u")
    os.system(
    'vina --config coordinatesandsetups1p9u.txt --receptor 1p9u.pdbqt --ligand {} --out {} --log {}'.format(
    str(i)+'.pdbqt', str(i)+'.1p9u.pdbqt', str(i)+'.1p9u.log')
    )
    print("Docking to 1p9u protein completed. Files saved as {} as docked complex and {} as a logfile".format(
        str(i)+'.1p9u.pdbqt',str(i)+'.1p9u.log'))
    

            
    
    print(url)
    
    
    input("Press any key to continue In Silico modelling for next lead compound")


# %%
"""Now postprocessing the results"""
import os, pandas as pd, requests
from biopandas.pdb import PandasPdb

for file in os.listdir():
    if file.endswith('.6w9c.pdbqt'): 
        os.system("obabel {} -O {} -l 1".format(file, file.strip('.6w9c.pdbqt')+'.best.pdb')) 
        prot_df = PandasPdb().read_pdb('6w9c.pdb') 
        flex_df = PandasPdb().read_pdb(file.strip('.6w9c.pdbqt')+'.best.pdb') 
        flex_df.df['ATOM']['chain_id'].replace(to_replace='A', value = "X", inplace = True)  
        prot_df.df['ATOM'] = prot_df.df['ATOM'].append(flex_df.df['ATOM']) 
        prot_df.to_pdb(path='6w9c-{}.cplx.pdb'.format(file.strip('.6w9c.pdbqt')),records=['ATOM'],gz=False,
                     append_newline=True) 
    if file.endswith('.1p9u.pdbqt'):
        os.system("obabel {} -O {} -l 1".format(file, file.strip('.1p9u.pdbqt')+'.best.pdb')) 
        prot_df = PandasPdb().read_pdb('1p9u.pdb') # read protein file
        flex_df = PandasPdb().read_pdb(file.strip('.1p9u.pdbqt')+'.best.pdb') 
        flex_df.df['ATOM']['chain_id'].replace(to_replace='A', value = "X", inplace = True)   
        prot_df.df['ATOM'] = prot_df.df['ATOM'].append(flex_df.df['ATOM']) 
        prot_df.to_pdb(path='1p9u-{}.cplx.pdb'.format(file.strip('.1p9u.pdbqt')),records=['ATOM'],gz=False,
                     append_newline=True) 
                
        

# %%
"""Upload pdb and download img"""
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options
import time

# %%
options=Options()
# options.headless = True
plip = webdriver.Firefox(options=options)
plip.get("https://projects.biotec.tu-dresden.de/plip-web/plip")
print("Connection successful")

for cplx in os.listdir():
    if cplx.endswith('.cplx.pdb'):
        print("Processing file {}".format(cplx))
    
        select_pdb_input = plip.find_element_by_xpath("//*[@id='select-pdb-by-file']").click() 

        browse = plip.find_element_by_xpath(
            '/html/body/div[1]/div[2]/div/form/div[1]/div[1]/div[3]/input'
        ).send_keys(                                                       
            os.getcwd()+'/{}'.format(cplx)
        )

        send_file = plip.find_element_by_xpath("//*[@id='submit']").click() 
        time.sleep(10) 
        try:
            try:
                open_interactions_1 = plip.find_element_by_xpath('/html/body/div/div[2]/div/div[1]/h2[2]').click()
                open_interactions_2 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[2]/h3').click()
                open_interactions_3 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[2]/div/h4').click()
                pngs = plip.find_elements_by_xpath("//a[contains(@href,'.png')]")
                pymolsessions = plip.find_elements_by_xpath("//a[contains(@href,'.pse')]")
                
            except:
                open_interactions_1 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/h2[1]').click()
                open_interactions_2 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[1]/h3').click()
                open_interactions_3 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[1]/div/h4').click()
                pngs = plip.find_elements_by_xpath("//a[contains(@href,'.png')]")
                pymolsessions = plip.find_elements_by_xpath("//a[contains(@href,'.pse')]")

            for image in pngs:
                print(image.get_attribute("href"))
                output_image = requests.get(image.get_attribute("href"))
                open(
                    os.getcwd()+'/{}'.format(cplx+'.png'), 'wb'
                ).write(output_image.content)
                print("Image saved as {}".format(cplx+'.png'))

            for pysession in pymolsessions:
                print(pysession.get_attribute("href"))
                pse = requests.get(pysession.get_attribute("href"))
                open(
                    os.getcwd()+'/{}'.format(cplx+'.pse'), 'wb'
                ).write(pse.content)
                print("Pymol sessions saved as {}".format(cplx+'.pse'))
                  
            restart_plip = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/p[3]/a').click()
            time.sleep(5)
        except:
            print("No interactions found for {} or damaged structure".format(cplx))
            try:
                restart_plip = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/p[3]/a').click()
                time.sleep(5)
            except:
                restart_plip = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div[2]/p/a').click() 

# %%
"""
 
"""

# %%
"""
PROJECT ENDS
"""