#!/usr/bin/env python
# coding: utf-8

# # Data Mining From PubChem

# In[139]:


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


# In[140]:


link1 = "https://pubchem.ncbi.nlm.nih.gov/sdq/cgi2rcgi.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22protacxn%22:%22notnull%22},{%22cid%22:%22notnull%22},{%22repacxn%22:%22P0C6X7%22}]},%22order%22:[%22activity,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PROTACXN_P0C6X7_bioactivity_protein%22}"


# In[141]:


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


# In[142]:


for i in range(2):
    try:
        req = requests.get(link1)
        url_content = req.content
        csv_file = open('downloaded1.csv', 'wb')
        csv_file.write(url_content)
        csv_file.close()
        print("Completed the Request")
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
            


# In[143]:


data = pd.read_csv("downloaded1.csv",error_bad_lines=False)

data_df = pd.DataFrame(data)

new_data = data[['cid','acvalue']]


new_data = new_data.dropna()


# In[144]:


#n_val=30
except_val=0


# In[145]:


cid_value = new_data['cid'].to_list()

PIC50_value = (-np.log10(new_data['acvalue']*10**-6)).to_list()


# In[146]:


PIC50_value = pd.DataFrame(PIC50_value,columns = ["y"])


# In[147]:


y_data = PIC50_value


# In[148]:


#len(y_data)


# In[149]:


link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
sub_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
sub_link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
#the link with the required properties.
counter = 1


# In[150]:


final_data_frame = pd.DataFrame()

# creating empty dataframe for attaching !

#final_data_frame


# In[151]:


x_data = pd.DataFrame()


cid_main_counter=0 


# 
# 

# In[152]:


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


# In[153]:


x_data_saved = x_data


# In[154]:


x_data = x_data_saved


# In[155]:


#x_data


# In[156]:


n_cid=[]


# In[157]:


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
    


# In[158]:


#len(n_cid)


# In[159]:


k_count = 0
t_count = 0
phase_count=0


# In[160]:


#n_cid_list = [44511706,25256828,44511706,25256828]
f_data = pd.DataFrame()


#replace


# In[162]:


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



    


# # Auto-QSAR And Drug-Lead Generation

# In[163]:


final_data_frame.shape


# In[164]:



#condition
#final_data_frame= pd.read_csv("Final_datacsv_test1.csv")


# In[165]:


#condition
#final_data_frame.drop(columns="Unnamed: 0",inplace = True)


# In[166]:


final_data_frame.drop_duplicates(inplace=True)


# In[167]:


print("The Final Training DataSet: X_data = {} , Y_data = {}".format(x_data.shape,y_data.shape))


# In[168]:


print("The Final Data For Preditcion is : {}".format(final_data_frame.shape))


# In[169]:


final_data_frame


# In[170]:


x_data = pd.DataFrame(x_data)
x_data.shape


# In[171]:


y_data = pd.DataFrame(y_data)
y_data.shape


# In[172]:


train_file = pd.DataFrame()
train_file


# In[173]:


x_data = x_data.astype("float64")
y_data = y_data.astype("float64")


# In[174]:


x_data.head()


# In[175]:


y_data.head()


# In[176]:


new = y_data['y'].to_list()


# In[177]:


train_file  = x_data
train_file['y'] = new


# In[178]:


train_file
cid_reg_list=train_file['CID'].to_list()
train_file.drop('CID',axis=1,inplace =True)


# In[179]:


train_file


# In[180]:


for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')


# In[ ]:





# In[181]:


len(cid_reg_list)


# In[ ]:





# In[ ]:





# In[182]:


train_file['y']


# In[183]:


x = list(train_file.columns)
x = x[:-1]

x_ = train_file[x].values
y_ = train_file['y'].values


# In[184]:


train_file.describe()


# In[ ]:





# In[185]:


x_ = x_.astype("float64")
y_ = y_.astype("float64")


# In[186]:


y_.shape


# In[187]:


x_.shape


# In[188]:


X_train, X_test, y_train, y_test = train_test_split(x_, y_, test_size=0.1, random_state=0)


# In[189]:


print("Regression Starts")


# In[190]:


regressor = LinearRegression()  
regressor.fit(X_train, y_train)


# In[191]:



coef_dict = dict(zip(x, regressor.coef_))

y_pred = regressor.predict(X_test)


# In[192]:


main_r2 = r2_score(y_test, y_pred, multioutput='uniform_average')
max_ = main_r2 



from itertools import combinations
comb_list = [[]]


def sub(arr,r):
    global comb_list
    for i in r:
        comb = list(combinations(arr,i))
        comb_list.append(comb)
    return comb_list


# In[ ]:





# In[193]:


newone =0

newone = sub(x , [2,3,4,5,6])
del newone[0]

#index_count = 0
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]


# In[194]:



for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        regressor = LinearRegression()  
        regressor.fit(X_train, y_train)
        y_pred = regressor.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]


# In[195]:



sec_index = max_r2.index(max(max_r2))


# In[196]:


#max_r2


# In[197]:


fir_index = index_r2[sec_index]


# In[198]:



if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
maxi_r2 = max(max_r2)


# In[199]:


#features


# In[200]:


reg_max_r2 = max(max_r2)
#reg_max_r2


# In[ ]:





# New Regression Part Starts [processing the best output to print out values/]

# In[201]:


#train_file


# In[202]:


for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')


# In[203]:


#train_file


# In[204]:


features  = list(features)


# In[205]:


x_trained =train_file[features].values


# In[206]:


x_trained.shape


# In[207]:


y_trained = train_file['y'].values


# In[208]:


y_trained.shape


# In[209]:


X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)


# In[210]:


regressor = LinearRegression()  
regressor.fit(X_train, y_train)


# In[211]:


coef_dict = dict(zip(features, regressor.coef_))
coef_dict = pd.DataFrame(coef_dict , index = [0])


# In[212]:


y_pred = regressor.predict(X_test)


# In[ ]:





# In[213]:


r2_score(y_test, y_pred, multioutput='uniform_average')


# In[214]:


features


# In[215]:


saved_features = features


# In[216]:


coef_dict

final_data_frame


# In[217]:


final_data_frame

final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]

final_data_frameF2


# In[218]:


pred_data_cid = final_data_frameF2["CID"].to_list()


# In[219]:


features.append('CID')


# In[220]:


new_features = features.append("CID")
features


# In[221]:


final_data_frame1 = final_data_frameF2[features]


# In[222]:


final_data_frame1


# In[223]:


final_data_frame1.shape


# In[224]:


for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')


# In[225]:


final_data_frame1.drop(columns='CID',inplace=True)


# In[226]:


final_data_frame1


# In[227]:


final_pred = regressor.predict(final_data_frame1)


# In[228]:


final_pred = list(final_pred)


# In[229]:


final_data_frame1['CID'] = pred_data_cid


# In[230]:


final_data_frame1['y_'] = final_pred


# In[231]:


final_data_frame1


# In[232]:


saved_final_data_frame1 = final_data_frame1


# In[233]:


sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)


# In[234]:


final_larg= sorted_final_df.head(30)


# In[235]:


final_larg


# In[236]:


final_larg = final_larg[['CID','y_']]


# In[237]:


top_cid = final_larg["CID"].to_list()


# In[238]:


for i in range(0,2):
    try:
        os.remove('downloaded1.csv')
        print("Deleted old File")
        break
    except Exception as e:
        print("No file ")
        break
    else:
        break
for i in range(0,2):
    try:
        os.remove('fdf1.csv')
        print("Deleted old File")
        break
    except Exception as e:
        print("No file ")
        break
    else:
        break


# In[239]:


print("The Top 30 Drug Leads Which are identified with PubChem cid's are: ")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))


# In[240]:


print("r2 = "+str(r2_score(y_test, y_pred, multioutput='uniform_average')))


# In[241]:


print("Intercept = "+str(regressor.intercept_))


# In[242]:


print("coefficients are :- ")


# In[243]:


print(coef_dict)


# In[246]:
# Automated In Silico modelling 

import urllib.request
from os.path import expanduser
import os
home = expanduser("~")
assert (os.path.isdir(home+'/MGLTools-1.5.6')), "AutoDockTools not found!"
#set preparation pathways
prepare_protein_path = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -A bonds_hydrogens -U nphs_lps_waters_nonstdres -r'
prepare_ligand_path = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -A bonds_hydrogens -U nphs_lps -l'



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


# In[248]:


"""Now postprocessing the results"""
import os, pandas as pd, requests
from biopandas.pdb import PandasPdb


# In[249]:


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
                
        


# In[253]:


"""Upload pdb and download img"""
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options
import time


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


# In[ ]:




