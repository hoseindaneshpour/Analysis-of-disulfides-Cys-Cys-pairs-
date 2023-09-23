#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import AlignIO
import numpy as np
import pandas as pd


# Analysis of disulfides (Cys-Cys pairs)

# In[2]:


alignment = AlignIO.read("data/UP hits 201-600 64 re-realigned.fasta", "fasta")
align_array = np.array(alignment) 
columns_of_interest = [363, 390, 451,529,560,584,605,709,730,742]
id_positions=[]
for i in range (0 , len(alignment)):
    count1= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[0]] == "-") #gaps before the first clmn of intrst
    count2= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[1]] == "-")
    count3= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[2]] == "-")
    count4= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[3]] == "-")
    count5= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[4]] == "-")
    count6= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[5]] == "-")
    count7= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[6]] == "-")
    count8= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[7]] == "-")
    count9= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[8]] == "-")
    count10= np.count_nonzero(align_array[i:i+1, 0:columns_of_interest[9]] == "-")
    lst1= [alignment[i].id.split("|")[1],columns_of_interest[0]-count1,columns_of_interest[1]-count2, columns_of_interest[2]-count3,
           columns_of_interest[3]-count4,columns_of_interest[4]-count5,
           columns_of_interest[5]-count6,
           columns_of_interest[6]-count7,
           columns_of_interest[7]-count8,
           columns_of_interest[8]-count9,
           columns_of_interest[9]-count10]
    id_positions.append(lst1)

    letters=[]
for record in alignment:
    amino_acids = ([record.seq[column-1] for column in columns_of_interest])
    lst2=amino_acids
    letters.append(lst2)

for j in range (0 ,len(id_positions)):
    print(f"{id_positions[j]} : {letters[j]}")


# In[3]:


counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results = pd.DataFrame(list(zip(values,counts)),columns=["value","count"])
df_results


# In[25]:


mylist=[]
for j in range (0 ,len(id_positions)):
    mylist.append(id_positions[j]+letters[j])
df_result = pd.DataFrame(mylist)
df_result
df_result.rename(columns={df_result.columns[11]: 'C1',df_result.columns[12]: 'C2',df_result.columns[13]: 'C3',
                          df_result.columns[14]: 'C4',df_result.columns[15]: 'C5', df_result.columns[16]: 'C6',
                          df_result.columns[17]: 'C7',df_result.columns[18]: 'C8',df_result.columns[19]: 'C9',
                          df_result.columns[20]: 'C10'}, inplace=True)
pd.set_option('display.max_columns', None)
df_result


df_result.to_html('result.html', index=False)


# In[5]:


results = []
for i in range(0, len(id_positions)):
    if df_result.loc[i, 'C1'].find('C') != -1 and df_result.loc[i, 'C7'].find('C') != -1:
        results.append([df_result.loc[i, 0], 1])
    elif df_result.loc[i, 'C1'].find('C') != -1 or df_result.loc[i, 'C7'].find('C') != -1:
        results.append([df_result.loc[i, 0], 0.5])
    else:
        results.append([df_result.loc[i, 0], 0])

df_results_DS_A = pd.DataFrame(results, columns=['ID', 'DS_A'])


# In[6]:


results = []
for i in range(0, len(id_positions)):
    if df_result.loc[i, 'C2'].find('C') != -1 and df_result.loc[i, 'C8'].find('C') != -1:
        results.append([df_result.loc[i, 0], 1])
    elif df_result.loc[i, 'C2'].find('C') != -1 or df_result.loc[i, 'C8'].find('C') != -1:
        results.append([df_result.loc[i, 0], 0.5])
    else:
        results.append([df_result.loc[i, 0], 0])

df_results_DS_B = pd.DataFrame(results, columns=['ID', 'DS_B'])


# In[7]:


results = []
for i in range(0, len(id_positions)):
    if df_result.loc[i, 'C3'].find('C') != -1 and df_result.loc[i, 'C6'].find('C') != -1:
        results.append([df_result.loc[i, 0], 1])
    elif df_result.loc[i, 'C3'].find('C') != -1 or df_result.loc[i, 'C6'].find('C') != -1:
        results.append([df_result.loc[i, 0], 0.5])
    else:
        results.append([df_result.loc[i, 0], 0])

df_results_DS_C = pd.DataFrame(results, columns=['ID', 'DS_C'])


# In[8]:


results = []
for i in range(0, len(id_positions)):
    if df_result.loc[i, 'C4'].find('C') != -1 and df_result.loc[i, 'C9'].find('C') != -1:
        results.append([df_result.loc[i, 0], 1])
    elif df_result.loc[i, 'C4'].find('C') != -1 or df_result.loc[i, 'C9'].find('C') != -1:
        results.append([df_result.loc[i, 0], 0.5])
    else:
        results.append([df_result.loc[i, 0], 0])

df_results_DS_D = pd.DataFrame(results, columns=['ID', 'DS_D'])


# In[9]:


results = []
for i in range(0, len(id_positions)):
    if df_result.loc[i, 'C5'].find('C') != -1 and df_result.loc[i, 'C10'].find('C') != -1:
        results.append([df_result.loc[i, 0], 1])
    elif df_result.loc[i, 'C5'].find('C') != -1 or df_result.loc[i, 'C10'].find('C') != -1:
        results.append([df_result.loc[i, 0], 0.5])
    else:
        results.append([df_result.loc[i, 0], 0])

df_results_DS_E = pd.DataFrame(results, columns=['ID', 'DS_E'])


# In[10]:


df_concat = pd.concat([df_results_DS_A, df_results_DS_B.drop('ID', axis=1), df_results_DS_C.drop('ID', axis=1), df_results_DS_D.drop('ID', axis=1), df_results_DS_E.drop('ID', axis=1)], axis=1)
df_concat


# In[27]:


df_concat['Row Sum'] = df_concat[['DS_A', 'DS_B', 'DS_C', 'DS_D', 'DS_E']].sum(axis=1)
df_concat.loc['Column Sum'] = df_concat[['DS_A', 'DS_B', 'DS_C', 'DS_D', 'DS_E']].sum()
df_concat

df_concat.to_html('result2.html', index=False)


# In[ ]:




