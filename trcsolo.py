#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
pd.set_option('display.unicode.east_asian_width',True)
tcr1 = pd.read_csv('D:/check/A12_CD40L_CMV_test_1.clonotypes.TRA.txt',sep='\t')
tcr2 = pd.read_csv('D:/check/A12_CD40L_CMV_test_2.clonotypes.TRA.txt',sep='\t')
sample1 = 'A12_CD40L_CMV_test_1'
sample2 = 'A12_CD40L_CMV_test_2'
tcr1['sample'] = sample1
tcr2['sample'] = sample2


# In[3]:


def set_group(tcr):
    sp = tcr['sample'][0]
    gp = sp[0:sp.rfind('_')]
    tcr['group'] = gp


# In[19]:


def clean(tcr):
    tcr = tcr[['aaSeqCDR3','sample','cloneFraction','nSeqCDR3','group']]
    tcr.columns = ['cdr3','sample','freq','dna','group']
    return tcr 
    


# In[68]:


def top20unsort(tcr1,tcr2):
    top20 = pd.DataFrame()
    tcr120 = tcr1.head(20)
    tcr220 = tcr2.head(20)
    '''tcr2_in1 = tcr1.loc[tcr2['dna'].isin(top20['dna'])]
    top20 = top20.append(tcr2_in1)
    tcr2top20 = tcr2.head(20)
    tcr2_no1 = tcr2top20.loc[~tcr2['dna'].isin( top20['dna'])]
    top20 = top20.append(tcr2_no1)'''
    tcr2in1 = tcr2.loc[tcr2['dna'].isin(tcr120['dna'])]
    tcr1in2 = tcr1.loc[tcr1['dna'].isin(tcr220['dna'])]
    top20 = tcr2in1.append(tcr1in2)
    top20 = top20.append(tcr120)
    top20 = top20.append(tcr220)
    top20 = top20.drop_duplicates()
    return top20.reset_index(drop=True)
        


# In[135]:


def top20sort(tcr):
    top20 = pd.DataFrame()
    sortby2 = pd.DataFrame()
    top20 = top20.append(tcr.loc[(tcr['freq']!=0)&(tcr['sample']==sample1)])
    top20 = top20.sort_values(by='freq',ascending=False)
    top20 = top20.reset_index(drop=True)
    for index,row in tcr.iterrows():
        if (row['sample']==sample2) & (row['dna'] in top20['dna'].tolist()):
            top20cut = pd.DataFrame()
            sameIndex = top20[top20['dna']==row['dna']].index.tolist()[0]
            top20cut = top20.iloc[0:(sameIndex+1)] 
            top20cut = top20cut.append(row)
            top20cut = top20cut.append(top20.iloc[(sameIndex+1)::])
            top20 = top20cut.reset_index(drop=True)
        if (row['sample']==sample2) & (row['dna'] not in top20['dna'].tolist()):
            sortby2 = sortby2.append(row)
    if not sortby2.empty:
        sortby2.sort_values(by='freq',ascending=False)
        top20 = top20.append(sortby2)                              
    for index,row in tcr.iterrows():    
        if (row['sample']==sample1) & (row['freq']==0):
            sameIndex = top20[top20['dna']==row['dna']].index.tolist[0]
            top20cut = top20.iloc[0:sameIndex]
            top20cut = top20cut.append(row)
            top20cut = top20cut.append(top20.iloc[(sameIndex+1)::])
            top20 = top20cut     
    return top20        
    
    
        
        


# In[136]:


set_group(tcr1)
set_group(tcr2)
t1Clean = clean(tcr1)
t2Clean = clean(tcr2)
top20Unsort = top20unsort(t1Clean,t2Clean)
top20Clean = top20sort(top20Unsort)
top20Clean.to_csv('D:/check/top20Clean.txt',sep='\t')


# In[137]:


top20Clean


# In[ ]:


import matplotlib.pyplot as plt

