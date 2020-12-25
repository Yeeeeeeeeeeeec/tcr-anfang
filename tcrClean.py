#!/usr/bin/env python
# coding: utf-8

# In[1]:


#list.append()会修改list，DataFrame新建副本 需重新赋值！
import pandas as pd
import os
import matplotlib.pyplot as plt
import re


# In[2]:


#遍历当前文件夹，抓取以‘clonotypes.TRA.txt'为结尾的文本文档
def pathList(base):
    pathlist = []
    for root,ds,fs in os.walk(base):
        for f in fs:
            if f.endswith('clonotypes.TRA.txt'):
                fullname = os.path.join(root,f)
                pathlist.append(fullname)
    return pathlist


# In[6]:


#将表头设置为['cdr3','sample','freq','dna','group']
def setFormat(pathlist):
    #加入sample和group列
    tcrlist = []
    for path in pathlist:
        samplename = path[(path.rfind('/')+1) : path.rfind('.clonotypes')]
        samplecsv = pd.read_csv(path,sep='\t')
        samplecsv['sample'] = samplename
        groupname = samplename[0:samplename.rfind('_')]
        samplecsv['group'] = groupname
        tcr = samplecsv[['aaSeqCDR3','sample','cloneFraction','nSeqCDR3','group']]
        tcr.columns = ['cdr3','sample','freq','dna','group']
        tcrlist = tcrlist.append(tcr)
    return tcrlist


# In[7]:


def splitTop20(tcrlist):
    tcralltop20 = pd.DataFrame()
    #抓取所有top20及其在其他样本的dna
    for tcr in tcrlist:
        top20 = tcr.head(20)
        for tcr in tcrlist:
            tcralltop20 = tcralltop20.append(tcr.loc[tcr['dna'].isin(top20['dna'])])
    #去除重复项
    tcralltop20.drop_duplicates(inplace=True)
    tcralltop20.sort_values(by='group',inplace=True)
    #tcralltop20.reset_index(drop=True)
    grouplist = tcralltop20['group'].unique()
    #按group分割tcr并放进list
    listbygroup = []
    for gp in grouplist:
        tcrbygroup = tcralltop20.loc[tcralltop20['group']==gp]
        tcrbygroup.reset_index(drop=True)
        listbygroup = listbygroup.append(tcrbygroup)
    return listbygroup


# In[8]:


def getLastDig(splist):
    return re.findall('(\d+)',splist)[-1]
def sortTcr(listbygroup):
    tcrclean = pd.DataFrame()
    for tcr in listbygroup:
        dnalist = tcr['dna'].unique()
        splist = tcr['sample'].unique()
        splist.sort(key=getLastDig)

        '''for dna in dnalist:
            listbydna.append(tcr.loc[tcr['dna']==dna])
        #将非0排序标准样本按sample分割并给予rank
        for dnadf in listbydna:'''
        
        rankdf = pd.DataFrame()
        
        for sp in splist():
            while not dnalist:
                temp = tcr.loc[(tcr['sample']==sp) & (tcr['freq']!=0) & (tcr['dna'].isin(dnalist))]
                temp.sort_values(by='freq',ascending=False,inplace=True)
                rankdf = rankdf.append(temp)
                for dna in temp['dna']:
                    dnalist.remove(dna)
        rankdf.reset_index(drop=True)
        rankdf['rank'] = rankdf.index #.tolist()
        dnalist3 = tcr['dna'].unique()
        
        for dna in dnalist3:
            id = rankdf.loc[rankdf['dna']==dna].index
            tcr.loc[tcr['dna']==dna]['rank'] = rankdf['rank'][id]
            
        tcr.sort_values(by=['rank','sample'],inplace=True)    
        
        tcrclean = tcrclean.append(tcr)
        tcrclean.to_csv('tcrTop20Clean.txt',sep='\t')
    return tcrclean    


# In[9]:


def tcrMain():
    base = '/data/Data/chenye/tcr/output_A/'
    pathlist = pathList(base)
    tcrlist = setFormat(pathlist)
    listbygroup = splitTop20(tcrlist)
    sortTcr(listbygroup)

