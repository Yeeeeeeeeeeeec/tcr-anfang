#!/usr/bin/env python
# coding: utf-8

# In[42]:


import matplotlib.pyplot as plt
import pandas as pd
tcr = pd.read_csv('D:/check/top20Clean.txt',sep='\t')
#y = range(0,0.15,0.05)
y1 = tcr[tcr['sample']=='A12_CD40L_CMV_test_1']['freq']
y2 = tcr[tcr['sample']=='A12_CD40L_CMV_test_2']['freq']
x = [i for i in range(0,32,1)]
plt.xticks(x,tcr[tcr['sample']=='A12_CD40L_CMV_test_1']['cdr3'].tolist())
plt.xticks(rotation=80)
plt.title('A12_CD40L_CMV_test')
plt.xlabel('CDR3(alpha)')
plt.ylabel('cloneFraciton')
width = 0.4
plt.bar(x,y1,width = width,color = 'red')
for i in x:
    x[i] = x[i] + 0.4
plt.bar(x,y2,width = width,color = 'blue')
plt.legend(['A12_CD40L_CMV_test_1','A12_CD40L_CMV_test_2'])

plt.show()

