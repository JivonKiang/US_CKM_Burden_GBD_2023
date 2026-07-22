from pathlib import Path
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT=Path(os.environ.get('CKM_PROJECT_ROOT', Path.cwd())); OUT=ROOT/'Figures_Tables'; DATA=ROOT/'data'
CAUSES=[('Ischemic heart disease','IHD'),('Stroke','Stroke'),('Atrial fibrillation and flutter','AF/flutter'),('Lower extremity peripheral arterial disease','PAD'),('Diabetes Mellitus','Diabetes'),('Chronic Kidney Disease','CKD')]
C={'IHD':'#153B5B','Stroke':'#3F78B3','AF/flutter':'#758EAE','PAD':'#A0957F','Diabetes':'#CE8B12','CKD':'#B21D35'}
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':8,'axes.linewidth':.8,'axes.spines.right':False,'axes.spines.top':False})
def save(f,n):
 for e,d in [('png',300),('pdf',300),('svg',300),('tiff',600)]:f.savefig(OUT/f'{n}.{e}',dpi=d,bbox_inches='tight')
 plt.close(f)
def lab(ax,s):ax.text(-.1,1.04,s,transform=ax.transAxes,fontweight='bold',fontsize=10)
def load(m):
 z=[]
 for fn,c in CAUSES:
  x=pd.read_csv(DATA/m/f'{fn}.csv');x['component']=c;z.append(x)
 return pd.concat(z)
def f1():
 j=pd.read_csv(DATA/'anchors/JAMA_NHANES_CKM_2011_2020_overall.csv').query("stage != 'Advanced stages 3 or 4'")
 w=pd.read_csv(DATA/'anchors/WHO_US_CKM_risk_indicators_2011_2024.csv'); t=load('DALY'); p=load('Prevalence')
 fig,axs=plt.subplots(2,2,figsize=(7.2,5.8)); a,b,c,d=axs.flat
 a.bar(j.stage,j.prevalence_percent,color=['#17807E','#5D94C8','#3F78B3','#CE8B12','#B21D35']);a.errorbar(range(5),j.prevalence_percent,yerr=[j.prevalence_percent-j.ci_lower_percent,j.ci_upper_percent-j.prevalence_percent],fmt='none',color='#222',capsize=2);a.set_ylabel('Prevalence (%)');a.set_ylim(0,60);lab(a,'a')
 for ax,x,ttl in [(b,p,'2023 prevalence'),(c,t,'2023 DALYs')]:
  q=x.query('Year==2023').sort_values('Value');ax.barh(q.component,q.Value,color=[C[i] for i in q.component]);ax.set_title(ttl,loc='left',fontweight='bold');ax.set_xlabel('Rate per 100,000');lab(ax,'b' if ax is b else 'c')
 for code,col,name in [('NCD_BMI_30A','#CE8B12','Obesity'),('NCD_DIABETES_PREVALENCE_AGESTD','#17807E','Diabetes'),('BP_04','#B21D35','Raised BP')]:
  q=w.query('indicator_code==@code').sort_values('year');d.plot(q.year,q.value_percent,color=col,lw=1.7,label=name)
 d.legend(fontsize=7);d.set_ylabel('%');d.set_title('WHO exposure trends',loc='left',fontweight='bold');lab(d,'d')
 fig.text(.5,.01,'JAMA stages, GBD components, and WHO exposures retain source-specific constructs and denominators.',ha='center',fontsize=7,color='#555');fig.tight_layout(rect=(0,0.04,1,1));save(fig,'Figure1_US_CKM_source_map')
def f2():
 fig,axs=plt.subplots(2,3,figsize=(7.2,5.4)); mets=[('Prevalence','Prevalence'),('Incidence','Incidence'),('Deaths','Deaths'),('DALY','DALYs')]
 for ax,(m,y),s in zip(axs.flat[:4],mets,'abcd'):
  z=load(m)
  for k in C:
   q=z[z.component==k];ax.plot(q.Year,q.Value,color=C[k],lw=1.35,label=k)
  ax.set_yscale('log');ax.set_ylabel(y+' per 100,000');ax.grid(axis='y',lw=.4,color='#d9dde3');lab(ax,s)
 axs[0,0].legend(ncol=3,fontsize=6,loc='upper left',bbox_to_anchor=(0,1.35))
 tr=pd.read_csv(OUT/'Table_GBD_US_CKM_component_trends_1990_2023.csv');mat=np.array([[tr[(tr.component==k)&tr.measure.str.contains(m,case=False)].percent_change_1990_2023.iloc[0] for m in ['Prevalence','Incidence','Deaths','DALYs']] for k in C]);ax=axs[1,1];im=ax.imshow(mat,cmap='RdBu_r',norm=TwoSlopeNorm(vcenter=0,vmin=-65,vmax=200));ax.set_xticks(range(4));ax.set_xticklabels(['Prev.','Inc.','Deaths','DALYs'],rotation=32,ha='right',rotation_mode='anchor',fontsize=6.5);ax.tick_params(axis='x',pad=2);ax.set_yticks(range(6),list(C));lab(ax,'e');fig.colorbar(im,ax=ax,fraction=.04,pad=.12)
 ax=axs[1,2];da=tr[tr.measure.str.contains('DALYs')].set_index('component').loc[list(C)].sort_values('rate_2023');ax.barh(da.index,da.rate_2023,color=[C[x] for x in da.index]);ax.set_xlabel('2023 DALYs/100,000');lab(ax,'f')
 fig.tight_layout();save(fig,'Figure2_US_CKM_GBD_trajectories_and_change')
def f3():
 fig,axs=plt.subplots(2,3,figsize=(7.2,4.7),sharex=True,sharey=True)
 for ax,(fn,k),s in zip(axs.flat,CAUSES,'abcdef'):
  for m,ls in [('Prevalence','-'),('Incidence','--'),('Deaths',':'),('DALY','-.')]:
   q=load(m).query('component==@k').sort_values('Year');ax.plot(q.Year,100*q.Value/q.Value.iloc[0],color=C[k],ls=ls,lw=1.4)
  ax.axhline(100,color='#999',lw=.6);ax.set_title(k,loc='left',fontweight='bold');ax.set_ylim(0,310);lab(ax,s)
 axs[1,0].set_ylabel('Index (1990=100)');axs[1,0].set_xlabel('Year');axs[1,1].set_xlabel('Year');axs[1,2].set_xlabel('Year')
 fig.legend(['Prevalence','Incidence','Deaths','DALYs'],ncol=4,loc='lower center',bbox_to_anchor=(.5,-.04),fontsize=7);fig.tight_layout();save(fig,'Figure3_US_CKM_component_index_trends')
def f4():
 fig,axs=plt.subplots(2,2,figsize=(7.2,5.3));
 for ax,(m,t),s in zip(axs.flat,[('Prevalence','Prevalence'),('Incidence','Incidence'),('Deaths','Deaths'),('DALY','DALYs')],'abcd'):
  q=load(m).query('Year==2023').sort_values('Value');v=q.Value.values;ax.barh(q.component,v,color=[C[x] for x in q.component]);ax.errorbar(v,range(6),xerr=[v-q['Lower bound'].values,q['Upper bound'].values-v],fmt='none',color='#222',capsize=2,lw=.7);ax.set_title(t+' (2023)',loc='left',fontweight='bold');ax.set_xlabel('Rate per 100,000');lab(ax,s)
 fig.tight_layout();save(fig,'Figure4_US_CKM_2023_absolute_burden')
def f5():
 s=pd.read_csv(DATA/'anchors/JAMA_NHANES_CKM_2011_2020_subgroups.csv');fig=plt.figure(figsize=(7.2,6.4));grid=fig.add_gridspec(2,3,height_ratios=[4.5,1.35],hspace=.52,wspace=.85);axs=[fig.add_subplot(grid[0,i]) for i in range(3)]
 stages=['stage0_percent','stage1_percent','stage2_percent','stage3_percent','stage4_percent'];colors=['#17807E','#5D94C8','#3F78B3','#CE8B12','#B21D35']
 for ax,(typ,title),letter in zip(axs,[('age','Age'),('sex','Sex'),('race_ethnicity','Race/ethnicity')],'abc'):
  q=s[s.subgroup_type==typ].copy();y=np.arange(len(q));left=np.zeros(len(q))
  for st,col in zip(stages,colors):ax.barh(y,q[st],left=left,color=col,height=.62);left+=q[st].to_numpy()
  ax.set_yticks(y,q.subgroup,fontsize=6.5);ax.invert_yaxis();ax.set_title(title,fontweight='bold',pad=14);ax.set_xlim(0,100);ax.set_xlabel('%');ax.text(-.12,1.05,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
 axs[0].legend(['Stage 0','Stage 1','Stage 2','Stage 3','Stage 4'],fontsize=6,ncol=5,loc='upper left',bbox_to_anchor=(0,1.27))
 ax=fig.add_subplot(grid[1,:]);ax.axis('off');ax.text(-.02,1.03,'d',transform=ax.transAxes,fontweight='bold',fontsize=10);ax.text(.5,1.04,'Evidence use: complementary, not interchangeable',ha='center',va='bottom',fontweight='bold',fontsize=8,color='#153B5B',transform=ax.transAxes)
 cards=[('JAMA/NHANES','Stage 0-4 are mutually\nexclusive persons','Direct population anchor','#e5f1ef'),('GBD','Components can coexist\nin one person','Use for burden trajectories','#e7edf5'),('WHO','Exposure indicators use\nsource-specific denominators','Use for upstream context','#f7eee0')]
 for i,(head,body,use,face) in enumerate(cards):
  x=.17+.33*i;ax.text(x,.58,head,ha='center',va='center',fontweight='bold',fontsize=7.5,color='#153B5B',transform=ax.transAxes);ax.text(x,.33,body,ha='center',va='center',fontsize=7,bbox=dict(boxstyle='round,pad=.35',fc=face,ec='#c9c9c9'),transform=ax.transAxes);ax.text(x,.05,use,ha='center',va='bottom',fontsize=6.5,color='#444',transform=ax.transAxes)
  if i<2:ax.annotate('',xy=(x+.13,.34),xytext=(x+.20,.34),xycoords=ax.transAxes,arrowprops=dict(arrowstyle='->',color='#777',lw=.8))
 ax.text(.5,-.18,'Do not convert, sum, ratio, or statistically validate these constructs as one common CKM-stage metric.',ha='center',va='top',fontsize=6.7,color='#B21D35',fontweight='bold',transform=ax.transAxes)
 fig.subplots_adjust(top=.84,bottom=.10,left=.09,right=.97);save(fig,'Figure5_US_CKM_NHANES_stage_subgroups')
def f6():
 fig,ax=plt.subplots(figsize=(7.2,3.8));ax.axis('off');cols=['JAMA/NHANES','GBD','WHO'];rows=[['Stage 0-4\nmutually exclusive','IHD, stroke, AF, PAD,\ndiabetes, CKD\nnonadditive','Obesity, diabetes, BP\nrisk exposures'],['Individual-level\nstage distribution','Long-term burden\nprevalence/incidence/deaths/DALYs','Upstream risk\ncontext'],['Direct anchor','Component trajectories','Independent exposure context']];
 for i,c in enumerate(cols):ax.text(.17+.33*i,.91,c,ha='center',va='center',fontweight='bold',color='#153B5B')
 for r,row in enumerate(rows):
  for i,txt in enumerate(row):ax.text(.17+.33*i,.70-.28*r,txt,ha='center',va='center',bbox=dict(boxstyle='round,pad=.5',fc=['#e5f1ef','#e7edf5','#f7eee0'][i],ec='#cccccc'))
 ax.text(.5,.06,'No conversion, summation, ratio, or validation statistic is permitted across the three construct families.',ha='center',color='#B21D35',fontweight='bold');save(fig,'Figure6_US_CKM_cross_source_interpretation')
if __name__=='__main__':
 f1();f2();f3();f4();f5();print('Final Figure 1-5 suite created')
