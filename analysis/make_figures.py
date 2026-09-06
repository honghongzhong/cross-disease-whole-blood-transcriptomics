"""Draw publication figures from the analysis result tables."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'lines.linewidth':1.1,'savefig.facecolor':'white'})
E=['hedges_g','log2fc','shrunken_log2fc'];N=["Hedges' g",'log2FC','Shrunken log2FC'];C=['#0072B2','#D55E00','#009E73']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();A=args.results;F=args.output;F.mkdir(parents=True,exist_ok=True)
 def read(n):return pd.read_csv(A/n,sep='\t')
 def save(fig,n):
  fig.savefig(F/(n+'.pdf'),bbox_inches='tight',pad_inches=.035);fig.savefig(F/(n+'.png'),dpi=450,bbox_inches='tight',pad_inches=.035);plt.close(fig)
 fig,ax=plt.subplots(figsize=(6.69,3.3));ax.set(xlim=(0,10),ylim=(0,5.2));ax.axis('off')
 boxes=[(.15,3.8,2.75,1.0,'Within-study effects\n48 disease–control contrasts\n43 studies · 12 categories'),(3.62,3.8,2.75,1.,'Learn shared directions\n36 comparisons · 6144 genes\nRank from training data'),(7.1,3.8,2.75,1.,'Represent other studies\n12 held-out comparisons\nFixed directions, fitted scores'),(.15,1.4,2.75,1.25,'Structure beyond chance?\n499 gene-identity permutations\nSelect consecutive directions'),(3.62,1.4,2.75,1.25,'Is the structure shared?\nExclude each disease category\nReselect rank and refit'),(7.1,1.4,2.75,1.25,'What does it contain?\nHallmark enrichment\nGroup and heme sensitivity')]
 for x,y,w,h,t in boxes:
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.07',facecolor='#F1F5F8',edgecolor='#56768C',lw=.8));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=6.8,linespacing=1.7)
 for x1,y1,x2,y2 in [(3.02,4.3,3.5,4.3),(6.5,4.3,6.98,4.3),(1.52,3.65,1.52,2.81),(5,3.65,5,2.81),(8.47,3.65,8.47,2.81)]:ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=10,color='#56768C'))
 ax.text(5,.55,'Many gene changes can be combinations of a few shared expression directions.',ha='center',fontsize=8.2);save(fig,'Figure_1')
 d=read('dimension_parallel_analysis.tsv');dec=read('dimension_decision.tsv');ks=dec[dec.scaling=='native'].set_index('effect_type').selected_K.to_dict()
 fig,axes=plt.subplots(1,3,figsize=(6.69,2.8),layout='constrained')
 for ax,e,n,col in zip(axes,E,N,C):
  x=d[(d.effect_type==e)&(d.scaling=='native')];ax.plot(x['rank'],x.excess_ratio,'o-',ms=3,c=col);ax.axhline(1,c='#666666',ls='--',lw=.8);ax.axvspan(.7,ks[e]+.25,color=col,alpha=.07);ax.set(title=f'{n}\nSelected K = {ks[e]}',xlabel='Direction',xticks=[1,5,10,15],xlim=(.5,15.5));ax.tick_params(labelsize=7)
 axes[0].set_ylabel('Eigenvalue / permutation threshold');save(fig,'Figure_2')
 m=read('selected_metrics.tsv');dom=read('selected_domain_transfer.tsv');nul=read('selected_null.tsv')
 diseases=['Malaria','SSc-associated ILD','Burn injury','Acute myocardial infarction','Severe asthma','Trisomy 21','Bipolar disorder','Chronic granulomatous disease','Diffuse large B-cell lymphoma','Ulcerative colitis','Colorectal cancer','Autoimmune thyroiditis']
 fig,axes=plt.subplots(2,1,figsize=(6.69,6.2),layout='constrained',gridspec_kw={'height_ratios':[.8,1.9]})
 for i,(e,n,col) in enumerate(zip(E,N,C)):
  vv=nul[nul.effect_type==e].median_R2;axes[0].scatter(np.full(len(vv),i-.10),vv,s=5,alpha=.13,c='gray');axes[0].scatter(i+.10,m[m.effect_type==e].R2.median(),s=35,c=col,zorder=3)
 axes[0].set(xticks=range(3),xticklabels=N,ylabel='Median test R²',ylim=(-.04,.52),title='a   Learned and gene-permuted spaces');axes[0].text(.02,.52,'Coloured: learned space\nGrey: 199 permuted spaces',transform=axes[0].transAxes,va='top',fontsize=7)
 order=dom[dom.effect_type==E[0]].contrast_id.tolist()
 for j,(e,n,col) in enumerate(zip(E,N,C)):
  x=dom[dom.effect_type==e].set_index('contrast_id').loc[order];axes[1].scatter(x.R2,np.arange(12)+(j-1)*.17,c=col,label=n,s=17)
 axes[1].set(yticks=range(12),yticklabels=diseases,xlabel='Category-excluded reconstruction R²',title='b   Learning without the test disease category',xlim=(0,.85));axes[1].invert_yaxis();axes[1].tick_params(axis='y',labelsize=7.5);axes[1].legend(loc='lower right',fontsize=7,frameon=False);axes[1].grid(axis='x',alpha=.15);save(fig,'Figure_3')
 en=read('selected_enrichment.tsv');paths=['HALLMARK_HEME_METABOLISM','HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_IL2_STAT5_SIGNALING','HALLMARK_OXIDATIVE_PHOSPHORYLATION'];labels=['Heme metabolism','Allograft rejection','TNF signalling via NF-κB','IL2–STAT5 signalling','Oxidative phosphorylation'];cols=[(e,k) for e in E for k in range(1,ks[e]+1)];z=np.full((5,len(cols)),np.nan)
 for j,(e,k) in enumerate(cols):
  for i,p in enumerate(paths):
   v=en[(en.effect_type==e)&(en.factor==k)&(en.pathway==p)].q_selected.min()
   if v<.05:z[i,j]=min(30,-np.log10(v))
 fig,ax=plt.subplots(figsize=(6.69,2.7),layout='constrained');cmap=plt.get_cmap('YlOrRd').copy();cmap.set_bad('#EEEEEE');im=ax.imshow(z,aspect='auto',cmap=cmap,vmin=0,vmax=30);ax.set(yticks=range(5),yticklabels=labels,xticks=range(len(cols)),xticklabels=[f'{dict(zip(E,["g","FC","sFC"]))[e]}{k}' for e,k in cols]);ax.tick_params(labelsize=7)
 for b in [4.5,8.5]:ax.axvline(b,c='white',lw=2)
 cb=fig.colorbar(im,ax=ax,shrink=.88,pad=.025);cb.set_label('−log₁₀ adjusted q\n(capped at 30)',fontsize=7);cb.ax.tick_params(labelsize=7);ax.set_xlabel('Direction within each effect representation');save(fig,'Figure_4')
 sens=read('selected_sensitivity.tsv');boot=read('selected_bootstrap.tsv');fig=plt.figure(figsize=(6.69,5.9),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1]);a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,:])
 a.boxplot([boot[boot.effect_type==e].overlap for e in E],tick_labels=['g','FC','sFC'],showfliers=True,flierprops={'marker':'.','markersize':3});a.set(ylabel='Subspace overlap',ylim=(0,1.05),title='a   Study-group resampling')
 for e,n,col in zip(E,N,C):
  yy=[m[m.effect_type==e].R2.median()]+[sens[(sens.effect_type==e)&(sens.analysis==s)].R2.median() for s in ['equal_block_weight','exclude_sickle_MPN']];b.plot(range(3),yy,'o-',c=col,ms=3,label=n)
 b.set(xticks=range(3),xticklabels=['Original','Equal group\nweight','Exclude selected\nblood diseases'],ylabel='Median test R²',ylim=(.28,.48),title='b   Same full-panel target');b.tick_params(axis='x',labelsize=6.5)
 for i,(e,n,col) in enumerate(zip(E,N,C)):
  y=[sens[(sens.effect_type==e)&(sens.analysis==s)].R2.median() for s in ['reference_on_nonheme','refit_without_heme']];c.plot([i-.13,i+.13],y,c=col,lw=1);c.scatter(i-.13,y[0],facecolors='white',edgecolors=col,s=45,zorder=3);c.scatter(i+.13,y[1],c=col,s=45,zorder=3)
 c.set(xticks=range(3),xticklabels=N,ylabel='Median test R²',ylim=(.28,.49),title='c   Matched non-heme target');c.text(.02,.97,'Open: original space projected onto non-heme genes\nFilled: space refitted without heme genes',transform=c.transAxes,va='top',fontsize=7);save(fig,'Figure_5')
if __name__=='__main__':main()
