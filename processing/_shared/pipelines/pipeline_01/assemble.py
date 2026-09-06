"""Audit all committed results and assemble matched masked E/U matrices."""
from pathlib import Path
import json,hashlib,sys,collections
import numpy as np,pandas as pd
from scipy.special import gammaln
sys.stdout.reconfigure(encoding='utf8')
ROOT=Path('E:/Biology');OUT=ROOT/'outputs/stage1_three_effects_20260903_v1';AUD=ROOT/'outputs/dataset_contrast_audit_20260903_v1'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def tab(p):return pd.read_csv(p,sep='\t',dtype={'gene_id':str})
plan=json.loads((OUT/'planned_contrasts.json').read_text(encoding='utf8'));inventory=json.loads((AUD/'audit_data.json').read_text(encoding='utf8'))['inventory']
views={'hedges_g':('hedges_g','se_hedges_g'),'log2fc':('log2fc','se_log2fc'),'shrunken_log2fc':('shrunken_log2fc','posterior_sd')}
canonical=set(pd.read_csv(OUT/'canonical_gene_universe.tsv',sep='\t',dtype=str).gene_id)
GZIP={'method':'gzip','compresslevel':1,'mtime':0}
rows=[];tables={};checks=[]
for r in plan:
 d=OUT/'contrasts'/r['contrast_id'];p=d/'qc.json'
 if not p.exists():rows.append(dict(r,build_status='PENDING',genes=0));continue
 q=json.loads(p.read_text());info=json.loads((d/'input.json').read_text(encoding='utf8'));t=tab(d/'three_effects.tsv.gz').set_index('gene_id');samples=tab(d/'samples.tsv')
 assert sha(d/'input_expression.tsv.gz')==info['input_sha256'],d
 assert t.index.is_unique and t.index.str.fullmatch(r'\d+').all(),d
 assert set(t.index)<=canonical,d
 assert np.isfinite(t.to_numpy()).all() and (t[['se_hedges_g','se_log2fc']]>0).all().all() and (t.posterior_sd>=0).all(),d
 assert ((t.shrunken_log2fc*t.log2fc)>=-1e-10).all() and (t.shrunken_log2fc.abs()<=t.log2fc.abs()+1e-8).all(),d
 assert len(samples)==q['n_case']+q['n_control'] and not samples.sample_id.duplicated().any(),d
 assert samples.group.eq('case').sum()==q['n_case'] and samples.group.eq('control').sum()==q['n_control'],d
 assert q['shrinkage']=='ashr_continuous_normal_mixture_normal_likelihood',d
 for name,(e,u) in views.items():
  z=tab(d/(name+'_E_U.tsv.gz')).set_index('gene_id');assert z.index.equals(t.index),d
  assert np.allclose(z.effect,t[e]) and np.allclose(z.uncertainty,t[u]),d
 # Independent algebraic checks on fixed, outcome-independent sampled gene rows.
 ix=np.linspace(0,len(t)-1,min(31,len(t)),dtype=int);want=set(t.index[ix]);x=tab(d/'normalized_expression.tsv.gz').set_index('gene_id').loc[t.index[ix]]
 assert list(x.columns)==samples.sample_id.astype(str).tolist(),d
 n1=q['n_case'];n0=q['n_control'];a=x.iloc[:,:n1].to_numpy();b=x.iloc[:,n1:].to_numpy();df=n1+n0-2
 va=a.var(1,ddof=1);vb=b.var(1,ddof=1);delta=(a.mean(1)-b.mean(1))/np.sqrt(((n1-1)*va+(n0-1)*vb)/df)
 J=np.exp(gammaln(df/2)-.5*np.log(df/2)-gammaln((df-1)/2));g=J*delta;se=J*np.sqrt((n1+n0)/(n1*n0)+delta**2/(2*df))
 assert np.allclose(g,t.iloc[ix].hedges_g,rtol=2e-7,atol=2e-7) and np.allclose(se,t.iloc[ix].se_hedges_g,rtol=2e-7,atol=2e-7),d
 if q['log2fc_model']=='OLS_on_log2_expression':
  assert np.allclose(a.mean(1)-b.mean(1),t.iloc[ix].log2fc,rtol=2e-7,atol=2e-7),d
  seols=np.sqrt(((n1-1)*va+(n0-1)*vb)/df*(1/n1+1/n0));assert np.allclose(seols,t.iloc[ix].se_log2fc,rtol=2e-7,atol=2e-7),d
 checks.append({'contrast_id':r['contrast_id'],'status':'PASS','algebraic_genes_checked':len(ix),'genes':len(t)})
 rows.append(dict(r,build_status='COMPLETE',genes=len(t),n_case=q['n_case'],n_control=q['n_control'],scope=info['scope'],transform=q['expression_transform'],model=q['log2fc_model'],posterior_sd_zero=q['posterior_sd_zero'],result_path=str(d/'three_effects.tsv.gz')));tables[r['contrast_id']]=t
 print('VERIFIED',r['contrast_id'],len(t),flush=True)
if '--validate-only' in sys.argv:
 pd.DataFrame(checks).to_csv(OUT/'validation_partial.tsv',sep='\t',index=False);print('PASS',len(checks),'completed contrasts');sys.exit(0)
ledger=pd.DataFrame(rows);ledger.to_csv(OUT/'contrast_build_ledger.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(OUT/'validation_checks.tsv',sep='\t',index=False)
sizes={}
for scope in ['core','supplemental','all']:
 ids=[r['contrast_id'] for r in rows if r['build_status']=='COMPLETE' and (scope=='all' or r['scope']==scope)]
 if not ids:continue
 md=OUT/'matrices'/scope;md.mkdir(parents=True,exist_ok=True)
 union=sorted(set.union(*(set(tables[c].index) for c in ids)),key=int);common=set.intersection(*(set(tables[c].index) for c in ids))
 pd.DataFrame({'gene_id':union}).to_csv(md/'gene_index.tsv',sep='\t',index=False)
 ledger.set_index('contrast_id').loc[ids].reset_index().to_csv(md/'contrast_index.tsv',sep='\t',index=False)
 for name,(e,u) in views.items():
  E=pd.DataFrame({c:tables[c][e] for c in ids}).reindex(union);U=pd.DataFrame({c:tables[c][u] for c in ids}).reindex(union);E.index.name=U.index.name='gene_id'
  E.to_csv(md/(name+'_E.tsv.gz'),sep='\t',na_rep='NA',compression=GZIP);U.to_csv(md/(name+'_U.tsv.gz'),sep='\t',na_rep='NA',compression=GZIP)
 mask=pd.DataFrame({c:pd.Series(True,index=tables[c].index) for c in ids}).reindex(union).notna().astype(int);mask.index.name='gene_id';mask.to_csv(md/'observed_mask.tsv.gz',sep='\t',compression=GZIP)
 pd.DataFrame({'gene_id':sorted(common,key=int)}).to_csv(md/'complete_case_gene_intersection.tsv',sep='\t',index=False)
 sizes[scope]={'contrasts':len(ids),'union_genes':len(union),'intersection_genes':len(common),'observed_entries':int(mask.to_numpy().sum())}
domains=['感染性疾病','系统性自身免疫与风湿性疾病','移植与急性组织损伤','心代谢与血管疾病','慢性呼吸系统疾病','遗传与先天性疾病','神经系统疾病','原发性免疫缺陷','血液系统肿瘤','炎症性肠病','实体肿瘤','内分泌自身免疫疾病']
counts=[]
for domain in domains:
 rr=[r for r in rows if r.get('domain')==domain and r.get('scope')=='core' and r['build_status']=='COMPLETE'];counts.append({'domain':domain,'contrasts':len(rr),'accessions':len(set(r['accession'] for r in rr))})
pd.DataFrame(counts).to_csv(OUT/'twelve_category_counts.tsv',sep='\t',index=False)
ex=[r for r in inventory if r['status'][0] not in 'ABC'];pd.DataFrame(ex).to_csv(OUT/'not_constructed_inventory.tsv',sep='\t',index=False)
pending=[r['contrast_id'] for r in rows if r['build_status']!='COMPLETE'];summary={'status':'COMPLETE' if not pending else 'IN_PROGRESS','planned':len(plan),'complete':len(checks),'pending':pending,'matrices':sizes,'category_counts':counts}
(OUT/'verification_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
report=['# 第一阶段：三种效应量构建','',f"状态：{summary['status']}；已完成 {len(checks)} / {len(plan)} 条预定义可构建 contrast。",'', '| 十二大类 | 已完成 contrasts | GEO accession 数 |','|---|---:|---:|']
report += [f"| {r['domain']} | {r['contrasts']} | {r['accessions']} |" for r in counts]
report += ['', '## 文件与统计口径','', '`contrast_build_ledger.tsv` 是逐条构建总表；`contrasts/<contrast_id>/` 包含实际样本、输入矩阵、基因映射、归一化表达、三套 E/U、ash 拟合及 QC。', '', '`matrices/core/` 是十二类候选集合，`matrices/supplemental/` 是类别尚未归入十二类的补充集合，`matrices/all/` 合并两者。每个集合的三种表示具有完全相同的行列顺序和缺失掩码。不同研究未覆盖的基因用 NA 表示，未填零；交集另行提供，不预先锁定下游训练基因面板。', '', '## 固定方法','', '原始计数先按可靠 Entrez 映射汇总，再用 edgeR filterByExpr、TMM 和 limma voom。芯片或已标准化表达按记录的量纲进行 log2 转换，线性芯片强度采用下限 1 的 log2 与分位数归一化。FPKM/TPM 等连续丰度采用 log2(x+1)，不冒充原始计数。基因汇总发生在效应估计之前，未按显著性或标准误挑探针。', '', 'Hedges g 使用归一化 log2 表达的病例减健康对照标准化均值差、精确 Gamma 比值小样本校正 J，以及常用近似抽样 SE。log2FC 使用两组设计的线性模型系数及未调节抽样 SE，计数数据使用 voom 权重。收缩 log2FC 使用 ashr 正态混合先验和正态似然近似，保存后验均值与后验 SD；后验 SD 不称为抽样 SE。', '', '三种表示均为当前病例/对照定义下的未协变量调整效应。年龄、性别、药物和批次等潜在混杂仍需在后续敏感性分析中处理；这批文件不表示因果疾病效应。小样本下的 Hedges 方差与 ashr 正态似然均有近似性质。', '', '## 已修正的样本定义','', 'GSE181228 仅纳入有计数的 132 个基线病例与 43 个健康对照；GSE94648 去除 4 个 bis 技术重复后拆分 CD 48/20 和 UC 25/20；GSE44314 将经典 1A 型 5/6 与暴发型 5/6 分开，后者保留在补充集合。GSE177044 的 UC 为 495/320；GSE18781 按真实 set 与疾病分组。GSE100156 与 GSE100150 重复样本、历史重复结果和审计目录不作为新增 contrasts。', '', '另构建本地证据已明确的 GSE42825 TB、GSE68004 不完全型川崎病、GSE69683 重度哮喘分支，保留其共享健康对照关系。', '', '## 使用限制与复现','', 'contrast 数不等于独立队列数。应结合 study_group、样本交集与既有来源证据进行分组划分；不得将共享对照、重复参与者或同一研究分支随机拆到训练和测试两侧。当前全基因经验先验属于第一阶段数据快照，严格盲测时须在训练折内重拟合，不能把已收缩的全基因矩阵当作未经处理的盲测输入。', '', '本次完成标记仅表示三种效应量的技术构建及本地验证完成，不覆盖后续低维模型、交叉验证、因果解释或中央 EU_locked 流程。', '', '原始 dataset 及旧结果保留。复现代码位于 `E:/Biology/modeling/stage1_three_effects_20260903_v1/`。依次运行 `cel.R`、`build.py`、`run_effects.py`、`assemble.py`；R 的固定库路径和实际 sessionInfo 已保存。', '', '[方案原文](https://app.notion.com/p/3d0d446bf21a81209bfae0642b792310)；[ashr 作者文档](https://stephens999.r-universe.dev/ashr/doc/manual.html)。']
if pending:report+=['','未完成：'+', '.join(pending)]
(OUT/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf8')
print(json.dumps(summary,ensure_ascii=False,indent=2))

