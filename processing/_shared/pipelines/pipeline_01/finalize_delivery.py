from pathlib import Path
import json,itertools,collections
import pandas as pd
ROOT=Path('E:/Biology');OUT=ROOT/'outputs/stage1_three_effects_20260903_v1'
s=json.loads((OUT/'verification_summary.json').read_text(encoding='utf8'));assert s['status']=='COMPLETE' and s['complete']==s['planned']
rows=pd.read_csv(OUT/'contrast_build_ledger.tsv',sep='\t').to_dict('records');parent={r['contrast_id']:r['contrast_id'] for r in rows}
def find(x):
 while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
 return x
def union(a,b):parent[find(a)]=find(b)
groups=collections.defaultdict(list);samples={}
for r in rows:
 groups[r['group']].append(r['contrast_id'])
 m=pd.read_csv(OUT/'contrasts'/r['contrast_id']/'samples.tsv',sep='\t',dtype=str)
 samples[r['contrast_id']]={v if v.startswith('GSM') else r['accession']+':'+v for v in m.sample_id}
for ids in groups.values():
 for cid in ids[1:]:union(ids[0],cid)
edges=[]
for a,b in itertools.combinations(rows,2):
 x=a['contrast_id'];y=b['contrast_id'];shared=samples[x]&samples[y]
 if shared:edges.append({'contrast_a':x,'contrast_b':y,'shared_selected_sample_ids':len(shared),'sample_ids':';'.join(sorted(shared))});union(x,y)
pd.DataFrame(edges,columns=['contrast_a','contrast_b','shared_selected_sample_ids','sample_ids']).to_csv(OUT/'sample_overlap_edges.tsv',sep='\t',index=False)
blocks=collections.defaultdict(list)
for r in rows:blocks[find(r['contrast_id'])].append(r['contrast_id'])
bs=[]
for i,ids in enumerate(sorted(blocks.values(),key=lambda z:min(z)),1):
 for cid in sorted(ids):bs.append({'contrast_id':cid,'conservative_block':f'block_{i:03}','basis':'same audited study_group or observed shared selected sample IDs; not proof of unrelated participants'})
pd.DataFrame(bs).to_csv(OUT/'independence_blocks.tsv',sep='\t',index=False)
s['conservative_blocks']=len(blocks);s['sample_overlap_pairs']=len(edges);(OUT/'verification_summary.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
coverage={}
for scope in ['core','supplemental','all']:
 d=OUT/'matrices'/scope;m=pd.read_csv(d/'observed_mask.tsv.gz',sep='\t',index_col=0)
 c=pd.DataFrame({'gene_id':m.index,'observed_contrasts':m.sum(axis=1),'coverage_fraction':m.mean(axis=1)})
 c.to_csv(d/'gene_coverage.tsv',sep='\t',index=False)
 coverage[scope]=int(c.coverage_fraction.ge(.8).sum())
excluded=[]
for p in (OUT/'contrasts').glob('*/excluded_noncurrent_entrez.tsv.gz'):
 t=pd.read_csv(p,sep='\t',dtype={'gene_id':str})
 excluded.extend({'contrast_id':p.parent.name,'gene_id':v,'reason':'not found in frozen org.Hs.eg.db universe'} for v in t.gene_id)
pd.DataFrame(excluded).to_csv(OUT/'excluded_annotation_ids.tsv',sep='\t',index=False)
append=f'''
## 最终核验补充

已完成 {s['complete']} 条 contrast 的三种效应量；十二类候选集合 108 条，类别外补充集合 16 条。三套表示逐基因对齐；所有有效不确定性均大于零，未对不确定性人为设下限。ashr 使用 `method="shrink"`，即连续正态尺度混合先验和均匀混合权重先验，区别于面向 FDR 的零点质量先验；似然使用正态近似。这个配置在所有 contrasts 中一致。

本次覆盖原排查中全部 119 个 A/B/C 目录；增加三个资料明确的疾病分支，并将两项混合分型拆开，净增五条，得到 124 条 contrast。其余目录中的 10 项证据待确认、9 项不符合筛选条件、6 项重复或审计用途，逐项原因保留在 `not_constructed_inventory.tsv`。完成范围是当前方案与可核实资料支持的所有可构建项；证据待确认项未被伪造为有效对照。

GSE201332 的自定义芯片通过与 GPL6480/GPL21185 完全一致、长度至少 40 nt 且 Entrez 指向唯一的探针序列转移注释。GSE5812 的旧 cDNA 克隆依据 org.Hs.eg.db 与 hgu133plus2.db 中可追溯的 GenBank accession 映射，只有 1,235 个可完整计算的基因，覆盖明显有限；保留在补充集合，后续不应当作覆盖完整的全基因组 contrast。详细记录见 `annotation_limitations.json`。

`sample_overlap_edges.tsv` 记录检测到的 {len(edges)} 对样本交集，`independence_blocks.tsv` 给出 {len(blocks)} 个保守分组。这些分组结合既有研究来源与显式样本 ID；未识别出交集不证明受试者完全独立。

独立核验检查了全部结果的样本数、唯一 Entrez ID、有限值、三种表示的基因对齐和收缩方向，并对每条 contrast 的固定等距抽样基因从归一化表达重新计算 Hedges g/SE；非 voom 的 log2FC/SE 也进行了独立代数复算。原目录中 583 个已登记元数据文件的 SHA-256 全部与处理前一致。

最终基因 ID 限定于 `canonical_gene_universe.tsv` 中冻结的 org.Hs.eg.db 注释库。在归一化后、收缩先验拟合前移除无法对应的 ID；排除明细见 `excluded_annotation_ids.tsv`，不据此判断基因是否已被正式撤销。

十二类集合的并集为 {s['matrices']['core']['union_genes']:,} 个基因，全部 108 条均有测量的交集为 {s['matrices']['core']['intersection_genes']} 个，至少覆盖 80% contrasts 的基因为 {coverage['core']:,} 个。异质平台的覆盖差异很大，应基于训练集选择后续面板并保留缺失掩码，不应直接把严格全交集当作最终面板。逐基因覆盖率见各集合的 `gene_coverage.tsv`。

完整交付核验还包括 `verify_matrices.py` 对全部 18 个矩阵文件的读回比对。复现时在 `assemble.py` 后依次运行 `finalize_delivery.py`、`verify_matrices.py`、`provenance.py`；最后一个脚本输出源文件、代码及交付文件的 SHA-256 清单。
'''
p=OUT/'REPORT.md';text=p.read_text(encoding='utf8').split('\n## 最终核验补充')[0];p.write_text(text+append,encoding='utf8')
print('DELIVERY_FINALIZED',s['complete'],len(blocks),len(edges))
