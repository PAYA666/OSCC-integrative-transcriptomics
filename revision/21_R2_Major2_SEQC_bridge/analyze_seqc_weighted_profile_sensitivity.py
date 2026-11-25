pass
import csv, math, statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
LONG = HERE / "SEQC_Final8_reference_expression_long.tsv"
REF = HERE.parent / "12_single_sample" / "TCGA_reference.tsv"
OUT = HERE / "SEQC_Final8_weighted_profile_sensitivity.tsv"
with REF.open(encoding="utf-8-sig", newline="") as f:
    beta={r['gene']:float(r['beta']) for r in csv.DictReader(f,delimiter='\t')}
with LONG.open(encoding="utf-8-sig", newline="") as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
cols={'RNAseq':'RNAseq_log2CPM_RefSeq_mapped_gene','qPCR':'qPCR_inverse_Cq_signal','microarray':'microarray_log2_mean_detected_probes'}
raw={a:{} for a in cols}
for r in rows:
    g=r['gene']; sample=r['reference_material']
    for assay,col in cols.items():
        try: raw[assay].setdefault(g,{})[sample]=float(r[col])
        except (ValueError,TypeError): pass

def corr(x,y,rank=False):
    if rank:
        def ranks(a):
            order=sorted(range(len(a)),key=a.__getitem__); out=[0.]*len(a); j=0
            while j<len(a):
                k=j
                while k+1<len(a) and a[order[k+1]]==a[order[j]]:k+=1
                for h in range(j,k+1):out[order[h]]=(j+k+2)/2
                j=k+1
            return out
        x,y=ranks(x),ranks(y)
    mx,my=statistics.mean(x),statistics.mean(y)
    den=math.sqrt(sum((v-mx)**2 for v in x)*sum((v-my)**2 for v in y))
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else float('nan')

out=[]
for a,b in [('RNAseq','qPCR'),('RNAseq','microarray'),('qPCR','microarray')]:
    common=sorted(s for s in 'ABCD' if all(s in raw[k].get(g,{}) for k in (a,b) for g in beta))
    pair_common={g:sorted(set(raw[a].get(g,{})) & set(raw[b].get(g,{}))) for g in beta}
    coverage_genes=[g for g in beta if len(pair_common[g])>=3]

    if len(coverage_genes)==len(beta) and len(common)>=3:
        profiles={}
        for assay in (a,b):
            prof=[]
            for g in beta:
                vals=[raw[assay][g][s] for s in common]
                sd=statistics.stdev(vals)
                zg={s:(raw[assay][g][s]-statistics.mean(vals))/sd for s in common}
                prof.append(zg)
            profiles[assay]=[sum(beta[g]*profiles_g[i][s] for i,g in enumerate(beta) for profiles_g in [prof]) for s in common]
        pearson=corr(profiles[a],profiles[b]); spearman=corr(profiles[a],profiles[b],True)
        status='All eight genes; assay-wise z standardization across reference mixtures'
    else:
        pearson=spearman=''
        status=f'Not estimable as a full Final8 composite; only {len(coverage_genes)}/8 genes have >=3 matched materials'
    refs=common if common else ['']
    for s in refs:
        out.append({'assay_pair':f'{a}_vs_{b}','reference_material':s,'n_common_full_profile_references':len(common),
                    'genes_with_ge3_pairwise_references':len(coverage_genes),
                    'score_a':profiles[a][common.index(s)] if status.startswith('All eight') else '',
                    'score_b':profiles[b][common.index(s)] if status.startswith('All eight') else '',
                    'pearson_profile_r':pearson,'spearman_profile_rho':spearman,'interpretation':status})
with OUT.open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
print('Wrote',OUT)
for pair in sorted({r['assay_pair'] for r in out}):
    r=next(x for x in out if x['assay_pair']==pair)
    print(pair,r['interpretation'],'r=',r['pearson_profile_r'],'rho=',r['spearman_profile_rho'])
