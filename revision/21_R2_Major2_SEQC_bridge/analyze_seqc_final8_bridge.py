import argparse, csv, gzip, pathlib, re, math, statistics
from collections import defaultdict
HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_REFERENCE = HERE.parents[1] / 'parameters' / 'TCGA_reference.tsv'
DEFAULT_OUTPUT = HERE.parent / 'analyses' / 'R2_Major2_SEQC_bridge'
parser = argparse.ArgumentParser(description='Reproduce the exploratory Final8 technical-platform bridge using public SEQC reference materials.')
parser.add_argument('--input-dir', required=True, type=pathlib.Path, help='Directory containing the four downloaded SEQC/GPL source files; see docs/DATA_SOURCES.md.')
parser.add_argument('--reference', type=pathlib.Path, default=DEFAULT_REFERENCE, help='Final8 reference TSV (default: repository parameters/TCGA_reference.tsv).')
parser.add_argument('--output-dir', type=pathlib.Path, default=DEFAULT_OUTPUT, help='Directory for aggregate bridge tables.')
args = parser.parse_args()
INPUT = args.input_dir.expanduser().resolve()
REF = args.reference.expanduser().resolve()
OUT = args.output_dir.expanduser().resolve()
OUT.mkdir(parents=True, exist_ok=True)
required = ['GPL10558_SOFT.txt', 'GSE56457_SEQC_qPCR_GEOSub.txt.gz', 'GSE47774_SEQC_ILM_COH.txt.gz', 'GSE56457_SEQC_ILM_GEOSub.txt.gz']
missing = [name for name in required if not (INPUT / name).is_file()]
if missing:
    parser.error('Missing input file(s) in --input-dir: ' + ', '.join(missing))
if not REF.is_file():
    parser.error('Reference TSV not found: ' + str(REF))
with REF.open('rt', encoding='utf8', newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
genes=[r['gene'] for r in rows]
geneids={r['gene']:r['gene_id'].split('.')[0] for r in rows}

gpl=INPUT/'GPL10558_SOFT.txt'; probe_symbols={}; refseq_symbols={}
with gpl.open('rt',encoding='utf8',errors='replace') as f:
 in_table=False; hdr=None
 for line in f:
  line=line.rstrip('\r\n')
  if line=='!platform_table_begin': in_table=True; hdr=None; continue
  if line=='!platform_table_end': break
  if in_table:
   parts=line.split('\t')
   if hdr is None: hdr=parts; continue
   if len(parts)>=len(hdr):
    d=dict(zip(hdr,parts)); sym=d.get('Symbol','').strip()
    if sym in genes:
     probe_symbols[parts[0]]=sym
     for rid in re.split(r'[;,| ]+',d.get('RefSeq_ID','').strip()):
      if rid: refseq_symbols[rid.split('.')[0]]=sym

qfile=INPUT/'GSE56457_SEQC_qPCR_GEOSub.txt.gz'
with gzip.open(qfile,'rt',encoding='utf8',errors='replace',newline='') as f:
 qrows=list(csv.DictReader(f,delimiter='\t'))

rfile=INPUT/'GSE47774_SEQC_ILM_COH.txt.gz'
with gzip.open(rfile,'rt',encoding='utf8',errors='replace',newline='') as f:
 reader=csv.reader(f,delimiter='\t'); header=next(reader); groups=defaultdict(list)
 for col in header[1:]:
  m=re.search(r'_([ABCD])_(\d+)_L\d+_',col)
  if m: groups[f'{m.group(1)}_{m.group(2)}'].append(header.index(col))
 libcounts=defaultdict(lambda:defaultdict(float)); totals=defaultdict(float); transcript_ids=set()
 for row in reader:
  tid=row[0].split('.')[0]; transcript_ids.add(tid); gene=refseq_symbols.get(tid)
  try: vals=[float(x or 0) for x in row[1:]]
  except ValueError: continue
  for lib,idxs in groups.items():
   val=sum(vals[i-1] for i in idxs)
   totals[lib]+=val
   if gene in genes: libcounts[lib][gene]+=val

assays=defaultdict(list)
for q in qrows:
 gene=q.get('official_symbol','').strip()
 if gene in genes:
  assays[gene].append((q.get('targeted',''),q))

rna={g:{} for g in genes}
for g in genes:
 for lib in libcounts:
  raw=libcounts[lib].get(g,0)
  cpm=(raw+0.5)/(totals[lib]+1)*1e6
  grp=lib.split('_')[0]; rna[g].setdefault(grp,[]).append(math.log2(cpm))
for g in genes:
 for grp,vs in list(rna[g].items()): rna[g][grp]=statistics.mean(vs)

qpcr={g:{} for g in genes}
for g in genes:
 for grp,col in zip('ABCD',['SEQC_RTPCR_A','SEQC_RTPCR_B','SEQC_RTPCR_C','SEQC_RTPCR_D']):
  vv=[]
  for tid,q in assays[g]:
   try: vv.append(float(q[col]))
   except (ValueError,TypeError): pass
  if vv: qpcr[g][grp]=statistics.mean(vv)
qPCR_signal={g:{grp:-v for grp,v in qpcr[g].items()} for g in genes}

mfile=INPUT/'GSE56457_SEQC_ILM_GEOSub.txt.gz'
arr={g:{} for g in genes}; detected=defaultdict(int)
with gzip.open(mfile,'rt',encoding='utf8',errors='replace',newline='') as f:
 rr=csv.reader(f,delimiter='\t'); h=next(rr)
 sample_cols={}
 for j,col in enumerate(h[1:],1):
  m=re.match(r'SEQC_ILM_COH_([ABCD])_(\d+)_MA$',col)
  if m: sample_cols[(m.group(1),m.group(2))]=(j,j+1)
 vals=defaultdict(lambda:defaultdict(list))
 for row in rr:
  sym=probe_symbols.get(row[0]);
  if not sym: continue
  for (grp,rep),(j,pj) in sample_cols.items():
   try: x=float(row[j]); p=float(row[pj])
   except (ValueError,IndexError): continue
   vals[(sym,grp)][rep].append((x,p))
 for (sym,grp),reps in vals.items():
  runs=[]
  for rep,pairs in reps.items():

   keep=[x for x,p in pairs if p<0.05]
   if keep: runs.append(math.log2(statistics.mean(keep)))
  if runs: arr[sym][grp]=statistics.mean(runs)

def corr(x,y):
 n=len(x)
 if n<3: return float('nan')
 mx=sum(x)/n; my=sum(y)/n
 den=math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
 return sum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else float('nan')
def ranks(v):
 order=sorted(range(len(v)),key=lambda i:v[i]); out=[0]*len(v); k=0
 while k<len(v):
  l=k
  while l+1<len(v) and v[order[l+1]]==v[order[k]]: l+=1
  rank=(k+l+2)/2
  for j in range(k,l+1):out[order[j]]=rank
  k=l+1
 return out

with (OUT/'SEQC_Final8_reference_expression_long.tsv').open('w',newline='',encoding='utf8') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['gene','reference_material','RNAseq_log2CPM_RefSeq_mapped_gene' ,'microarray_log2_mean_detected_probes','qPCR_raw_Cq','qPCR_inverse_Cq_signal','n_RNAseq_RefSeq_accessions_mapped' ,'n_microarray_probes','n_qPCR_matched_assays'])
 for g in genes:
  nprobe=sum(1 for p,s in probe_symbols.items() if s==g)
  for grp in 'ABCD':
   w.writerow([g,grp,rna[g].get(grp,''),arr[g].get(grp,''),qpcr[g].get(grp,''),qPCR_signal[g].get(grp,''),len(assays[g]),nprobe,len(assays[g])])
with (OUT/'SEQC_Final8_pairwise_profile_correlations.tsv').open('w',newline='',encoding='utf8') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['gene','platform_pair','n_reference_materials','Pearson_r','Spearman_rho','RNAseq_RefSeq_accessions' ,'microarray_probes','qPCR_assays'])
 for g in genes:
  for p1,p2,d1,d2 in [('RNAseq_vs_microarray','RNAseq','RNAseq','microarray'),('RNAseq_vs_qPCR','RNAseq','RNAseq','qPCR'),('microarray_vs_qPCR','microarray','microarray','qPCR')]:
   common=[c for c in 'ABCD' if c in ({'RNAseq':rna,'microarray':arr,'qPCR':qpcr}[d1][g]) and c in ({'RNAseq':rna,'microarray':arr,'qPCR':qPCR_signal}[d2][g])]
   x=[{'RNAseq':rna,'microarray':arr,'qPCR':qpcr}[d1][g][c] for c in common]; y=[{'RNAseq':rna,'microarray':arr,'qPCR':qPCR_signal}[d2][g][c] for c in common]
   w.writerow([g,p1,len(common),corr(x,y) if len(common)>=3 else '',corr(ranks(x),ranks(y)) if len(common)>=3 else '',len(assays[g]),sum(1 for p,s in probe_symbols.items() if s==g),len(assays[g])])
with (OUT/'SEQC_Final8_mapping_QC.tsv').open('w',newline='',encoding='utf8') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['gene','Ensembl_gene_ID_reference','qPCR_assays_found','RNAseq_RefSeq_accessions_mapped','GPL10558_probes','microarray_detectable_any_reference'])
 for g in genes:
  w.writerow([g,geneids[g],len(assays[g]),len(assays[g]),sum(1 for p,s in probe_symbols.items() if s==g),int(bool(arr[g]))])
print('genes',genes)
print('RNA libraries',len(libcounts),'qPCR mapped genes',sum(bool(assays[g]) for g in genes),'array-mapped genes',sum(bool([p for p,s in probe_symbols.items() if s==g]) for g in genes))
print('wrote bridge tables')

