import argparse, csv, gzip, math, random, statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_MATRIX = HERE.parent / "data" / "GSE84846_series_matrix.txt.gz"
DEFAULT_REFERENCE = HERE.parents[1] / "parameters" / "TCGA_reference.tsv"
BOOTSTRAPS = 5000
SEED = 20261125

parser = argparse.ArgumentParser(description='GSE84846 binary-status analysis')
parser.add_argument("--matrix", type=Path, required=True,
                    help="GSE84846 processed GEO series matrix (.txt.gz)")
parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE,
                    help="Frozen Final8 reference table with gene and beta columns")
parser.add_argument("--output-dir", type=Path, default=HERE,
                    help="Directory for per-sample scores and aggregate report")
args = parser.parse_args()
MATRIX = args.matrix.resolve()
REFERENCE = args.reference.resolve()
OUT_DIR = args.output_dir.resolve()
if not MATRIX.is_file():
    parser.error(f"GSE84846 input matrix not found: {MATRIX}")
if not REFERENCE.is_file():
    parser.error(f"Final8 reference table not found: {REFERENCE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_SCORES = OUT_DIR / "GSE84846_Final8_binary_audit_recomputed.tsv"
OUT_SUMMARY = OUT_DIR / "GSE84846_frozen_score_binary_analysis.txt"
PROBES = {
    "CTSG": ["A_23_P140384"],
    "P4HA1": ["A_24_P406693"],
    "NKAIN1": ["A_23_P51376"],
    "LRRC38": ["A_23_P103877", "A_24_P25794"],
    "SPSB4": ["A_23_P258381"],
    "SPINK6": ["A_23_P92909"],
    "TNFRSF4": ["A_32_P26092"],
    "DEFB1": ["A_23_P71480"],
}

def read_matrix(path):
    meta = {}
    data = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = None
        for row in reader:
            if not row:
                continue
            key = row[0].strip()
            if key.lower() == "!series_matrix_table_begin":
                header = next(reader)
                samples = [x.strip().strip('"') for x in header[1:]]
                break
            if key.startswith("!Sample_"):
                meta.setdefault(key, []).append([x.strip().strip('"') for x in row[1:]])
        if header is None:
            raise ValueError("Series matrix table start was not found")
        for row in reader:
            if not row or row[0].strip().lower() == "!series_matrix_table_end":
                break
            probe = row[0].strip().strip('"')
            if probe in {p for ps in PROBES.values() for p in ps}:
                vals = [float(x.strip().strip('"')) for x in row[1:]]
                if len(vals) != len(samples):
                    raise ValueError(f"Sample count mismatch for probe {probe}")
                data[probe] = vals

    annotations = {s: {} for s in samples}
    titles = meta.get("!Sample_title", [[]])[0]
    for i, s in enumerate(samples):
        if i < len(titles):
            annotations[s]["title"] = titles[i]
    for key, rows in meta.items():
        if key.startswith("!Sample_characteristics_ch1"):
            for vals in rows:
                for i, val in enumerate(vals[:len(samples)]):
                    if ":" in val:
                        label, value = val.split(":", 1)
                        annotations[samples[i]][label.strip().lower()] = value.strip()
    return samples, data, annotations

def auc(labels, scores):
    pos = [s for y, s in zip(labels, scores) if y == 1]
    neg = [s for y, s in zip(labels, scores) if y == 0]
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))

def percentile(vals, q):
    vals = sorted(vals)
    k = (len(vals) - 1) * q
    lo, hi = math.floor(k), math.ceil(k)
    return vals[lo] if lo == hi else vals[lo] * (hi - k) + vals[hi] * (k - lo)

with REFERENCE.open(encoding="utf-8-sig", newline="") as fh:
    ref = {r["gene"]: float(r["beta"]) for r in csv.DictReader(fh, delimiter="\t")}
if set(ref) != set(PROBES):
    raise ValueError("Frozen coefficient/reference genes do not match the platform mapping")
samples, probe_data, metadata = read_matrix(MATRIX)
missing = [p for ps in PROBES.values() for p in ps if p not in probe_data]
if missing:
    raise ValueError("Mapped Final8 probes missing from matrix: " + ", ".join(missing))
if len(samples) != 99:
    raise ValueError(f"Expected 99 samples from GSE84846; found {len(samples)}")

gene_values = {}
for gene, ps in PROBES.items():
    gene_values[gene] = [statistics.median(probe_data[p][i] for p in ps) for i in range(len(samples))]
means = {g: statistics.mean(v) for g, v in gene_values.items()}
sds = {g: statistics.stdev(v) for g, v in gene_values.items()}
if any(sd == 0 or not math.isfinite(sd) for sd in sds.values()):
    raise ValueError("A gene has zero/nonfinite cohort SD")
z = {g: [(x - means[g]) / sds[g] for x in gene_values[g]] for g in PROBES}
scores = [sum(z[g][i] * ref[g] for g in PROBES) for i in range(len(samples))]

fields = []
for s in samples:
    m = metadata[s]
    status = m.get("suvival", "")
    recurrence = m.get("recurrence", "")
    fields.append((status, recurrence))

with OUT_SCORES.open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["GSM", "suvival_literal_code", "recurrence_literal_code", "Final8_score_cohort_z"])
    for i, s in enumerate(samples):
        w.writerow([s, fields[i][0], fields[i][1], f"{scores[i]:.12g}"])

rng = random.Random(SEED)
lines = [
    "GSE84846 Final8 frozen-coefficient binary-code screen",
    f"Input matrix: {MATRIX}",
    f"Samples: {len(samples)}; all 8 genes covered; probe rows={sum(map(len, PROBES.values()))}",
    "Score: within-cohort sample-SD gene z scores × frozen TCGA coefficients; LRRC38 duplicate probes summarized by per-sample median.",
    "Endpoints retain literal GEO codes; no code meaning is assigned; no time-to-event analysis.",
    f"Bootstrap: stratified percentile CI; replicates={BOOTSTRAPS}; seed={SEED}.",
    f"Per-sample output: {OUT_SCORES.name}",
]
for label, field in (("suvival", "suvival"), ("recurrence", "recurrence")):
    valid = [(int(fields[i][0 if field == "suvival" else 1]), scores[i]) for i in range(len(samples))
             if fields[i][0 if field == "suvival" else 1] in ("0", "1")]
    ys, ss = zip(*valid)
    by = {g: [s for y, s in valid if y == g] for g in (0, 1)}
    point = auc(ys, ss)
    bs = []
    for _ in range(BOOTSTRAPS):
        draw = [(g, rng.choice(by[g])) for g in (0, 1) for __ in range(len(by[g]))]
        bs.append(auc([y for y, _ in draw], [s for _, s in draw]))
    lines += ["", f"Field {label}: N={len(valid)}; code0={len(by[0])}; code1={len(by[1])}",
              f"AUC(code1 vs code0)={point:.6f}; bootstrap 95% CI [{percentile(bs,.025):.6f}, {percentile(bs,.975):.6f}]"]
lines += ["", "No conclusion about deaths, recurrence-free survival, discrimination at a time horizon, calibration, or net benefit is supported by this analysis."]
OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(OUT_SUMMARY.read_text(encoding="utf-8"))
