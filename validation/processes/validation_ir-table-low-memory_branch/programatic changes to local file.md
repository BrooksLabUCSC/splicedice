change local file



```bash
python3 - <<'EOF'
path = "ir_table.py"
with open(path) as f:
    src = f.read()

old = '''def getJunctions(samples, coverageDirectory, annotated, args):
    from multiprocessing import Pool
    worker_args = [(sample, coverageDirectory, annotated, args.allJunctions) for sample in samples]
    with Pool(args.numThreads) as pool:
        results = pool.map(_getJunctionsForSample, worker_args)
    return set().union(*results)'''

new = '''def getJunctions(samples, coverageDirectory, annotated, args):
    from multiprocessing import Pool
    worker_args = [(sample, coverageDirectory, annotated, args.allJunctions) for sample in samples]
    with Pool(args.numThreads) as pool:
        junctions = set()
        for result in pool.imap_unordered(_getJunctionsForSample, worker_args):
            junctions |= result
    return junctions'''
    
    
if old in src:
    src = src.replace(old, new)
    with open(path, "w") as f:
        f.write(src)
    print("SUCCESS - getFilteredJunctions replaced")
else:
    print("FAILED - getFilteredJunctions not found")
EOF

```



```bash
python3 - <<'EOF'
path = "ir_table.py"
with open(path) as f:
    src = f.read()

old2 = '''def getFilteredJunctions(samples, coverageDirectory, annotated, args):
    import time
    from multiprocessing import Pool
    t = time.time()
    junctions = getJunctions(samples, coverageDirectory, annotated, args)
    print(f"getJunctions complete: {len(junctions)} junctions. {time.time()-t:.1f}s")
    worker_args = [(sample, coverageDirectory, junctions, args.RSDthreshold) for sample in samples]
    with Pool(args.numThreads) as pool:
        results = pool.map(_filterJunctionsForSample, worker_args)
    filtered_junctions = set().union(*results)
    print(f"RSD filtering complete: {len(filtered_junctions)} junctions retained. {time.time()-t:.1f}s")
    return filtered_junctions'''

new2 = '''def getFilteredJunctions(samples, coverageDirectory, annotated, args):
    import time
    from multiprocessing import Pool
    t = time.time()
    junctions = getJunctions(samples, coverageDirectory, annotated, args)
    print(f"getJunctions complete: {len(junctions)} junctions. {time.time()-t:.1f}s")
    worker_args = [(sample, coverageDirectory, junctions, args.RSDthreshold) for sample in samples]
    with Pool(args.numThreads) as pool:
        filtered_junctions = set()
        for result in pool.imap_unordered(_filterJunctionsForSample, worker_args):
            filtered_junctions |= result
    print(f"RSD filtering complete: {len(filtered_junctions)} junctions retained. {time.time()-t:.1f}s")
    return filtered_junctions'''

if old2 in src:
    src = src.replace(old2, new2)
    with open(path, "w") as f:
        f.write(src)
    print("SUCCESS - getFilteredJunctions replaced")
else:
    print("FAILED - getFilteredJunctions not found")
EOF
```

sparse `getInclusionCounts`:

```bash
python3 - <<'EOF'
path = "ir_table.py"
with open(path) as f:
    src = f.read()

old = '''def getInclusionCounts(filename, annotated=None):
    import pandas as pd
    df = pd.read_csv(filename, sep="\\t", index_col=0)
    counts = df.to_dict(orient="index")
    counts = {sample: {junction: counts[junction][sample] for junction in counts} for sample in df.columns}
    return counts'''

new = '''def getInclusionCounts(filename, annotated=None):
    counts = {}
    with open(filename) as f:
        samples = f.readline().rstrip().split("\\t")[1:]
        for sample in samples:
            counts[sample] = {}
        for line in f:
            row = line.rstrip().split("\\t")
            junction = row[0]
            for i, val in enumerate(row[1:]):
                v = float(val)
                if v != 0.0:
                    counts[samples[i]][junction] = v
    return counts'''

if old in src:
    src = src.replace(old, new)
    with open(path, "w") as f:
        f.write(src)
    print("SUCCESS - getInclusionCounts replaced")
else:
    print("FAILED - getInclusionCounts not found")
EOF

```

```bash
python3 - <<'EOF'
path = "ir_table.py"
with open(path) as f:
    src = f.read()

old = '''def getInclusionCounts(filename, annotated=None):
    counts = {}
    with open(filename) as f:
        samples = f.readline().rstrip().split("\\t")[1:]
        for sample in samples:
            counts[sample] = {}
        for line in f:
            row = line.rstrip().split("\\t")
            junction = row[0]
            for i, val in enumerate(row[1:]):
                v = float(val)
                if v != 0.0:
                    counts[samples[i]][junction] = v
    return counts'''

new = '''def getInclusionCounts(filename, annotated=None):
    counts = {}
    with open(filename) as f:
        samples = f.readline().rstrip().split("\\t")[1:]
        for sample in samples:
            counts[sample] = {}
        for line in f:
            row = line.rstrip().split("\\t")
            junction = row[0]
            for i, val in enumerate(row[1:]):
                counts[samples[i]][junction] = float(val)
    return counts'''

if old in src:
    src = src.replace(old, new)
    with open(path, "w") as f:
        f.write(src)
    print("SUCCESS - getInclusionCounts updated")
else:
    print("FAILED - text not found")
EOF
```

