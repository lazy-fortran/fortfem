#!/usr/bin/env python3
"""Run unchanged original FortFEM tests one process at a time; no aggregate RSS."""
import argparse, csv, hashlib, json, os, pathlib, statistics, subprocess
p=argparse.ArgumentParser()
p.add_argument('--baseline-build',type=pathlib.Path,required=True)
p.add_argument('--native-build',type=pathlib.Path,required=True)
p.add_argument('--measure',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--repeats',type=int,default=9)
a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
logs=a.output/'logs'; logs.mkdir(exist_ok=True)
cases=['sparse_direct_backend','sparse_direct_derivatives','sparse_direct_complex_derivatives','sparse_constrained_derivatives','capi_solver']
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C')
rows=[]
with (a.output/'trials.jsonl').open('w') as raw:
    for trial in range(-1,a.repeats):
        for case in cases:
            variants=['baseline','native'] if trial%2==0 else ['native','baseline']
            for variant in variants:
                build=getattr(a,variant+'_build'); binary=build/('test_'+case)
                log=logs/f'{case}-{variant}-{trial}.txt'
                result=subprocess.run([str(a.measure),str(binary),str(log)],env=env,capture_output=True,text=True)
                record=json.loads(result.stdout)
                record.update(case=case,variant=variant,trial=trial,warmup=trial<0,
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
                raw.write(json.dumps(record,sort_keys=True)+'\n'); raw.flush()
                if result.returncode or record['exit_status'] or record['wall_s']>=1:
                    raise SystemExit(f'Failed or >=1s case: {record}; log={log}')
                if trial>=0: rows.append(record)
summary=[]
for case in cases:
    record={'case':case,'n_per_variant':a.repeats,'status':'PASS'}
    for metric in ['wall_s','user_s','system_s','cpu_s','maxrss_kib']:
        values={variant:[r[metric] for r in rows if r['case']==case and r['variant']==variant] for variant in ['baseline','native']}
        for variant in values:
            record[f'{variant}_{metric}']=statistics.median(values[variant])
            record[f'{variant}_{metric}_min']=min(values[variant]);record[f'{variant}_{metric}_max']=max(values[variant])
        baseline=record[f'baseline_{metric}'];native=record[f'native_{metric}']
        record[f'delta_{metric}']=native-baseline
        record[f'delta_{metric}_pct']=100*(native-baseline)/baseline if baseline else None
    summary.append(record)
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (a.output/'summary.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
print(json.dumps({'measured_trials':len(rows),'warmups':len(cases)*2,'all_pass':True,'summary':str(a.output/'summary.json')}))
