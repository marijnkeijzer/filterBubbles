"""Portable reimplementation of the retained notebook's one-topic ABM.
See README.md for legacy coverage and published/notebook design differences.
No defSim dependency; RNG sequences are not identical to the historical runs.
"""
import argparse, csv, gzip, json, math, random, statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def coverage_legacy(x):
    # Deliberately preserve the original CoverageReporter indexing, including its bug.
    x=sorted(x);r=.5/len(x);lo=[v-r for v in x];hi=[v+r for v in x];uncovered=0
    for i in range(len(x)):
        if i==1:
            if lo[i]>0: uncovered+=lo[i]
        elif i<len(x)-1:
            if lo[i]>hi[i-1]: uncovered+=lo[i]-hi[i-1]
        elif hi[i]<1: uncovered+=1-hi[i]
    return 1-uncovered

def coverage_interval(x):
    # Union of opinion-centered intervals of total width 1/N, clipped to [0,1].
    r=.5/len(x);end=0.;total=0.
    for v in sorted(x):
        left=max(0.,v-r);right=min(1.,v+r)
        total+=max(0.,right-max(end,left));end=max(end,right)
    return total

def simulate(seed,bubble,steps,alpha=.5,gamma=2.,calibration=None,n=50):
    rng=random.Random(seed);x=[rng.random() for _ in range(n)]
    for _ in range(steps):
        i=rng.randrange(n)
        neighbors=sorted((j for j in range(n) if j!=i),key=lambda j:abs(x[j]-x[i]))[:bubble]
        if not neighbors: continue
        j=rng.choice(neighbors);delta=x[j]-x[i]
        if calibration is None: a=alpha;b=alpha*gamma
        else:
            a=rng.gauss(calibration['alpha_mean'],calibration['alpha_sd'])
            b=rng.gauss(calibration['beta_mean'],calibration['beta_sd'])
        # Equivalent to a*(1-(b/a)*abs(delta))*delta, without division by a.
        x[i]=min(1.,max(0.,x[i]+a*delta-b*abs(delta)*delta))
    mu=statistics.mean(x)
    return {'spread':max(x)-min(x),'dispersion':2*statistics.mean(abs(v-mu) for v in x),
            'coverage_legacy':coverage_legacy(x),'coverage_interval':coverage_interval(x)}

def calibration_from_draws(folder):
    with gzip.open(folder/'m11.csv.gz','rt') as f: rows=list(csv.DictReader(f))
    a=[float(r['alpha']) for r in rows];b=[float(r['gamma']) for r in rows]
    return dict(alpha_mean=statistics.mean(a),alpha_sd=statistics.stdev(a),beta_mean=statistics.mean(b),beta_sd=statistics.stdev(b))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile',choices=['illustration','fixed','calibrated','caption'],default='calibrated')
    p.add_argument('--seed',type=int,default=20261117)
    p.add_argument('--smoke',action='store_true')
    p.add_argument('--draws',type=Path,default=ROOT/'data'/'posterior')
    args=p.parse_args()
    if args.profile=='illustration': grid=[(b,1.5+.1*k) for b in [20,40] for k in range(16)];steps=1000
    elif args.profile=='fixed':grid=[(b,1.58) for b in range(51)];steps=1000
    elif args.profile=='calibrated':grid=[(b,None) for b in [10,20,30,40]];steps=5000
    else:grid=[(b,None) for b in range(51)];steps=1000
    reps=10
    if args.smoke:grid=grid[:2];steps=25;reps=2
    cal=calibration_from_draws(args.draws)
    seed_rng=random.Random(args.seed);rows=[]
    for bubble,gamma in grid:
        for rep in range(reps):
            seed=seed_rng.randrange(2**63)
            values=simulate(seed,bubble,steps,gamma=gamma,calibration=cal if gamma is None else None)
            rows.append(dict(bubble=bubble,gamma=gamma,replicate=rep+1,seed=seed,steps=steps,**values))
    out=ROOT/'output';out.mkdir(exist_ok=True)
    name='abm_'+args.profile+('_smoke' if args.smoke else '')
    with (out/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (out/(name+'_settings.json')).write_text(json.dumps(dict(profile=args.profile,seed=args.seed,steps=steps,repetitions=reps,agents=50,calibration=cal,coverage='both legacy notebook and interval union'),indent=2)+'\n')
    print(f'Wrote {len(rows)} runs: output/{name}.csv')

if __name__=='__main__': main()
