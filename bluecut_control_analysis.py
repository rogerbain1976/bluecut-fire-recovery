# Blue Cut Fire Recovery: Burned vs Matched Unburned Control
# Author: Roger Bain
# Input: data/BlueCut_Recovery_Points_v2.csv from BlueCut_Recovery_v2 (Earth Engine)
# Steps: drop developed/water/barren/ag land and control points burned since 1996,
# trim controls to the burn area 5th-95th percentile of elevation, slope, and pre-fire NDVI,
# weight controls to match burned points across pre-fire NDVI x elevation x slope strata,
# then compare NDVI relative to pre-fire. CIs use a 1 km spatial block bootstrap (1,000 reps, seed 1).

import pandas as pd, numpy as np
d=pd.read_csv('data/BlueCut_Recovery_Points_v2.csv').drop(columns=['system:index','.geo'])
d=d[~d.NLCD2016.isin([21,22,23,24,11,31,81,82])]
b=d[d.group=='burned'].copy()
c=d[(d.group=='control')&(d.Burned_Since_1996==0)].copy()
print('burned',len(b),'control after filters',len(c))
# overlap trim to burned 5-95 pct
for v in ['Elevation_m','Slope_deg','Prefire_NDVI']:
    lo,hi=b[v].quantile([.05,.95]); c=c[c[v].between(lo,hi)]
print('control after overlap trim',len(c))
# CEM-style weights: bins from burned quantiles
def binz(s,ref,q): return pd.cut(s,np.unique(ref.quantile(np.linspace(0,1,q+1))),include_lowest=True)
for df in (b,c):
    df['k']=(binz(df.Prefire_NDVI,b.Prefire_NDVI,5).astype(str)+'|'+binz(df.Elevation_m,b.Elevation_m,4).astype(str)+'|'+binz(df.Slope_deg,b.Slope_deg,3).astype(str))
bc=b.k.value_counts(); cc=c.k.value_counts()
common=bc.index.intersection(cc.index)
c=c[c.k.isin(common)].copy(); c['w']=c.k.map(bc/cc)
bm=b[b.k.isin(common)].copy(); bm['w']=1.0
print('matched strata',len(common),'burned in matched',len(bm),'of',len(b),'controls used',len(c))
print('balance (mean) burned vs weighted control:')
for v in ['Prefire_NDVI','Elevation_m','Slope_deg']:
    print(' ',v, round(bm[v].mean(),3), round(np.average(c[v],weights=c.w),3), '| unweighted ctrl',round(c[v].mean(),3))
P=['Prefire','Postfire','Year1','Year5','Year9','Year10']
rows=[]
for p in P:
    col=p+'_NDVI'
    bb=bm[bm[col]>-1]; cq=c[c[col]>-1]
    rows.append([p, bb[col].mean(), np.average(cq[col],weights=cq.w)])
t=pd.DataFrame(rows,columns=['period','burned','control'])
t['burned_%pre']=t.burned/t.burned[0]*100; t['control_%pre']=t.control/t.control[0]*100
t['relative_recovery']=t['burned_%pre']/t['control_%pre']*100
print(t.round(3).to_string())
# bootstrap CI by spatial 1km blocks for relative recovery Year9, Year10
rng=np.random.default_rng(1)
for df in (bm,c):
    df['blk']=(np.floor((df.lon+117.5)*92.4)).astype(int).astype(str)+'_'+(np.floor((df.lat-34.3)*111)).astype(int).astype(str)
def rr(B,C,y):
    return (B[y+'_NDVI'].mean()/B.Prefire_NDVI.mean())/(np.average(C[y+'_NDVI'],weights=C.w)/np.average(C.Prefire_NDVI,weights=C.w))*100
bb_blocks=bm.groupby('blk'); cb_blocks=c.groupby('blk')
bk=list(bb_blocks.groups); ck=list(cb_blocks.groups)
print('blocks burned',len(bk),'control',len(ck))
for y in ['Year1','Year5','Year9','Year10']:
    vals=[]
    for _ in range(1000):
        B=pd.concat([bb_blocks.get_group(k) for k in rng.choice(bk,len(bk))]); C=pd.concat([cb_blocks.get_group(k) for k in rng.choice(ck,len(ck))])
        vals.append(rr(B,C,y))
    print(y,'relative recovery',round(rr(bm,c,y),1),'95% CI',np.percentile(vals,[2.5,97.5]).round(1))
print('control NDVI change Year9->Year10 %', round((np.average(c.Year10_NDVI,weights=c.w)/np.average(c.Year9_NDVI,weights=c.w)-1)*100,1))
