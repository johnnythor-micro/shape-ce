import sys; sys.path.insert(0,'.')
import numpy as np, glob, os
from shapece.io import read_abif, get_channel
from shapece import preprocess as pp
from shapece import peaks, align, reactivity as rx
U="/mnt/user-data/uploads"
F={os.path.basename(f).split('_')[0]:f for f in glob.glob(U+"/*_20260701_*.fsa")}
# RX=1M7 lanes, BG=vehicle lanes (37C); SIZE=DATA205, data=DATA9
def load(w): 
    e=read_abif(F[w]); return get_channel(e,9), get_channel(e,205)
lanes=['D12','E12','F12','A12','B12','C12']  # 3x 1M7 + 3x vehicle (37C)
raw={}; size={}
for w in lanes:
    d,s=load(w)
    d=pp.baseline(pp.smooth(pp.correct_saturation(d)))
    raw[w]={'RX':d}; size[w]=s
# size-standard alignment onto D12
aligned=align.size_standard_align(raw,size,reference='D12',roles=('RX',))
ref=np.mean([aligned[w]['RX'] for w in ['D12','E12','F12']],axis=0)
pos=peaks.detect_peaks(ref,lo=300,hi=3000,min_spacing=12,prominence_frac=0.03)
print("peaks detected:",len(pos))
q={w:peaks.quantify(aligned[w]['RX'],pos,mode='gaussian') for w in lanes}
# area-difference reactivity per replicate, boxplot normalized
pairs=[('D12','A12'),('E12','B12'),('F12','C12')]
profs=[]
for p,m in pairs:
    react=rx.area_difference(q[p]['area'],q[m]['area'],scale=True)
    norm,factor=rx.boxplot_normalize(react)
    profs.append(norm)
profs=np.array(profs)
from scipy.stats import pearsonr
print("replicate correlations (area-difference):",
      [f"{pearsonr(profs[i],profs[j])[0]:.2f}" for i,j in [(0,1),(0,2),(1,2)]])
print("mean normalized reactivity:",round(profs.mean(),3),"| range",round(profs.min(),2),"-",round(profs.max(),2))
print("CORE AREA-DIFFERENCE PIPELINE OK")
