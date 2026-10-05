"""Synthetic end-to-end demonstration; not market evidence."""
from fractions import Fraction as F
from pipeline import Pipeline
SAFE=[[0,8,4],[12,0,18],[19,15,0]]
CRITICAL=[[0,5,7],[2,0,20],[6,9,0]]

def run():
    model=Pipeline();out=[]
    for t,kind in enumerate(['safe','safe','critical','critical','critical','invalid','safe','safe']):
        w=SAFE if kind=='safe' else CRITICAL
        eps=F(10000) if kind=='invalid' else F(1,1000000)
        r=model.step(w,eps,t,t,model.version)
        out.append({'t':t,'input':kind,'label':r['evidence']['label'],'output':r['output'],'state':r['state']})
    return out
if __name__=='__main__':
    import json
    print(json.dumps(run(),indent=2))
