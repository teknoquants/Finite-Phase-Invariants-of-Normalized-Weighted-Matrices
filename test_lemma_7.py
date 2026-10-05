"""Lemma 7: weighted quantiles, integrated ES, and composed gate scaling."""
from fractions import Fraction as F
from pipeline import Pipeline
from demo import SAFE,CRITICAL

def law(rows):
    out={}
    for x,p in rows:
        x,p=F(x),F(p)
        if p<=0:raise ValueError('Weights must be positive')
        out[x]=out.get(x,F(0))+p
    if sum(out.values(),F(0))!=1:raise ValueError('Weights must sum to one')
    return sorted(out.items())
def quantile(rows,alpha):
    alpha=F(alpha)
    if not 0<alpha<1:raise ValueError('Require 0 < alpha < 1')
    c=F(0)
    for x,p in law(rows):
        c+=p
        if c>=alpha:return x
    raise AssertionError('Unreachable')
def es(rows,alpha):
    alpha=F(alpha)
    if not 0<alpha<1:raise ValueError('Require 0 < alpha < 1')
    c=F(0);integral=F(0)
    for x,p in law(rows):
        nxt=c+p
        integral+=x*max(F(0),nxt-max(c,alpha));c=nxt
    return integral/(1-alpha)
def threshold_es(rows,alpha):
    q=quantile(rows,alpha)
    return q+sum(p*max(x-q,F(0)) for x,p in law(rows))/(1-alpha)
def run():
    laws=[[(1,F(1,4)),(2,F(1,4)),(4,F(1,4)),(8,F(1,4))],
          [(-9,F(1,10)),(-2,F(2,5)),(0,F(1,5)),(7,F(3,10))],
          [(0,F(9,10)),(10,F(1,10))],[(5,F(1,3)),(5,F(2,3))],
          [(-3,F(1,7)),(2,F(2,7)),(2,F(1,7)),(11,F(3,7))],[(0,F(1))]]
    scales=[F(1,1000),F(1,7),F(1),F(3),F(1000)];count=0;independent=0
    for rows in laws:
        alphas={F(1,100),F(1,4),F(1,2),F(9,10),F(99,100)};c=F(0)
        for _,p in law(rows):
            c+=p
            for a in [c-F(1,1000),c,c+F(1,1000)]:
                if 0<a<1:alphas.add(a)
        for a in sorted(alphas):
            q=quantile(rows,a);tail=es(rows,a)
            assert sum(p for x,p in law(rows) if x<q)<a<=sum(p for x,p in law(rows) if x<=q)
            assert tail==threshold_es(rows,a) and tail>=q;independent+=1
            assert quantile(list(reversed(rows)),a)==q and es(list(reversed(rows)),a)==tail
            for s in scales:
                scaled=[(s*x,p) for x,p in rows]
                assert quantile(scaled,a)==s*q and es(scaled,a)==s*tail;count+=1
    # An atom at VaR requires fractional tail allocation, not conditioning >= VaR.
    atom=laws[2];assert quantile(atom,F(3,4))==0 and es(atom,F(3,4))==4
    naive=sum(F(x)*p for x,p in atom);assert naive==1 and naive!=es(atom,F(3,4))
    # Exact analytic uniform-family values provide non-discrete examples.
    uniform=0
    for left,right in [(F(-4),F(7)),(F(0),F(2)),(F(-8),F(-2))]:
        for a in [F(1,10),F(1,2),F(99,100)]:
            q=left+(right-left)*a;tail=left+(right-left)*(1+a)/2
            for s in scales:
                assert s*left+(s*right-s*left)*a==s*q
                assert s*left+(s*right-s*left)*(1+a)/2==s*tail;uniform+=1
    # End-to-end scale change: graph, uncertainty and magnitude rescale together.
    seq=[SAFE,SAFE,CRITICAL,CRITICAL,CRITICAL,SAFE,SAFE];pipeline_cases=0
    for s in scales:
        a=Pipeline();b=Pipeline()
        for t,w in enumerate(seq):
            u=a.step(w,F(1,1000000),t,t,a.version)
            v=b.step([[s*x for x in row] for row in w],s*F(1,1000000),t,t,b.version)
            assert u['output']==v['output'] and u['evidence']['label']==v['evidence']['label']
            payout=u['output']*es(laws[0],F(3,4))
            scaled_payout=v['output']*es([(s*x,p) for x,p in laws[0]],F(3,4))
            assert scaled_payout==s*payout;pipeline_cases+=1
    return {'weighted_scale_cases':count,'independent_tail_identity_checks':independent,
            'uniform_scale_cases':uniform,'pipeline_scale_steps':pipeline_cases,
            'atom_counterexample':{'integrated_ES':'4','naive_conditional_mean':'1'},
            'scope':'Finite weighted laws and analytic uniform examples; not all probability laws.'}
if __name__=='__main__':print(run())
