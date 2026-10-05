from fractions import Fraction as F
from itertools import product

def run():
    # Exact finite-horizon nuisance-run probability, iid law as one admissible case.
    T=8;h=3;p=F(1,5);prob=F(0)
    for bits in product([0,1],repeat=T):
        if any(all(bits[k:k+h]) for k in range(T-h+1)):
            s=sum(bits);prob+=p**s*(1-p)**(T-s)
    bound=(T-h+1)*p**h
    assert prob<=bound
    dependent=F(0)
    for bits in product([0,1],repeat=T):
        mass=F(1);previous=0
        for b in bits:
            q=F(1,5) if previous else F(1,10)
            mass*=q if b else (1-q);previous=b
        if any(all(bits[k:k+h]) for k in range(T-h+1)):dependent+=mass
    assert dependent<=bound
    # Conditional mean risk identity on a finite joint distribution (Z,C,Y).
    rows=[(0,0,0,F(1,4)),(0,1,1,F(1,4)),(1,0,0,F(1,8)),(1,0,1,F(1,8)),(1,1,1,F(1,4))]
    def means(key):
        a={};b={}
        for z,c,y,w in rows:
            k=key(z,c);a[k]=a.get(k,F(0))+w*y;b[k]=b.get(k,F(0))+w
        return {k:a[k]/b[k] for k in a}
    m=means(lambda z,c:z);n=means(lambda z,c:(z,c))
    r0=sum(w*(y-m[z])**2 for z,c,y,w in rows)
    r1=sum(w*(y-n[z,c])**2 for z,c,y,w in rows)
    gain=sum(w*(n[z,c]-m[z])**2 for z,c,y,w in rows)
    assert r0-r1==gain and gain>0
    # C is determined by full inputs X=(Z,C): adding it cannot change the oracle.
    full={(z,c):n[z,c] for z,c,y,w in rows};assert all(full[z,c]==n[z,c] for z,c,y,w in rows)
    # Detection delay bound for every event offset on a rational grid.
    Delta=F(1);latency=F(1,4);h=3
    for i in range(20):
        tau=F(i,20);first=int(tau)+1;t_arm=F(first+h-1)+latency
        assert t_arm<=tau+h*Delta+latency
    return {'run_probability':str(prob),'run_union_bound':str(bound),'dependent_run_probability':str(dependent),
            'bayes_risk_before':str(r0),'bayes_risk_after':str(r1),'risk_gain':str(gain),
            'delay_offsets':20,'universal_proof_claimed':False}
