from fractions import Fraction as F
from random import Random
from detector import Q,normalize,embed,certify,charge

def run():
    rng=Random(62);w=[[0,2,1],[3,0,4],[2,1,0]];eps=F(1,100000)
    cert=certify(w,eps);assert cert['accepted'];n=3;checks=0
    for _ in range(40):
        # Perturb a random directed edge; all data and path points are exact.
        i,j=rng.sample(range(n),2);sign=rng.choice([-1,1]);v=[list(r) for r in w];v[i][j]+=sign*eps
        for k in range(11):
            u=[[F(w[a][b])+F(k,10)*(v[a][b]-w[a][b]) for b in range(n)] for a in range(n)]
            zx,zy,_=embed(u);assert charge(zx,zy)==cert['label'];checks+=1
    assert not certify(w,13)['accepted'] # epsilon >= observed mass
    assert not certify(w,F(5))['accepted'] # margin fails well below total mass
    assert certify(w,0)['accepted']
    for s in [F(1,7),F(9)]:
        c=certify([[s*x for x in row] for row in w],s*eps)
        assert c['accepted'] and c['label']==cert['label'] and c['slack']==cert['slack']
    # Attached Proposition 9: vertices of the normalized simplex and interior fixtures
    # all retain the baseline charge when kappa < 1/8.
    constant_cases=0
    for i in range(3):
        for j in range(3):
            if i==j:continue
            v=[[0]*3 for _ in range(3)];v[i][j]=1
            assert certify(v,0,kappa=F(1,16))['label']==1;constant_cases+=1
    for _ in range(12):
        v=[[0 if i==j else rng.randint(1,10) for j in range(3)] for i in range(3)]
        assert certify(v,0,kappa=F(1,16))['label']==1;constant_cases+=1
    # Search fixed exact fixtures demonstrating sensitivity beyond total mass.
    examples={}
    for _ in range(2000):
        v=[[0 if i==j else rng.randint(1,20) for j in range(3)] for i in range(3)]
        c=certify(v,F(1,1000000))
        if c['accepted']:examples.setdefault(c['label'],v)
        if 0 in examples and any(k!=0 for k in examples):break
    assert 0 in examples and any(k!=0 for k in examples)
    return {'certified_path_samples':checks,'label':cert['label'],
            'rho':str(cert['rho']),'beta_lower':str(cert['beta']),
            'exact_certificate_slack':str(cert['slack']),'scale_certificate_checks':2,
            'rejected_bounds':2,'constant_phase_cases':constant_cases,'distinct_graph_labels':{str(k):v for k,v in examples.items()}}
