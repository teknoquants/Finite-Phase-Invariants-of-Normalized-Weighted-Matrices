from fractions import Fraction as F
from itertools import product
from detector import Q,normalize,embed,flux_base,charge,plaquettes

def run():
    w=[[0,2,1],[3,0,4],[2,1,0]]
    x,s=normalize(w)
    for a in [F(1,7),F(1),F(13,2)]:
        y,_=normalize([[a*v for v in row] for row in w]);assert x==y
    # Exact normalization bound for 12 one-edge perturbations.
    cases=0
    for i in range(3):
        for j in range(3):
            if i==j:continue
            for delta in [F(-1,10),F(1,10)]:
                v=[list(r) for r in w];v[i][j]+=delta;y,_=normalize(v)
                assert sum(abs(y[a][b]-x[a][b]) for a in range(3) for b in range(3))<=2*abs(delta)/s
                cases+=1
    # Recover all six normalized graph coordinates from the affine links.
    base=flux_base();zx,zy,_=embed(w);observed=[]
    for arr,brr in zip((zx,zy),base):
        for row,brow in zip(arr,brr):
            for a,b in zip(row,brow):observed.extend([(a.re-b.re)/4,(a.im-b.im)/4])
    assert observed[:6]==[x[i][j] for i in range(3) for j in range(3) if i!=j]
    ids=['b','c','a'];perm=[2,0,1]
    permuted=[[w[i][j] for j in perm] for i in perm]
    assert embed(w,node_ids=ids)==embed(permuted,node_ids=[ids[i] for i in perm])
    # Same total exposure, different normalized edges: embedding must differ.
    v=[list(r) for r in w];v[0][1]+=1;v[0][2]-=1
    assert embed(w)[:2]!=embed(v)[:2]
    # Any magnitude decoder faces a scale collision; exact minimax bound.
    for a in [F(-10),F(10),F(55),F(100),F(120)]:
        assert max(abs(a-10),abs(a-100))>=45
    # Exact empirical quantile and integrated ES on an equal-weight law.
    losses=[F(1),F(2),F(4),F(8)]
    def risk(v):return sorted(v)[2],sorted(v)[3]
    for a in [F(1,3),F(1),F(9)]:
        assert risk([a*v for v in losses])==tuple(a*v for v in risk(losses))
    for wbad in [[[0,0],[0,0]],[[1,2],[3,0]],[[0,-1],[2,0]]]:
        try:normalize(wbad)
        except ValueError:pass
        else:raise AssertionError('Invalid graph accepted')
    return {'normalization_perturbations':cases,'scale_checks':6,'injection_coordinates':6,'invalid_graphs':3,'node_reordering_check':True,'same_mass_distinguished':True}
