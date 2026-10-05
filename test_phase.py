from random import Random
from fractions import Fraction as F
import cmath,math
from detector import Q,charge,plaquettes,flux_base

def random_q(rng):
    while True:
        z=Q(rng.randint(-5,5),rng.randint(-5,5))
        if not z.zero():return z

def run():
    rng=Random(10205);count=0;maxres=0.
    for n in [2,3,5]:
        for _ in range(15):
            zx=[[random_q(rng) for j in range(n)] for i in range(n)]
            zy=[[random_q(rng) for j in range(n)] for i in range(n)]
            try:c=charge(zx,zy)
            except ValueError:continue
            assert 2*abs(c)<n*n
            val=sum(cmath.phase(p.complex()) for p in plaquettes(zx,zy))/(2*math.pi)
            assert abs(val-c)<1e-12;maxres=max(maxres,abs(val-c))
            h=[[random_q(rng) for j in range(n)] for i in range(n)]
            ax=[[h[i][j].conj()*zx[i][j]*h[(i+1)%n][j] for j in range(n)] for i in range(n)]
            ay=[[h[i][j].conj()*zy[i][j]*h[i][(j+1)%n] for j in range(n)] for i in range(n)]
            assert charge(ax,ay)==c
            assert charge([[z.conj() for z in row] for row in zx],[[z.conj() for z in row] for row in zy])==-c
            count+=1
    zx,zy=flux_base();assert charge(zx,zy)==1
    assert charge([[Q(1)]*2 for _ in range(2)],[[Q(1)]*2 for _ in range(2)])==0
    z=[[Q(1)]*2 for _ in range(2)];z[0][0]=Q(-1)
    try:charge(z,[[Q(1)]*2 for _ in range(2)])
    except ValueError:pass
    else:raise AssertionError('Branch boundary accepted')
    # Intermediate running product reaches negative real for the flux fixture.
    ps=plaquettes(zx,zy);assert (ps[0]*ps[1]).im==0 and (ps[0]*ps[1]).re<0
    return {'admissible_random_lattices':count,'gauge_checks':count,'orientation_checks':count,
            'independent_float_max_residual':maxres,'exact_flux_charge':1,'boundary_rejected':True}
