"""Lemma 15: exact full reduced-state equality for local TP channels."""
from fractions import Fraction as F
from random import Random
from detector import Q
from exact_matrix import *

def density(rng,n):
    b=[[Q(F(rng.randint(-3,3),5),F(rng.randint(-3,3),7)) for _ in range(n)] for _ in range(n)]
    r=mul(b,adj(b));t=trace(r);assert t.im==0 and t.re>0
    return scale(r,1/t.re) # positive semidefinite by construction

def run():
    rng=Random(1514);cases=0;completeness=0
    for da,db in [(2,2),(2,3),(3,2)]:
        channels=[[eye(da)]]
        phase=eye(da);phase[0][0]=Q(F(3,5),F(4,5));channels.append([phase])
        dephase=[]
        for a in range(da):
            k=zeros(da,da);k[a][a]=Q(1);dephase.append(k)
        channels.append(dephase)
        if da==2:
            channels.append([[[Q(1),Q(0)],[Q(0),Q(F(4,5))]],[[Q(0),Q(F(3,5))],[Q(0),Q(0)]]])
            paulis=[eye(2),[[Q(0),Q(1)],[Q(1),Q(0)]],[[Q(0),Q(0,-1)],[Q(0,1),Q(0)]],[[Q(1),Q(0)],[Q(0),Q(-1)]]]
            channels.append([scale(k,F(1,2)) for k in paulis])
        for ks in channels:
            total=zeros(da,da)
            for k in ks:total=add(total,mul(adj(k),k))
            assert total==eye(da);completeness+=1
            for _ in range(5):
                r=density(rng,da*db);out=channel(r,ks,db)
                assert out==adj(out) and trace(out)==Q(1)
                assert ptrace_a(out,da,db)==ptrace_a(r,da,db);cases+=1
    bell=zeros(4,4)
    for i in [0,3]:
        for j in [0,3]:bell[i][j]=Q(F(1,2))
    project=[[Q(1),Q(0)],[Q(0),Q(0)]]
    conditional=channel(bell,[project],2);prob=trace(conditional).re
    assert prob==F(1,2)
    posterior=scale(conditional,1/prob)
    assert ptrace_a(posterior,2,2)!=ptrace_a(bell,2,2)
    # Restoring both outcomes is trace-preserving and restores the same marginal.
    other=[[Q(0),Q(0)],[Q(0),Q(1)]]
    assert ptrace_a(channel(bell,[project,other],2),2,2)==ptrace_a(bell,2,2)
    return {'exact_joint_state_channel_cases':cases,'kraus_completeness_checks':completeness,
            'bipartite_dimensions':[[2,2],[2,3],[3,2]],'all_marginal_residuals':'exactly zero',
            'postselection_negative_control':True,'unconditional_recovery_check':True,
            'scope':'Finite-dimensional examples; not an exhaustive or infinite-dimensional proof.'}
if __name__=='__main__':print(run())
