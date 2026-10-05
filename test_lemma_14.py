"""Lemma 14: formal local polynomial check of metric independence.
Validates the bare CS integrand algebra, not global geometry or a path integral.
"""
from fractions import Fraction as F
from itertools import permutations
class Poly:
    def __init__(self,d=None):self.d={k:F(v) for k,v in (d or {}).items() if v}
    @classmethod
    def var(cls,name):return cls({(name,):1})
    def __add__(self,b):
        if not isinstance(b,Poly):b=Poly({():b})
        d=self.d.copy()
        for k,v in b.d.items():d[k]=d.get(k,F(0))+v
        return Poly(d)
    __radd__=__add__
    def __mul__(self,b):
        if not isinstance(b,Poly):b=Poly({():b})
        d={}
        for a,x in self.d.items():
            for c,y in b.d.items():
                k=tuple(sorted(a+c));d[k]=d.get(k,F(0))+x*y
        return Poly(d)
    __rmul__=__mul__
    def diff(self,name):
        d={}
        for k,v in self.d.items():
            n=k.count(name)
            if n:
                a=list(k);a.remove(name);a=tuple(a);d[a]=d.get(a,F(0))+n*v
        return Poly(d)
    def evaluate(self,values):
        total=F(0)
        for names,c in self.d.items():
            for n in names:c*=F(values[n])
            total+=c
        return total

def matmul(a,b):return [[sum((a[i][k]*b[k][j] for k in range(2)),Poly()) for j in range(2)] for i in range(2)]
def tr(a):return a[0][0]+a[1][1]
def sign(p):return (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
def run():
    A=[[[Poly.var(f'A{m}{i}{j}') for j in range(2)] for i in range(2)] for m in range(3)]
    dA=[[[[Poly.var(f'D{u}{v}{i}{j}') for j in range(2)] for i in range(2)] for v in range(3)] for u in range(3)]
    bilinear=Poly();cubic=Poly()
    for i,j,k in permutations(range(3)):
        bilinear+=sign((i,j,k))*tr(matmul(A[i],dA[j][k]))
        cubic+=F(2,3)*sign((i,j,k))*tr(matmul(matmul(A[i],A[j]),A[k]))
    integrand=bilinear+cubic
    assert bilinear.d and cubic.d and integrand.d
    metrics=[f'g{i}{j}' for i in range(3) for j in range(i,3)]
    for g in metrics:assert not integrand.diff(g).d
    names={x for k in integrand.d for x in k};values={name:F((i%7)-3) for i,name in enumerate(sorted(names))}
    for scale in [1,2,7]:
        v={**values,**{g:F(0) for g in metrics},'g00':F(scale),'g11':F(2*scale),'g22':F(3*scale),'g01':F(1,10)}
        assert integrand.evaluate(v)==integrand.evaluate(values)
    # Nonzero constant non-Abelian connection on the unit periodic three-torus.
    # Ax=E12, Ay=E21, Az=diag(1,-1); dA=0. Cubic coefficient integrates to 4.
    fixture={name:F(0) for name in names}
    fixture.update({'A001':F(1),'A110':F(1),'A200':F(1),'A211':F(-1)})
    assert bilinear.evaluate(fixture)==0 and cubic.evaluate(fixture)==4
    # Negative control: a metric-dependent density does not pass this check.
    bad=Poly.var('g00')*integrand
    assert bad.diff('g00').d==integrand.d and bad.diff('g00').d
    return {'symbolic_monomials':len(integrand.d),'metric_derivatives_zero':len(metrics),
            'metric_substitutions':3,'nonabelian_torus_integrand':'4',
            'negative_control_detected':True,
            'scope':'Formal local algebra and a constant torus example; no quantum metric-independence claim.'}
if __name__=='__main__':print(run())
