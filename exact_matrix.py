"""Small exact rational complex matrices for finite-dimensional channel checks."""
from fractions import Fraction as F
from detector import Q

def zeros(n,m):return [[Q(0) for _ in range(m)] for _ in range(n)]
def eye(n):return [[Q(int(i==j)) for j in range(n)] for i in range(n)]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(a,s):return [[x*s for x in r] for r in a]
def adj(a):return [[a[j][i].conj() for j in range(len(a))] for i in range(len(a[0]))]
def mul(a,b):
    if len(a[0])!=len(b):raise ValueError('Dimension mismatch')
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q(0)) for j in range(len(b[0]))] for i in range(len(a))]
def trace(a):return sum((a[i][i] for i in range(len(a))),Q(0))
def kron(a,b):return [[x*y for x in rowa for y in rowb] for rowa in a for rowb in b]
def ptrace_a(r,da,db):return [[sum((r[a*db+i][a*db+j] for a in range(da)),Q(0)) for j in range(db)] for i in range(db)]
def channel(r,ks,db):
    out=zeros(len(r),len(r))
    for k in ks:
        u=kron(k,eye(db));out=add(out,mul(mul(u,r),adj(u)))
    return out
