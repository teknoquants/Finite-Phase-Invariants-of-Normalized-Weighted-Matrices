"""Exact rational finite-phase detector. No floating-point arithmetic in certificates.
An accepted result is CONDITIONAL on the user-supplied measurement error bound.
It certifies a finite model label, never a financial event or continuum invariant.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from math import gcd


def rat(x):
    if isinstance(x, (int,F)): return F(x)
    if isinstance(x,str): return F(x)
    raise TypeError('Use integers, Fractions or decimal strings, not binary floats')

@dataclass(frozen=True)
class Q:
    re:F
    im:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'re',rat(self.re));object.__setattr__(self,'im',rat(self.im))
    def __add__(self,b):
        if not isinstance(b,Q):b=Q(b)
        return Q(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self):return Q(-self.re,-self.im)
    def __sub__(self,b):return self+-b
    def __mul__(self,b):
        if not isinstance(b,Q):b=Q(b)
        return Q(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def conj(self):return Q(self.re,-self.im)
    def zero(self):return self.re==0 and self.im==0
    def primitive(self):
        if self.zero():return self
        den=self.re.denominator*self.im.denominator
        a=int(self.re*den);b=int(self.im*den);g=gcd(a,b)
        return Q(a//g,b//g) # positive scaling preserves the argument
    def complex(self):return complex(float(self.re),float(self.im))


def dimensions(zx,zy):
    n=len(zx)
    if n<2 or len(zy)!=n or any(len(r)!=n for r in zx+zy):
        raise ValueError('Require equal square lattices with n >= 2')
    if any(not isinstance(z,Q) or z.zero() for a in (zx,zy) for r in a for z in r):
        raise ValueError('All links must be nonzero rational complex numbers')
    return n


def plaquettes(zx,zy):
    n=dimensions(zx,zy)
    return [zx[i][j]*zy[(i+1)%n][j]*zx[i][(j+1)%n].conj()*zy[i][j].conj()
            for i in range(n) for j in range(n)]


def side(z):
    # Sign of principal argument in (-pi,pi], with Arg(negative real)=pi.
    if z.im>0 or (z.im==0 and z.re<0):return 1
    if z.im<0:return -1
    return 0


def charge(zx,zy):
    """Exact winding accumulation; no atan, pi, eigenvectors or rounding."""
    ps=plaquettes(zx,zy);s=Q(1);w=0
    for p in ps:
        if p.im==0 and p.re<0:raise ValueError('Plaquette on excluded branch boundary')
        p=p.primitive();v=(s*p).primitive()
        sa,sb,sv=side(s),side(p),side(v)
        if sa>0 and sb>0 and sv<=0:w+=1
        if sa<0 and sb<0 and sv>=0:w-=1
        s=v
    if s.im!=0 or s.re<=0:raise ArithmeticError('Periodic cancellation identity failed')
    return w


def normalize(weights):
    n=len(weights)
    if n<2 or any(len(r)!=n for r in weights):raise ValueError('Square graph, N >= 2 required')
    w=[[rat(x) for x in r] for r in weights]
    if any(x<0 for r in w for x in r):raise ValueError('Negative weights')
    if any(w[i][i]!=0 for i in range(n)):raise ValueError('Diagonal must be zero')
    s=sum((x for r in w for x in r),F(0))
    if s<=0:raise ValueError('Zero total exposure')
    return [[x/s for x in r] for r in w],s


def flux_base(n=2):
    if n!=2:raise ValueError('The provided exact fixture uses n=2')
    return [[Q(1),Q(0,-1)],[Q(1),Q(0,-1)]],[[Q(1),Q(1)],[Q(1),Q(-1)]]


def embed(weights,kappa=F(4),base=None,node_ids=None):
    x,total=normalize(weights);N=len(x);kappa=rat(kappa)
    if node_ids is not None:
        if len(node_ids)!=N or any(not isinstance(v,str) or not v for v in node_ids) or len(set(node_ids))!=N:
            raise ValueError('Unique nonempty string node IDs required')
        order=sorted(range(N),key=lambda i:node_ids[i])
        x=[[x[i][j] for j in order] for i in order]
    if kappa<=0:raise ValueError('kappa must be positive')
    if base is None:base=flux_base()
    zx=[[Q(z.re,z.im) for z in r] for r in base[0]]
    zy=[[Q(z.re,z.im) for z in r] for r in base[1]]
    n=dimensions(zx,zy)
    if N*(N-1)>4*n*n:raise ValueError('Embedding capacity exceeded; specify a larger fixed base')
    # Inject each directed edge into a distinct real coordinate of the link array.
    values=[x[i][j] for i in range(N) for j in range(N) if i!=j]
    for slot,value in enumerate(values):
        link,part=divmod(slot,2);axis,pos=divmod(link,n*n);i,j=divmod(pos,n)
        arr=zx if axis==0 else zy
        arr[i][j]=arr[i][j]+(Q(kappa*value) if part==0 else Q(0,kappa*value))
    return zx,zy,total


def certify(weights,epsilon,kappa=F(4),base=None,node_ids=None):
    """epsilon bounds TOTAL entrywise l1 error against unknown true nonnegative W.
    Returns UNKNOWN if assumptions or sufficient margins fail. accepted is
    conditional on that externally supplied bound being valid, not estimated.
    """
    epsilon=rat(epsilon);kappa=rat(kappa)
    if epsilon<0:raise ValueError('Negative error bound')
    try:
        zx,zy,total=embed(weights,kappa,base,node_ids)
        if epsilon>=total:return {'accepted':False,'reason':'mass bound','label':None}
        c=charge(zx,zy)
        rho=min(max(abs(z.re),abs(z.im)) for a in (zx,zy) for row in a for z in row)
        beta=min(F(1) if p.re>=0 else abs(p.im)/(abs(p.re)+abs(p.im)) for p in plaquettes(zx,zy))
        d=2*kappa*epsilon/total
        slack=beta*rho-8*d
        return {'accepted':slack>0,'reason':'conditional bound' if slack>0 else 'margin bound',
                'label':c if slack>0 else None,'observed_label':c,
                'rho':rho,'beta':beta,'link_error_bound':d,'slack':slack,'total':total}
    except ValueError as e:
        return {'accepted':False,'reason':str(e),'label':None}
