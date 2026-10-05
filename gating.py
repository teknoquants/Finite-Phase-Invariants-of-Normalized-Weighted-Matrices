"""Persistent causal gate with explicit UNKNOWN and frozen latent protection."""
from dataclasses import dataclass
from fractions import Fraction as F
from detector import rat

@dataclass(frozen=True)
class Policy:
    version:str='v2-fixed'
    h_base:int=2
    h_on:int=3
    h_off:int=2
    min_dwell:F=F(2)
    max_age:F=F(1)
    max_gap:F=F(2)
    critical:frozenset=frozenset({0})
    def __post_init__(self):
        if any(type(h)!=int or h<1 for h in (self.h_base,self.h_on,self.h_off)):
            raise ValueError('Persistence counts must be positive integers')
        for a in ('min_dwell','max_age','max_gap'):
            object.__setattr__(self,a,rat(getattr(self,a)))
            if getattr(self,a)<0:raise ValueError('Negative time parameter')
        if self.max_gap<=0:raise ValueError('max_gap must be positive')
        if not isinstance(self.critical,frozenset) or any(type(c)!=int for c in self.critical):
            raise ValueError('Critical labels must be an immutable set of integers')
        if not self.version:raise ValueError('Version required')

class Gate:
    def __init__(self,policy=Policy()):
        self.p=policy;self.state='BOOT';self.baseline=self.on=self.off=0
        self.active_since=None;self.last_decision=None;self.last_sample=None
    def unknown(self):
        self.baseline=self.on=self.off=0
        if self.state!='ARMED':self.state='BOOT'
        return None
    def step(self,t,sample_time,label,accepted,version):
        t=rat(t);sample_time=rat(sample_time)
        if self.last_decision is not None and t<=self.last_decision:
            raise ValueError('Decision times must strictly increase; no state change')
        gap=self.last_decision is not None and t-self.last_decision>self.p.max_gap
        self.last_decision=t
        fresh=(0<=t-sample_time<=self.p.max_age and
               (self.last_sample is None or sample_time>self.last_sample))
        if gap or not fresh or (type(accepted)!=bool or not accepted) or type(label)!=int or version!=self.p.version:
            return self.unknown()
        self.last_sample=sample_time
        stress=label in self.p.critical
        if self.state=='BOOT':
            self.baseline=self.baseline+1 if not stress else 0
            if self.baseline>=self.p.h_base:self.state='READY';self.baseline=0
        elif self.state=='READY':
            self.on=self.on+1 if stress else 0
            if self.on>=self.p.h_on:
                self.state='ARMED';self.active_since=t;self.on=0;self.off=0
        else:
            self.off=min(self.off+1,self.p.h_off) if not stress else 0
            if self.off>=self.p.h_off and t-self.active_since>=self.p.min_dwell:
                self.state='READY';self.off=0;self.active_since=None
        return int(self.state=='ARMED')
