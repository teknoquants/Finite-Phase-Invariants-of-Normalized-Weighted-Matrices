"""Bind exact certification to the gate under a versioned fixed configuration."""
from dataclasses import replace
from hashlib import sha256
from detector import certify,flux_base,rat
from gating import Gate,Policy

class Pipeline:
    def __init__(self,node_ids=('a','b','c'),kappa=4,base=None,policy=Policy(),bound_contract='external-l1-v1'):
        self.node_ids=tuple(node_ids);self.kappa=rat(kappa)
        base=flux_base() if base is None else base
        self.base=tuple(tuple(tuple(row) for row in a) for a in base)
        if len(set(self.node_ids))!=len(self.node_ids) or any(not isinstance(x,str) or not x for x in self.node_ids):
            raise ValueError('Unique node IDs required')
        if not bound_contract:raise ValueError('Measurement bound contract must be versioned')
        material=(self.node_ids,str(self.kappa),repr(self.base),policy.h_base,policy.h_on,policy.h_off,
                  str(policy.min_dwell),str(policy.max_age),str(policy.max_gap),tuple(sorted(policy.critical)),
                  policy.version,bound_contract,'finite-phase-v2')
        self.version=sha256(repr(material).encode()).hexdigest()
        self.gate=Gate(replace(policy,version=self.version))
    def step(self,weights,epsilon,t,sample_time,declared_version):
        try:
            evidence=certify(weights,epsilon,self.kappa,self.base,self.node_ids)
        except (TypeError,ValueError) as e:
            evidence={'accepted':False,'label':None,'reason':str(e)}
        b=self.gate.step(t,sample_time,evidence['label'],evidence['accepted'],declared_version)
        return {'version':self.version,'evidence':evidence,'output':b,'state':self.gate.state,
                'scope':'conditional on external measurement-error bound'}
