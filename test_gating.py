from fractions import Fraction as F
from itertools import product
from gating import Gate,Policy

def execute(seq,policy=Policy()):
    g=Gate(policy);out=[]
    for t,x in enumerate(seq):out.append(g.step(t,t,x,x is not None,policy.version))
    return out,g

def run():
    p=Policy(h_base=2,h_on=3,h_off=2,min_dwell=2)
    out,g=execute([0,0,1,1,0,0,0,None,1,1],p)
    assert out==[0,0,0,0,0,0,1,None,1,0]
    # Exhaustive short streams test no activation without enough stress evidence.
    cases=0
    for seq in product([None,0,1],repeat=8):
        out,g=execute(seq,p);run=0
        for k,(x,b) in enumerate(zip(seq,out)):
            run=run+1 if x==0 else 0
            if b==1 and (k==0 or out[k-1]==0):assert run>=p.h_on
        oracle=Gate(p);safe_run=0;armed_at=None
        for k,x in enumerate(seq):
            old=oracle.state
            safe_run=safe_run+1 if x==1 else 0
            oracle.step(k,k,x,x is not None,p.version)
            if old!='ARMED' and oracle.state=='ARMED':armed_at=k
            if old=='ARMED' and oracle.state!='ARMED':
                assert safe_run>=p.h_off and k-armed_at>=p.min_dwell
            if old=='ARMED' and x is None:assert oracle.state=='ARMED'
        cases+=1
    # Invalid evidence while armed does not disarm latent protection.
    out,g=execute([1,1,0,0,0,None],p);assert out[-1] is None and g.state=='ARMED'
    assert execute([0]*20,p)[1].state=='BOOT'
    assert execute([1,1,0,0,1,0,0,1],p)[1].state=='READY'
    # Recovery is also consecutive, with a minimum elapsed-time hold.
    p2=Policy(h_base=1,h_on=1,h_off=2,min_dwell=5)
    out,g=execute([1,0,1,1,1,1,1],p2);assert out==[0,1,1,1,1,1,0]
    checks=0
    for kind in ['stale','future','duplicate','version','unaccepted','gap']:
        g=Gate(p);g.step(0,0,1,True,p.version)
        args={'t':F(1),'sample_time':F(1),'label':1,'accepted':True,'version':p.version}
        if kind=='stale':args['sample_time']=F(-2)
        if kind=='future':args['sample_time']=F(2)
        if kind=='duplicate':args['sample_time']=F(0)
        if kind=='version':args['version']='changed'
        if kind=='unaccepted':args['accepted']=False
        if kind=='gap':args.update(t=F(10),sample_time=F(10))
        assert g.step(**args) is None and g.state=='BOOT';checks+=1
    g=Gate(p);g.step(1,1,1,True,p.version)
    try:g.step(1,1,1,True,p.version)
    except ValueError:pass
    else:raise AssertionError('Duplicate decision time accepted')
    prefix=[1,1,0,0]
    assert execute(prefix+[0,1],p)[0][:4]==execute(prefix+[None,0],p)[0][:4]
    from demo import run as demo_run
    from pipeline import Pipeline
    demo=demo_run();assert [r['output'] for r in demo]==[0,0,0,0,1,None,1,0]
    assert Pipeline(kappa=4).version!=Pipeline(kappa=5).version
    return {'end_to_end_steps':len(demo),'configuration_hash_check':True,'exhaustive_streams':cases,'stream_length':8,'invalid_evidence_cases':checks,
            'initial_critical_abstention':True,'latent_hold_verified':True,'dwell_verified':True}
