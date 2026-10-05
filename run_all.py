"""Reproducible proof checks. Refuses optimized execution with disabled assertions."""
import argparse,json,platform,sys
from pathlib import Path
if not __debug__:
    raise SystemExit('Run without -O; assertions are required.')
import test_algebra,test_phase,test_certificate,test_gating,test_risk
import test_lemma_5,test_lemma_7,test_lemma_14,test_lemma_15

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    report={'python':platform.python_version(),'arithmetic':'Exact rational core; floating point used only for an independent phase comparison',
            'paper':'Finite Phase Invariants of Normalized Weighted Matrices',
            'package_version':'3.0',
            'scope':'Finite exact checks mapped to all 15 numbered results; not universal proofs',
            'results':{}}
    for mod in [test_algebra,test_phase,test_certificate,test_gating,test_risk,test_lemma_5,test_lemma_7,test_lemma_14,test_lemma_15]:
        report['results'][mod.__name__]={'status':'PASS',**mod.run()}
    s=json.dumps(report,indent=2)
    if args.output:args.output.write_text(s+'\n')
    print(s)
if __name__=='__main__':main()
