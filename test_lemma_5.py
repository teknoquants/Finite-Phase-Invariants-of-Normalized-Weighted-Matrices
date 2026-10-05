"""Lemma 5: exact half-integer rounding, including negative integers."""
from fractions import Fraction as F
from math import floor

def recover(x):return floor(x+F(1,2))
def run():
    count=0
    integers=list(range(-30,31))+[-10**50,10**50]
    errors=[F(-1,2)+F(1,10**20),F(-1,3),F(0),F(1,3),F(1,2)-F(1,10**20)]
    for c in integers:
        for e in errors:
            assert abs(e)<F(1,2) and recover(F(c)+e)==c;count+=1
        assert recover(F(c)+F(1,2))==c+1
        assert recover(F(c)-F(1,2))==c # asymmetric tie convention
        assert recover(F(c)+F(1,2)+F(1,1000))!=c
        assert recover(F(c)-F(1,2)-F(1,1000))!=c
    # Banker's rounding is not the displayed floor convention at ties.
    assert round(F(1,2))!=recover(F(1,2))
    return {'strict_interior_cases':count,'boundary_and_outside_cases':4*len(integers),
            'arithmetic':'exact rational','note':'This auxiliary lemma is not used by the exact winding detector.'}
if __name__=='__main__':print(run())
