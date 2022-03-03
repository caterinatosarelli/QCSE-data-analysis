from hypothesis import given, strategies as st
from Functions import polynom2

@given(x=st.floats(min_value=(-8),max_value=(8),allow_nan=None), a0=st.floats(allow_nan=None), 
       a1=st.floats(allow_nan=None), a2=st.floats(min_value=(0),exclude_min=True),allow_nan=None) # coeff. of x^2 can't be 0
def test_symmetry(x,a0,a1,a2):
    # Convex parabola
    par_r = polynom2(x, a0, a1, a2)
    # Symmetric points with respect to axis :
    x_l = -(a1/a2) - x 
    par_l = polynom2(x_l, a0, a1, a2)
    # Test symmetry
    assert par_r == par_l 
    # Consider also concave parabola
    a2n = - a2 
    x_l = -(a1/a2n) - x 
    par_1 = polynom2(x, a0, a1, a2n)
    par_2 = polynom2(x, a0, a1, a2n)
    assert par_1 == par_2
    
@given(x=st.floats(min_value=(-8),max_value=(8),allow_nan=None), a0=st.floats(allow_nan=None), 
       a1=st.floats(allow_nan=None), a2=st.floats(min_value=(0),exclude_min=True,allow_nan=None))
def test_sign(x,a0,a1,a2):
    # Test that all the points are positive(negative) if a2>0(<0) and vertex y coord>0(<0)
    y_vertex = - (a1**2 - 4*a2*a0)/(4*a2)
    par = polynom2(x,a0,a1,a2)
    if y_vertex > 0:
        assert par > 0
    a2n = - a2
    par_n = polynom2(x,a0,a1,a2n)
    if y_vertex < 0:
        assert par_n < 0
        

