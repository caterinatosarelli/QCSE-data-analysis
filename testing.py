# from hypothesis import given, strategies as st
from Functions import polynom2
from Functions import gaussian_gen
from Functions import peak_lim

#%% 
def test_polynom2_points():
    """Test that the quadratic function passes from 
    - intercept, when x=0
    - origin, when c=0 and x=0
    - a generic point. """
    assert polynom2(0,1.2,0.005,0.0001)== 1.2 # intercept
    assert polynom2(0,0,0.001,0.5) == 0 # origin
    assert polynom2(2,1,1,2) == 11 # random point

def test_polynom2_sign():
    """ Tests sign of function when there are no intersection with x axis,
    for a>0 and a>0 parabola."""
    assert polynom2(1,3,2,1) > 0 # a>0, delta<0 
    assert polynom2(1,-1,5,-8) < 0 # a<0 , delta <0
    
def test_polynom2_symmetry():
    """ Test symmetry of the function w.r.t.
    - origin, when y = a x^2
    - it's axis, where axis has equation x=-b/2a """
    assert polynom2(-2,0,0,1) == polynom2(2,0,0,1) # wrt origin
    assert polynom2(0.7,0.1,-2,1) == polynom2(1.3,0.1,-2,1) # wrt axis 

#%%        
def test_gaussgen_symmetry():
    """Test symmetry of the function wrt to
    - origin, if mean = 0
    - x = mean """ 
    assert gaussian_gen(0.1,50,0.05,0,1) == gaussian_gen(-0.1,50,0.05,0,1) #wrt to origin if mu=0
    assert gaussian_gen(0.1,50,0.05,0.5,1)== gaussian_gen(0.9,50,0.05,0.5,1) #wrt to mu
    
#%% 
def test_peak_lim():
    """ Tests that the conversion between wavelength and energy works."""
    wav = [930.6106180,930.6550800,930.6995420,930.7440040,930.7884660,
           930.8329280,930.8773900]
    en = [1.3322887,1.3322250,1.3321614,1.3320977,1.3320341,1.3319705,1.3319068]
    assert peak_lim(wav,en,1.33212,0.00011)[0] == 1.3321614
      
#%%
# @given(x=st.floats(min_value=(-8),max_value=(8),allow_nan=None), a0=st.floats(allow_nan=None), 
#        a1=st.floats(allow_nan=None), a2=st.floats(min_value=(0),exclude_min=True),allow_nan=None) # coeff. of x^2 can't be 0
# def test_symmetry(x,a0,a1,a2):
#     # Convex parabola
#     par_r = polynom2(x, a0, a1, a2)
#     # Symmetric points with respect to axis :
#     x_l = -(a1/a2) - x 
#     par_l = polynom2(x_l, a0, a1, a2)
#     # Test symmetry
#     assert par_r == par_l 
#     # Consider also concave parabola
#     a2n = - a2 
#     x_l = -(a1/a2n) - x 
#     par_1 = polynom2(x, a0, a1, a2n)
#     par_2 = polynom2(x, a0, a1, a2n)
#     assert par_1 == par_2
    
# @given(x=st.floats(min_value=(-8),max_value=(8),allow_nan=None), a0=st.floats(allow_nan=None), 
#        a1=st.floats(allow_nan=None), a2=st.floats(min_value=(0),exclude_min=True,allow_nan=None))
# def test_sign(x,a0,a1,a2):
#     # Test that all the points are positive(negative) if a2>0(<0) and vertex y coord>0(<0)
#     y_vertex = - (a1**2 - 4*a2*a0)/(4*a2)
#     par = polynom2(x,a0,a1,a2)
#     if y_vertex > 0:
#         assert par > 0
#     a2n = - a2
#     par_n = polynom2(x,a0,a1,a2n)
#     if y_vertex < 0:
#         assert par_n < 0
        

    
