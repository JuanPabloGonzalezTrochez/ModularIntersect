# tst11
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:= -75*x^4*z-38*x*y^4-9*x*y^3*z+25*y^3+7*x^2-43*x*y;
b:= 93*y^4*z-37*y^3*z^2+65*y*z^2-87*y*z+46*z^2+46;
f:= -17*x*y^3+68*x*y*z^2-27*y*z^3+42*y^3+23*x^2+76*x*z;

#d=(5,4,5) 
# Time Intersect: 0.279
# Time Modular Intersect: 0.722