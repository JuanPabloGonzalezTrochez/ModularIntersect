# tst14
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:= -41*x^6*y+44*x*y^6+18*y^3*z^3+28*y^2*z^4-56*x*z^4-90*x^3*z;
b:= 45*y^5*z+91*z^6+79*y^2*z^3-93*z^5+62*y*z^3-3;
f := 8*x^6*z+5*x^5*z^2-8*x^2*y*z^4-70*x^4*y^2-6*y^5*z-88*y*z^3;

#d=(7,6,7) 
# Time Intersect: 13.110
# Time Modular Intersect: 11.313