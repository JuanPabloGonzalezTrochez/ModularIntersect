# tst8
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=-50*x^3*y*z+33*x*z^4-10*x*y^2*z+36*y^2*z^2+54*x*y^2-5*y^3;
b:= 86*y^3*z-17*z^4+51*y^2*z+43*y^2+19*y*z+44*z;
f:= 89*y^5-43*x*y^3+28*x^2*y-23*z^3+59*x*z+82*y^2;

#d=(5,4,5) 
# Time Intersect: 0.113
# Time Modular Intersect: 1.342
