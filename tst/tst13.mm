# tst13
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:= 4*x^9-40*x^5*y^2*z+6*x^3*y^3*z+27*x*y^6+68*x*y^3*z^2-11*z^5;
b:= -33*y^8*z+8*y^5*z^2-69*y^4*z^2-34*z^6-58*y^5-53*y*z^2;
f:= -7*x^3*y^2*z^4-50*y^4*z^5-70*x^3*y^5+19*x*y^5-5*y^3*z+48*x;

#d=(8,7,8) 
# Time Intersect: 301.062
# Time Modular Intersect: 27.999
