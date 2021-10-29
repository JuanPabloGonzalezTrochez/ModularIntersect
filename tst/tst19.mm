# tst15
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=6*x^7*z-23*x^2*y^2*z^4+88*x*y^6*z+41*z^8+83*x^2*y^4*z-23*y*z^3; 
b:=-92*y^7*z-85*y*z^7+60*y^3*z^2-81*y*z^2+54*z^3+22*z; 
f:=-13*x^4*y^2*z^2+88*x*y^7-7*y^7*z+16*z^7+31*x^2*y-67*y^3;

#d=(8,7,8) 
# Time Intersect: 43.324
# Time Modular Intersect: 18.497
