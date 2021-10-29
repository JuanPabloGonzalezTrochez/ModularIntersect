# tst5
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=-93*x*y^2*z^2-86*y^3*z^2-79*x^4-67*x^3*y-24*x*z^3+63*y^3; 
b:=50*y^3+33*z^3-76*y^2+26*y*z+71*z^2+30; 
f:=-67*x^3*y^2-75*x*y^2*z^2-73*x^3*z+70*x*y^2*z+52*y^3*z-14*x^3;

#d=(5,3,5) 
# Time Intersect: 0.173
# Time Modular Intersect: 0.745
