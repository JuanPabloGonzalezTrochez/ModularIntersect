# tst15
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=63*x^5*z^3+54*x*y^5*z+28*x*y^5+13*y^2*z^3-77*x^2*y^2+71; 
b:=41*y^4*z^4+2*z^8-46*y^7-94*y^2*z^5-3*y^4*z+76*z^5; 
f:=63*x^3*y^5-92*x^3*y*z^3+38*y^4*z^2+13*z^6-71*x^2*y*z^2-96*x^4;

#d=(8,7,8) 
# Time Intersect: 24.418
# Time Modular Intersect: 12.920