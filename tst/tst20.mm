# tst15
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=90*x^4*y*z^3-42*y^8-95*x^4*y^3+98*x^4*z^3+86*y*z^5-47*y^2*z; 
b:=-28*y^8+71*y^6*z^2+48*y*z^5+31*y^2-52*z+87; 
f:=-6*x*y^4*z^4+32*x^7-57*x^4*y*z-69*y^4*z^2-87*y^5-84*x*z^3;

#d=(8,7,8) 
# Time Intersect: 63.850
# Time Modular Intersect: 23.085
