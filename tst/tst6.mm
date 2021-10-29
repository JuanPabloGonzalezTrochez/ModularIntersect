# tst6
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t :=-37*x^4*y+61*x^3*y*z+5*x^2*z^3-31*x^2*z^2-22*y^2*z^2-47*y;
b:= 57*y^3*z-78*y^2*z^2+46*y*z^3+83*z^4-73*z^2+52*y;
f:= -46*x^4*z+45*y^5-2*y^4*z+39*x^4-70*x^3-53*z^2;

#d=(5,4,5) 
# Time Intersect: 0.368
# Time Modular Intersect: 1.820
