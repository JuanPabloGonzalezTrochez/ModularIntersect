# tst9
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=-70*x^5-4*x^2*z^3-16*y^3*z^2-18*x^4+96*x*z-11;
b:= 94*y^3*z+56*y^2*z-13*y^2-52*y*z-81*z^2+59;
f := 69*x^2*y^3+23*x^2*z^3+56*x*y^2*z^2-9*x*y*z^3+71*y^4+62*y^2*z;

#d=(5,4,5) 
# Time Intersect: 0.378
# Time Modular Intersect: 1.087
