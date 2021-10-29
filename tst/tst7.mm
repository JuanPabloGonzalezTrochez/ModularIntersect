# tst7
######################
m := 469762049;
R := PolynomialRing([x,y,z], m);

t:=-6*x^3*z^2-40*x*y^2*z^2+79*x*y+73*z^2+42*x-6*z;
b:= -91*y^4-97*y^2*z^2-63*z^4+81*y^2*z+89*z^3-81*z^2;
f:=-38*y^5-28*x^2*y^2-31*x^3-63*x*y^2-15*y*z+57*x;

#d=(5,4,5) 
# Time Intersect: 0.148
# Time Modular Intersect: 0.823
