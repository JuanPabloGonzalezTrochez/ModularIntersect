## lrc is a list of regular chains with polynomials
## in the polynomial ring R

constructibleSet := proc(lrc, R, h:=1)
local lrs, cs, hs;
    lrs := map(RegularSystem, lrc, R);
    cs := ConstructibleSet(lrs, R);
    hs := GeneralConstruct([h], [], R);
    return Difference(cs,hs,R);
end proc;


AreEqual := proc(cs1, cs2, R)
   return IsContained(cs1, cs2, R) and IsContained(cs2, cs1, R);
end proc;

Normalize_lrc := proc(lrc, R)
   local lnrc :=  map(TRDnormalize_zerodim_regularchain, lrc, R);
end proc;

(*
nor_inter := Normalize_lrc(Inter, R);
inter_una := EquiprojectableDecomposition(nor_inter, R);
*)

(*
R := PolynomialRing([x, y, z]);

sys := {-x^5+y^5-3*y-1, 5*y^4-3, -20*x+y-z};
dec := Triangularize(sys, R);
dec1 := Triangularize(sys, R, normalized=strongly);

AreEqual(constructibleSet(dec, R), constructibleSet(dec1, R), R);
*)

## GeneralConstruct(F, H, R)
