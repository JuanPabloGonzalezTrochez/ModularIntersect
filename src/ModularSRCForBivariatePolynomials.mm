#<-----------------------------------------------------------------------------
# Function: { subresultant_chain_2_var}
# Briefly: { return the subresultant chain of two bivariate polynomials}
# Calling sequence:
# { TRDsub_resultant_chain_bivariate_modP(A, B, u, v, m) }
# Input:
# { A : polynomial or a modp1(polynomial)}
# { B : polynomial or a modp1(polynomial)}
# { u : symbol }
# { v : symbol }
# { m : prime }
# Output:
# { return the subresultant chain of `A` and `B` in Zm[x] }
#>-----------------------------------------------------------------------------
ModularSRCForBivariatePolynomials := module()
local ModuleApply := proc(A::polynom, 
                          B::polynom, 
                          u::symbol, 
                          v::symbol, 
                          m::prime,
                          {mode :: identical(:-polynom, :-modp1) := ':-polynom'},
                          $)

   local a,ai,am,da,b,bi,bm,db,d,x,y,i,j,l,n,r,t,la,lb,D,L,M,N,R,
         cofBound,rowBound,colBound,bezBound,delta,Delta,alpha,beta,
         src, S1, numCoeS1, idx, z, S1Int, lcoeff_bi_to_power;

   if not type([A,B],list(':-polynom'(':-integer',[u,v]))) or A=0 or B=0 then 
      error "inputs must be non-zero polynomials in Z[%1,%2]",u,v; 
   end if;

   if degree(A,u)=0 and degree(B,u)=0 then 
      return subresultant_chain(A, B, v, m);
   end if;

   if degree(A,v)=0 and degree(B,v)=0 then 
      return [A]; 
   end if;

   # MBM: Syl(A,B,x) must not have a zero column for bound computations, hence
   if ldegree(A,u)>0 and ldegree(B,u)>0 then #change, A and B mul of u
      return 0;
   end if;

   userinfo(2,resultant,"bivariate case: using the modular method");
   cofBound := 2*`resultant/bound`(A,B,u); # 2 is for positive/negative coefficients

   da := degree(A,u); a := [seq(coeff(A,u,i), i=0..da)]; la := lcoeff(a[da+1],v); 
   db := degree(B,u); b := [seq(coeff(B,u,i), i=0..db)]; lb := lcoeff(b[db+1],v);

   bezBound := degree(A,{u,v})*degree(B,{u,v}); # Bezout bound   
   rowBound := da*degree(B,v) + db*degree(A,v);
   colBound :=
   add( max( 0, seq( degree(a[1+d],v), d=max(0,i-db+1)..min(da,i) ),
              seq( degree(b[1+d],v), d=max(0,i-da+1)..min(db,i) ) ), i=0..da+db-1 );
   D := min(rowBound,colBound,bezBound); # degree bound: D >= degree(R,v)

   bezBound := ldegree(A,{u,v})*ldegree(B,{u,v});
   rowBound := da*ldegree(B,v) + db*ldegree(A,v);
   colBound :=
   add( min( D, seq( ldegree(a[1+d],v), d=max(0,i-db+1)..min(da,i) ),
              seq( ldegree(b[1+d],v), d=max(0,i-da+1)..min(db,i) ) ),
      i=0..da+db-1 );
   L := max(rowBound,colBound,bezBound); # low degree bound: L <= ldegree(R,v)

   userinfo(3,RegularChains,
      "Coefficient bound: log[10](Height(R)) <= %1 digits.",1+ilog10(cofBound));
   userinfo(3,RegularChains,
      "Degree bounds: deg[%1](R) <= %2 and ldeg[%1] >= %3.",v,D,L);

   delta := 0;
   d := D; # initial degree bound
   l := L; # and low degree bound


   if modp(la,m)=0 or modp(lb,m)=0 then
      error "m must not divide the leading Coefficient of A and B";
   end if;

   am := [seq( modp1(ConvertIn(t,v), m), t=a )]; # am = a mod m
   bm := [seq( modp1(ConvertIn(t,v), m), t=b )]; # bm = b mod m 

   N := d-l+delta;
   alpha :=1;
   beta := NEXTPRIMELEM(m,10);
   i := 0;
   numCoeS1 := 0;

   # Evaluate at N+1 "random" points from Zm
   while i <= N do 
      alpha := modp(beta*alpha,m); 
      if alpha=1 then 
         break 
      end if;

      ai := [seq(modp1(Eval(am[i], alpha), m), i = 1 .. nops(am))];
      ai := modp1(ConvertIn(ai,u),m);
      if modp1(Degree(ai), m) < da then 
         next 
      end if;

      bi := [seq(modp1(Eval(bm[i], alpha), m), i = 1 .. nops(bm))];
      bi := modp1(ConvertIn(bi,u),m); 
      if modp1(Degree(bi), m) < db then 
         next 
      end if;

      src := subresultant_chain(ai, bi, u, m);

      r:= src[1];
      S1 := subresultant_chain_index_one(src,m);
    
      if l>0 then 
         r := modp(r/(alpha &^ l),m); 
      else 
         r := modp(r,m); 
      end if;

      x[i], y[i] := alpha, r;

      z[i] := [coeffs(S1, u)];

      if numCoeS1 < nops(z[i]) then 
         ASSERT( numCoeS1 = 0, "Number of coeficients of S1 changed");
         numCoeS1 := nops(z[i]);
      end if;


      i := i+1;
   end do;

   if alpha=1 then # we ran out of points 
      userinfo(3,resultant,"The algorithm ran out of points while Interpolating."); 
      return FAIL;
   end if;

   r := modp1(Interp( [seq(x[i], i=0..N)], [seq(y[i], i=0..N)], v ), m);
   S1Int := InterMult(x, z, v, numCoeS1, m, N);

   if mode = ':-polynom' then
      r := modp1(ConvertOut(Shift(r,l),v), m); 
      S1Int := ConvOutList(S1Int, numCoeS1, r, v, m);
   else
      r := modp1(Shift(r,l), m);
   end if;

   if r=0 then 
      (d,l) := (L-1,L) 
   else 
      (d,l) := (degree(r,v),ldegree(r,v)) 
   end if;

   userinfo(3,resultant,sprintf("Degree correction from %d,%d to %d,%d.",L,D,l,d));
   (j,delta,R,M) := (1,1,r,m); # M=m


   # Termination
   if degree(A,u)<>1 or degree(B,u)<>1 then
      if mode = ':-polynom' then
         S1 := PolynomialTools:-FromCoefficientList(S1Int, u);
      else
         S1 := S1Int;
      end if;
      return [R, S1]; 
   else
      return [R, undefined];
   end if;
   
end proc;

# Computes and remembers the first primitive element > b in Zm.
local NEXTPRIMELEM := proc(m,b) 
   local f,p,beta,good; option remember;

   f := ifactors(m-1);
   f := [seq(p[1],p=op(2,f))];   
   beta := b+1;

   do good := true;
      for p in f while good do
         if modp(power(beta,iquo(m-1,p)),m)=1 then 
            good := false;
         end if;
      end do;

      if good then 
         return beta; 
      end if;
      beta := beta+1;
   end do;
end proc;

# Interpolate a list of modp1 polynomials.
local InterMult := proc(x, z, v::name, n::nonnegint, m, N, $)
   local inter, j, i;
   local x_list := [seq(x[i], i=0..N)];
   inter := [seq(modp1(Interp( 
               x_list, [seq(z[i][j], i=0..N)], v ), m), j=1..n)];

   return inter;

end proc;

# Convert out a list of modp1 polynomials.
local ConvOutList := proc(S1Int::list, n::nonnegint, r, v, m, $)
   local S1Conv, i;
   S1Conv := [];

   for i from 1 to n do
      S1Conv :=[op(S1Conv), modp1(ConvertOut(S1Int[i],v), m)];
   end do;

   return S1Conv;
end proc;


local subresultant_chain;
local subresultant_chain_index_one;
$include "src/resultant_chain.mm"

end module;