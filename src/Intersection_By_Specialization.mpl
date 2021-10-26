#  read("src/Intersection_By_Specialization.mpl");

kernelopts(opaquemodules=false):
TRDinit := RegularChains:-TRDinit;
TRDempty_regular_chain := RegularChains:-TRDempty_regular_chain:
TRDzerodim_normal_form := RegularChains:-TRDzerodim_normal_form;
TRDnormalize_zerodim_regularchain := RegularChains:-TRDnormalize_zerodim_regularchain;
kernelopts(opaquemodules=true):

#Bound of the system of equations
Bound := proc(list, $) 
	local bnd, p; 
	bnd := 1; 

	for p in list do 
		bnd := bnd*degree(p); 
	end do; 

	return 2*bnd + 1; 
end proc;

ResultantChainBySpecialization := module()
#Discriminant computed by specialization
## Assume rc is a one-dimensional chain over the integers
## in the polynomial ring R.
## f a polynomial over R with n variables
local ModuleApply := proc(f::polynom, rc::TRDrc, 
                          R::TRDring, bnd::integer,
                          $)

	local cont, lp, spe_rc, rc_normalized, f_specialized, s_v,
		  var, p, g1_v, g2_v, r_v, rational_reconstruction_without_gcd,
          pnt, deg_s_v, adi_comp1, adi_comp2, src2_v, src;

    local m := R[characteristic];
    if m = 0 then
        error "the characteristic of the field must be prime";
    end if;

    local termOrder := R[variables]; 
    local v := 0;

    local sw := 0; 
    local G1 := Array(1..0);
    local G2 := Array(1..0);
    local Res := Array(1..0);
    local Src2 := Array(1..0);
    local Res_hG1_G2:= Array(1..0);
    local Res_ht_G2:= Array(1..0);
    local good_points := Array(1..0);
    local num_coef_G2 := 0;
    local num_coef_G1 := 0;
    local deg_s := 1;
    local cont_degs_0 := 0;

	var := termOrder[nops(termOrder)]; 
	lp := ListTools:-Reverse(Equations(rc,R));

    local q := Polynomial(termOrder[1], rc, R);

    if MainVariable(f, R)<>termOrder[1] then
        error "case not covered, main variable of %1 must be %2",
                                                 f, termOrder[1];
    end if;

    local deg_f := degree(f, termOrder[1]);

    if deg_f=1 and degree(q, termOrder[1])=1 then
        return DegOneCase(f, rc, R);
    end if;

    # Initials
    local ini_1 := TRDinit(Polynomial(termOrder[1], rc, R), R);
    local ini_2 := TRDinit(Polynomial(termOrder[2], rc, R), R);

	do 
        # points to be evaluated.
        v := v+1;

        # We check if v is a good specialization point.
		spe_rc := IsGoodSpecializationPoint([var = v], lp, rc, R); 
        if spe_rc[1] = false then
            next;
        end if;

		if spe_rc[1] = true then 
            # We normalize the regular chain rc specialized in v. 
            #rc_normalized := TRDnormalize_zerodim_regularchain(spe_rc[2], R);
            rc_normalized := spe_rc[2];
     
            # We specialize f. 
			f_specialized := eval(f, [var = v]); 
            # We check our assumptions.
            if MainVariable(f_specialized, R)<>termOrder[1]    
               or degree(f_specialized, termOrder[1]) <> deg_f then
                next;
            end if;

            # We compute the subresultants of index zero and one. 
			s_v, g1_v, src2_v, adi_comp1, adi_comp2:= SubresultantOfIndexZeroAndOne(f_specialized, 
                                                             rc_normalized, 
                                                             deg_s, 
                                                             eval(ini_1, [var = v]) mod m,
                                                             R);

            r_v := src2_v[1]; 
            g2_v := subresultant_chain_index_one(src2_v,m);
            deg_s_v := modp1(Degree(s_v), m);            

            # We check our degree assumptions of the resultant
            # of f and b, s_v.
            if  deg_s_v < deg_s then
                # We ignore v.
                next;
            elif deg_s_v > deg_s then
                # We restart the computations.
                deg_s := deg_s_v;
                G1 := Array(1..0);
                G2 := Array(1..0);
                Res_hG1_G2:= Array(1..0);
                Res_ht_G2:= Array(1..0);
                Res := Array(1..0);
                good_points := Array(1..0);
                sw := 0;
                num_coef_G2 := 0;
                num_coef_G1 := 0;
                next;
            end if;
            
            # We save the coefficients of our polynomials.
            G1 ,= map(Coefficients, g1_v, m);

            G2 ,= [Coefficients(g2_v, m)];

            Res ,= Coefficients(r_v, m);

            Res_hG1_G2,=adi_comp1;
            Res_ht_G2,=adi_comp2;

            # Src2 saved.
            Src2 ,= src2_v;       

            # We save the good specialization point v.
			good_points ,= modp(v, m);

			sw := sw + 1; 
		end if; 
	until sw = bnd; 
    
    # number of coefficients
    num_coef_G1 := map(modp1, map(Degree, g1_v), m);
    num_coef_G2 := numelems(G2[-1]);

    # Border polynomial.
    local bder_2 := Resultant(Polynomial(termOrder[2], rc, R), 
                                         ini_1, termOrder[2]) mod m;    

    good_points := convert(good_points, ':-list');
    
    # We use M for the rational function reconstruction.
    local M := mul(seq(var-pnt, pnt in good_points));

    # We interpolate the resultant.
	local resultant_interpolated := modp1(':-Interp'( good_points,
                                            convert(Res, ':-list'), var), m);

    # We Convert out from modp1.
	local interpolant := PolynomialTools:-FromCoefficientList(
                                        modp1(ConvertOut(resultant_interpolated), m), var);

    # We do rational function reconstruction.
    local rational_reconstruction := Ratrecon(interpolant, M, var) mod m;   

    # We get the numerator.
    rational_reconstruction := numer(rational_reconstruction); 

    # We make the rational_reconstruction squarefree.
    local rational_reconstruction_sqfree := SquarefreePart(rational_reconstruction,
                                                            var) mod m;

    (*# We clean common factors with the border polynomial and the initials.
    rational_reconstruction_without_gcd := SimplifyResultant(bder_2,
                                                             rational_reconstruction_sqfree);

    rational_reconstruction_without_gcd := SimplifyResultant(ini_1, 
                                                             rational_reconstruction_without_gcd);
    
    rational_reconstruction_without_gcd := SimplifyResultant(ini_2, 
                                                             rational_reconstruction_without_gcd);
    *)
    rational_reconstruction_without_gcd := rational_reconstruction_sqfree;
    # We interpolate the g1_v.
    local G1_interpolated := InterpolatePolynomials(good_points, 
                                                    G1, termOrder[1..2], 
                                                    var, num_coef_G1, 
                                                    M, m,
                                                    ':-mode'=modp1);

    # We interpolate the g2_v.
    local G2_interpolated := InterpolatePolynomials(good_points, 
                                                    G2, termOrder[2], 
                                                    var, num_coef_G2, 
                                                    M, m);

    # We verify that our polynomials form a regular chain.                                                
    #local local_rc := Chain([rational_reconstruction_without_gcd], R);
    if MainVariable(G2_interpolated, R) <> termOrder[2]then
       error "case not covered, the main variable, %1, of %2 is wrong",
              MainVariable(G2_interpolated, R), G2_interpolated;
    end if; 

    if MainVariable(G1_interpolated, R) <> termOrder[1]then
       error "case not covered, the main variable, %1, of %2 is wrong",
              MainVariable(G1_interpolated, R), G1_interpolated;
    end if; 

    local ini := TRDinit(G2_interpolated, R);
    if not IsRegularMod(ini, rational_reconstruction_without_gcd, R) then
        error "case not covered, the initial %1 is not regular", ini;
    end if;

    # We check if the initial of l is regular mod {s, g}
    local Res_hG1_G2_interp := modp1(':-Interp'( good_points,
                                            convert(Res_hG1_G2, ':-list'), var), m);

    Res_hG1_G2_interp := PolynomialTools:-FromCoefficientList(
                                        modp1(ConvertOut(Res_hG1_G2_interp), m), var);


    if not IsRegularMod(Res_hG1_G2_interp, rational_reconstruction_without_gcd, R) then
        error "case not covered, the initial %1 is not regular", ini;
    end if;

    # We check if the initial of t is regular mod {s, g}
    local Res_ht_G2_interp := modp1(':-Interp'( good_points,
                                            convert(Res_ht_G2, ':-list'), var), m);

    Res_ht_G2_interp := PolynomialTools:-FromCoefficientList(
                                        modp1(ConvertOut(Res_ht_G2_interp), m), var);


    if not IsRegularMod(Res_ht_G2_interp, rational_reconstruction_without_gcd, R) then
        error "case not covered, the initial %1 is not regular", ini_1;
    end if;


    return rational_reconstruction_without_gcd, G2_interpolated, G1_interpolated; 
end proc;

# Coefficients of a modp1 polynomial.
local Coefficients::static := proc(p::modp1,
                                   m::prime,
                                   $)
    local i;
    local dg := modp1(Degree(p), m);                                   
    return  seq(modp1(Coeff(p, i), m), i in 0..dg);                                 
end proc;

# We compute all the subresultants of index 0 and 1 between
# f and rc.
local SubresultantOfIndexZeroAndOne::static := proc(f::polynom, 
                                                    rc::TRDrc,
                                                    deg_s::posint, 
                                                    ini_in::polynom,
                                                    R::TRDring, 
                                                    $)
    local src1, src2, deg_r1;
    local m := R[characteristic];

    local u := R[variables][2];
    local v := R[variables][1];
    local p := Polynomial(v, rc, R);
    local h := TRDinit(p, R);

    # We make p monic.
    local q := p/h mod m;

    # We compute the subresultant chain between q and f.
    if MainDegree(f, R) <= MainDegree(q, R) then
        src1 := ModularSRCForBivariatePolynomials(q, f mod m, v, u, m,
                                         ':-mode' = modp1);
    else 
        src1 := ModularSRCForBivariatePolynomials(f mod m, q, v, u, m,
                                         ':-mode' = modp1);
    end if;

    # We get the resultant.
    local r1 := src1[1];
    # We get the subresultant of index 1.
    local g1 := src1[2];

    p := Polynomial(u, rc, R);
    h := TRDinit(p, R);

    # We make p monic.
    q := p/h mod m;

    deg_r1 := modp1(Degree(r1), m);

    # We check our assumptions. 
    # r1 must be a polynomial in y.
    if deg_r1 < 1 then
        error "case not covered, %1 is constant", r1;
    end if;     

    # We compute the subresultant chain between r1 and q.
    if MainDegree(q, R) <= deg_r1 then
        src2 := subresultant_chain(r1, q, u, m,
                                         ':-mode' = modp1);
    else 
        src2 := subresultant_chain(q, r1, u, m,
                                         ':-mode' = modp1);
    end if;
    
    # We get the subresultant of index 1.
    local g2 := subresultant_chain_index_one(src2,m);

    # Resultant between the initial of g1 and g2.
    local adi_comp1 := modp1(Resultant(g1[-1],g2), m);

    # Resultant between the initial of t and g2.
    local ini := modp1(ConvertIn(ini_in, u), m);
    local adi_comp2 := modp1(Resultant(ini,g2), m);

    return r1, g1, src2, adi_comp1, adi_comp2;

end proc:

local SimplifyResultant::static := proc(p::polynom, res::polynom, $)
    local result, dummy;

    Gcd(p, res, 'dummy', 'result') mod m;
    return result;
end proc;

# Univariate squarefree part of p.
local SquarefreePart::static := proc(p::polynom, v::name, $)
    local result, dummy;
    local p_prime := diff(p, v) mod m;
    
    Gcd(p, p_prime, 'result', 'dummy') mod m;

    return Normal(result) mod m;
end proc;

# We interpolate the coefficients of a polynomial.
local InterpolatePolynomials::static := proc(X::list(integer), 
                                             G::Array,
                                             vars:: {list(name), name},
                                             v::name,
                                             in_num_coef::{nonnegint, list(nonnegint)},
                                             M,
                                             m::prime,
                                             {mode :: identical(:-polynom, :-modp1) := ':-polynom'},
                                             $)

    local i, c, p, num_coef, poly;
    if mode = ':-polynom' then
        num_coef := in_num_coef;
    else 
        num_coef := add(in_num_coef) +2;
    end if;
    # We resize all the vectors in G.
    local Y := map2(Vector[':-row'], num_coef, G);
    # We save the columns of Y mod m.
    Y := [seq([seq(v[i] mod m , v in Y)] , i=1..num_coef)];
    # We interpolate the columns of G.
    local interpolants := [seq(modp1(':-Interp'( X,
                                            c, v), m) , c in Y)];

    # We ConvertOut from modp1.
    interpolants := [seq(PolynomialTools:-FromCoefficientList(
                                        modp1(ConvertOut(p), m), v), 
                                        p in interpolants)];
    
    # We apply rational recontruction to each p.
    local rational_reconstruction := [seq( Ratrecon(p, M, v) mod m, 
                                           p in interpolants)];

    #This should not be necessary.
    rational_reconstruction := [seq(ifelse(p<>FAIL, p, 0),
                                     p  in rational_reconstruction)];

    # We compute the cofactors.
    local numerators := map(numer, rational_reconstruction);
    local denominators := map(denom, rational_reconstruction);
    local lcm_denoms := lcm(op(denominators)) mod m;
    local cofactors := map(normal, map2(`/`, lcm_denoms, denominators)) mod m;
  
    # We clean the denominators.
    local cleanning_denominators := map(normal, numerators *~ cofactors) mod m;

    # We reconstruct our polynomial using the order given by vars
    if mode= ':-polynom' then
        poly := PolynomialTools:-FromCoefficientList(cleanning_denominators, vars);
    else 
        poly := [PolynomialTools:-FromCoefficientList(cleanning_denominators[1..(in_num_coef[1]+1)], vars[2]),
                 PolynomialTools:-FromCoefficientList(cleanning_denominators[(in_num_coef[1]+2)..], vars[2])];
        #print(poly);
        poly := PolynomialTools:-FromCoefficientList(poly, vars[1]);
    end if;

    return poly;
end proc;


#GOOD SPECIALIZATION
# We return rc normalized and specialized to ls.
local IsGoodSpecializationPoint::static := proc(ls::list, 
                                                lp::list, 
                                                rc::TRDrc, 
                                                R::TRDring, 
                                                $)
	local p, h, h_specialized, rc_specialized, m,
          chk, inv;
    rc_specialized := TRDempty_regular_chain(R);
    m := R[characteristic];

    # We check the first element of lp.
    h := TRDinit(lp[1], R);
    h_specialized := eval(h, ls) mod m;

    if h_specialized = 0 then 
        return [false]; 
    end if;

    local lp_v := h_specialized^(-1)*eval(lp[1],ls) mod m;
    rc_specialized := Chain([lp_v], rc_specialized, R);
    

   # We check the second element of lp.
    h := TRDinit(lp[2], R);
    h_specialized := eval(h, ls) mod m;

    if h_specialized = 0 then 
        return [false]; 
    end if;

    chk := IsRegularMod(h_specialized, lp_v, R, mode=':-inverse');

    if chk[1] then
        lp_v := Rem(chk[2]*eval(lp[2],ls), lp_v, R[variables][2]) mod m; 
        rc_specialized := Chain([lp_v], rc_specialized, R);
    else
        return [false];
    end if; 

    return [true, rc_specialized];
end proc:


# Check that a poly is regular mod a polynomial rc.
# The inverse mode, returns also the inverse of f mod rc.
local IsRegularMod::static := proc(f::polynom,
                                   rc::polynom, 
                                   R::TRDring,
                                   {mode :: identical(:-normal, :-inverse) := ':-normal'},
                                   $)
    local g, inv, dummy;
    local m := R[characteristic];

    if type(f, constant) and f<>0 then
        if mode=':-normal' then
            return true;
        else 
            return [true, f^(-1)];
        end if;
    end if;

    if type(rc, polynom) then
        if mode=':-normal' then
            g:=Gcd(f, rc) mod m;
        else 
            g:=Algebraic:-ExtendedEuclideanAlgorithm(f, rc, R[variables][2],
                                                    'inv', 'dummy') mod m;
        end if;

        if degree(g) = 0 then
            if mode=':-normal' then
                return true;
            else 
                return [true, inv];
            end if;
        else 
            if mode=':-normal' then
                return false;
            else 
                return [false];
            end if;
        end if;
    end if;

end proc;

local DegOneCase::static := proc(f::polynom, 
                                 rc::TRDrc, 
                                 R::TRDring,
                                 $)

    local m := R[characteristic];
    local termOrder := R[variables]; 
    local t := Polynomial(termOrder[1], rc, R);
    local b := Polynomial(termOrder[2], rc, R);

    # Inital and tail of f.
    local h_f := TRDinit(f, R);
    local t_f := Tail(f, R);

    # Inital and tail of t.
    local h_t := TRDinit(t, R);
    local t_t := Tail(t, R);

    # A polynomial in y and z.
    local r := Normal(h_f*t_t -  h_t*t_f) mod m;

    local src := ModularSRCForBivariatePolynomials(r, b, 
                                          termOrder[2],
                                          termOrder[3],
                                          m); 

    local s := src[1];
    local g := src[2];

    return s, g, undefined;
end proc;

local subresultant_chain;
local subresultant_chain_index_one;
local ModularSRCForBivariatePolynomials;

$include "/src/resultant_chain.mm"
$include "/src/ModularSRCForBivariatePolynomials.mm"
end module;

















