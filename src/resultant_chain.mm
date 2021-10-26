##########
# Ducos' Subresutlant Chain Algorithm without Ducos' optimization! 
# given polynomials p_in, q_in
# reutrns the list of subresutlant chain w.r.t variable v
##########
## This needs to be made a module.
subresultant_chain := proc(p_in ::depends({modp1, polynom(integer,[v])}),
                           q_in ::depends({modp1, polynom(integer,[v])}), 
                           v :: symbol,
                           m::prime, 
                           {mode :: identical(:-polynom, :-modp1) := ':-polynom'},
                           $)
    local p, q, delta, dp, dq, s, t, src, z,
          p1, p2, quoCoef, i;

    # Minus one.
    local nOne := modp1(Constant(-1, v), m);
    # Zero.
    local Zero := modp1(Constant(0, v), m);

    #Create the "C objects" in case it is necessary
    if p_in :: ':-polynom'( ':-integer', [v]) then
        p := modp1(ConvertIn(p_in, v), m);
    else
        p := p_in;
    end if;

    if q_in :: ':-polynom'( ':-integer', [v] ) then
        q := modp1(ConvertIn(q_in, v), m);
    else
        q := q_in;
    end if;

    # Degree of p and q.
    dp := modp1(Degree(p), m);
    dq := modp1(Degree(q), m);

    if dp < 1 or dq < 1 then 
        return [];
    elif dp < dq then 
        p, q := q, p;
        # Degree of p and q.
        dp := modp1(Degree(p), m);
        dq := modp1(Degree(q), m);
    end if;

    if mode = ':-polynom' then
        src := [modp1(ConvertOut(q, v), m), 
                modp1(ConvertOut(p, v), m)];
    else
        src := [q, p];
    end if;

    delta := dp - dq;

    if delta <> 0 then
        s :=  modp1(Lcoeff(q), m)^(delta-1);
        s := modp1(Constant(s, v), m);
        if mode = ':-polynom' then
            src := [modp1(ConvertOut(Multiply(s, q), v), m), 
                    op(src)];
        else
            src := [modp1(Multiply(s, q) , m), op(src)];
        end if;
    end if;

    s :=  modp1(Lcoeff(q), m)^(delta);
    s := modp1(Constant(s, v), m);

    # First subresultant: 
    t := q;
    q := modp1(Prem(p, Multiply(nOne,q)), m);
    if mode = ':-polynom' then
        src := [modp1(ConvertOut(q, v), m), op(src)];
    else
        src := [q, op(src)];
    end if;
    p := t;

    # Degree of p.
    dp := modp1(Degree(p), m);

    # To compute other subresutlants in the chain:
    do 
        if q = Zero then 
            src := [seq(0, 1..dp-1), op(src)]; 
            return src;
        end if;

        # Degree of q.
        dq := modp1(Degree(q), m);
        # Difference between degrees.
        delta := dp - dq;
        
        # Coefficient of q power dl-1.
        p1 :=  modp1(Lcoeff(q), m)^(delta-1);
        p1 := modp1(Constant(p1, v), m);

        # Coefficient of p power dl-1.
        #p2 := modp1(Lcoeff(p), m)^(delta-1);
        #p2 := modp1(Constant(p2, v), m);

        # s power delta. 
        p2 := modp1(Power(s,delta-1), m);

        # Quotient between p1 and p2
        quoCoef := modp1(Quo(p1, p2), m);

        z := modp1(Multiply(quoCoef, q) , m);
        ### 
        if delta >= 2 then 
            if mode = ':-polynom' then 
                src := [modp1(ConvertOut(z , v), m), 
                        seq(0, i=1..delta-2), op(src)];
            else
                src := [z, seq(0, i=1..delta-2), op(src)];
            end if;
        end if;

        if modp1(Degree(z), m) = 0 then
            return src;
        end if;

        p1 := modp1(Prem(p, Multiply(nOne,q)), m);
        p2 := modp1(Multiply(Power(s,delta), 
                             Constant(Lcoeff(p), v)), m);
        q := modp1(Quo(p1, p2), m);

        ####
        if mode = ':-polynom' then 
            src := [modp1(ConvertOut(q , v), m),
                    op(src)];
        else
            src := [q, op(src)];
        end if;
        p := z;
        s := modp1(Lcoeff(p), m);
        s := modp1(Constant(s, v), m);

        # Degree of p.
        dp := modp1(Degree(p), m);
    end do;
end proc:

# Returns the subresultant of index 1, whenever it exists.
subresultant_chain_index_one := proc(src :: {list(polynom), list(modp1)}, 
                                     m::prime,
                                    {mode :: identical(:-polynom, :-modp1) := ':-polynom'},
                                     $)
    local dp, dq;
    local num_elem := numelems(src);

    if mode = ':-polynom' then
        dp := degree(src[num_elem]);
        dp := degree(src[num_elem -1 ]);
    else
        dp := modp1(Degree(src[num_elem]), m);
        dq := modp1(Degree(src[num_elem-1]), m);
    end if;

    if dp = 1 and dq = 1 then
        return undefined;
    else
        return src[2];
    end if;
end proc;
