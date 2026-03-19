//Magma Code to (hopefully) compute the Zeta Function for the Elliptic Surface y^2 = x^3 + t^5 + t^-5
time for j in [1..11] do;

Gp:=GF(5^j);
Gt<t>:=FunctionField(Gp);
E:=EllipticCurve([2*(t^8+14*t^4+1),4*t^2*(t^8+6*t^4+1)]); //From Kloosterman's Paper, has MW rank 15 over CC

pts:=0;

disTop:=Numerator(Discriminant(E));
disBot:=Denominator(Discriminant(E));

Bads:=[Roots(f[1])[1][1] : f in Factorization(Numerator(Discriminant(E)))| Degree(f[1]) eq 1];

Goods:=[i : i in Gp | i notin Bads]; 
Exclude(~Goods,0);
time for i in Goods do;
pts+:=#Specialization(E,i);
end for;

k<x,y>:=AffineSpace(Gp,2);
for kk in Bads do;
f:=y^2-x^3-(kk^5+kk^-5);
S:=Scheme(k,f);
pts+:=#Points(S)+1;
end for;

[j,pts];
end for;
