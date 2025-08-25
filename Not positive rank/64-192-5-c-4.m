/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 64.192.5.c.4.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. It has rank 0.
NOT bielliptic. 
******************************************************************************/
P<[x]> := ProjectiveSpace(Rationals(),4);
/*
Model obtained from David Zywina's Github:
*/

C := Curve(P,[-x[1]*x[5] + x[4]^2,
x[1]^2 + x[2]*x[3],
-4*x[2]^2 + x[3]^2 + x[5]^2]);

S := AutomorphismGroup(C); 
Pt:=PointSearch(C,100)[1];
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
assert #auts eq #S;


//There is one genus one quotient by an involution

for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
CG,prj := CurveQuotient(AG);
if Genus(CG) eq 1 then
E:=EllipticCurve(CG,prj(Pt));
Rank(E);

end if;

end if;

end for;

