/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.5.c.3.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. It is an elliptic curve of rank 0.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[x^2 - y*z, y^2 - w^2 + t^2, 2*z^2 + w*t]);

Pt:=PointSearch(C,100)[1];
S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
assert #auts eq #S;


//There is one genus one quotient by an involution which is an elliptic curve of rank 0

for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
CG,prj := CurveQuotient(AG);
if Genus(CG) eq 1 then
E:=EllipticCurve(CG,prj(Pt));
Rank(E); // 0 true

end if;

end if;

end for;

