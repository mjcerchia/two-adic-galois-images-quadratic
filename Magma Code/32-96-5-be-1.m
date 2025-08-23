/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.be.1. 
We find there are three genus one quotients by an involution. 
We are able to find a point on one of these. From this, we are able to construct an elliptic curve, 
which we find to be rank 1. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*w + z*t, 2*x^2 + y*w + y*t - z*t, y^2 + 4*y*z - 4*z^2 + w^2 + w*t]);
pt:=PointSearch(C,1000)[1];
S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
assert #auts eq #S;




for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
CG,prj := CurveQuotient(AG);
if Genus(CG) eq 1 then
E := EllipticCurve(CG,CG!prj(pt)); 
Rank(E); // one of these has rank 1
end if;

end if;

end for;


