/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.f.2. 
We find there are three genus one quotients by an involution. 
We are able to find a point on the third curve by intersecting with a hyperplane. From this,
we are able to construct an elliptic curve, which
we find to have rank 1. 
******************************************************************************/

//Model obtained from the LMFDB
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[x*t + y^2, 2*x^2 - z*w, z^2 + 16*w^2 + 2*t^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//We find three genus one quotients by involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
CG,prj := CurveQuotient(AG);prj;
if Genus(CG) eq 1 then
try
l := Append(l,CG);
catch e
m := Append(m,CG);
end try;
end if;
CG; Genus(CG);
end if;
print ".........";
end for;

P<[x]> := ProjectiveSpace(Rationals(),7);
C1 :=Curve(P,[16*x[1]^2 - x[3]^2 - 2*x[7]^2,
x[1]*x[2] - x[7]^2,
x[1]*x[4] - x[3]*x[6],
x[1]*x[5] - x[4]*x[6],
-x[3]*x[4] + 16*x[1]*x[6] - 2*x[7]*x[8],
-x[6]^2 + x[1]*x[7],
-x[6]*x[7] + x[1]*x[8],
2*x[2]^2 + x[5]^2 - 16*x[7]^2,
x[2]*x[3] - x[5]*x[7],
x[2]*x[4] - x[5]*x[8],
x[2]*x[6] - x[7]*x[8],
x[2]*x[7] - x[8]^2,
x[4]*x[5] - 16*x[6]*x[7] + 2*x[2]*x[8],
x[3]*x[5] - 16*x[6]^2 + 2*x[8]^2,
-x[4]*x[6] + x[3]*x[7],
-x[5]*x[6] + x[3]*x[8],
x[4]^2 - 16*x[6]^2 + 2*x[8]^2,
-x[5]*x[6] + x[4]*x[7],
-x[5]*x[7] + x[4]*x[8],
-x[7]^2 + x[6]*x[8]]);

//By intersecting with a hyperplane, we find a point
Points(C1 meet Scheme(AmbientSpace(C1),x[2]));

//This produces a point
pt:=C1![-1/4 , 0 , 1 , 0 , 0 , 0 , 0 , 0];
E := EllipticCurve(C1,pt);
Rank(E); //1
