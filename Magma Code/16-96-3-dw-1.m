/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.h.1. 
We find there are five genus one quotients by an involution. 
One of them is a rank 1 elliptic curve
******************************************************************************/
P<x,y,z,u,t,w>:=ProjectiveSpace(Rationals(),5);
 C := Curve(P,[x^2 - y*z, 2*y*w + 2*w^2 - t*u + u^2, 2*y*w - 2*w^2 - t^2 - t*\
u, y*w + 4*z*w - t*u, -y*u + 4*z*t - w*t + w*u, y*t + 4*z*u + w*t + w*u, y^2 +\
 16*z^2 - 2*w^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are five genus one quotients by an involution
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

l;
/*
[
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    -16*x[1]^2 - 32*x[1]*x[2] - 40*x[2]^2 - 24*x[2]*x[3] - 7*x[3]^2 + x[4]^2,
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    16*x[1]^2 + 32*x[2]^2 + 8*x[3]^2 + x[4]^2,
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    8*x[1]^2 + 32*x[2]^2 + 16*x[3]^2 + x[4]^2,
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    16*x[1]^2 - 96*x[2]^2 + 16*x[3]^2 + x[4]^2,
    Curve over Rational Field defined by
    x[2]^2 - x[1]*x[3],
    -16*x[1]^2 - 32*x[1]*x[2] - 40*x[2]^2 - 24*x[2]*x[3] - 7*x[3]^2 + x[4]^2
]
*/

//One of these quotients has the following model:
P<[x]>:=ProjectiveSpace(Rationals(),3);
C := Curve(P,[x[2]^2 - x[1]*x[3], 16*x[1]^2 + 32*x[2]^2 + 8*x[3]^2 + x[4]^2]);
Rank(EllipticCurve(H)); //1
