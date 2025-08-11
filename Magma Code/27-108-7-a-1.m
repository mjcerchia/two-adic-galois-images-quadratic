/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 27.108.7.a.1.
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve by intersecting with a hyperplane. From this, we are able to construct an elliptic curve, 
which we find to have rank 1. 
******************************************************************************/
P<x,y,z,u,v,w,t> := ProjectiveSpace(Rationals(),6);
C := Curve(P,[x*y - w*t, x*t - x*u + x*v + y*z, x*z + w^2 + w*u + w*v + t*u + u*v, x*y + x*z - 2*w^2 - w*u - t*u + u^2, 2*x*u + x*v - z^2, 2*y*u + y*v + z*t - z*u + z*v, 3*x^2 - y*w - y*u + z*w + z*u, y*t - y*u + y*v + 3*z*w + z*t, 3*x*w + x*t - y^2, 2*x*z + 2*w^2 - 4*w*u - w*v + t^2 - 2*t*u + t*v + u^2 + v^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//We compute genus one quotients by an involution.
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

//One of these quotients has the following model:
P<[x]> := ProjectiveSpace(Rationals(),11);
C1 := Curve(P,[x[1]^2 - x[2]*x[4] - 9*x[10]*x[11] - 81*x[12]^2,
x[1]*x[2] - x[4]^2,
x[1]*x[3] - x[4]*x[7] - 9*x[6]*x[8],
-x[2]^2 + x[1]*x[4] - 9*x[7]*x[8],
x[1]*x[5] - 9*x[8]^2 - x[4]*x[10],
x[1]*x[6] - x[7]^2 - 9*x[8]*x[9],
-x[4]*x[5] + x[1]*x[7],
x[1]*x[8] - x[7]*x[10] - 9*x[9]*x[11],
x[1]*x[9] - x[8]*x[10] - 9*x[9]*x[12],
-x[4]*x[7] + x[1]*x[10],
-x[10]^2 + x[1]*x[11] - 9*x[11]*x[12],
-x[10]*x[11] + x[1]*x[12] - 9*x[12]^2,
x[2]*x[3] - x[4]*x[5],
x[2]*x[5] - x[4]*x[7],
x[2]*x[6] - x[10]^2 - 9*x[11]*x[12],
x[2]*x[7] - x[4]*x[10],
-x[7]^2 + x[2]*x[8],
-x[7]*x[8] + x[2]*x[9],
-x[4]*x[5] + x[2]*x[10] + 9*x[9]*x[10],
-x[7]*x[10] + x[2]*x[11],
-x[8]*x[10] + x[2]*x[12],
x[3]^2 - 9*x[6]*x[9] - x[10]^2 - 9*x[11]*x[12],
x[3]*x[4] - 9*x[8]^2 - x[4]*x[10],
x[3]*x[5] - x[7]^2 - 9*x[8]*x[9],
x[3]*x[6] - x[7]*x[8] - 9*x[9]^2,
x[3]*x[7] - x[7]*x[10] - 9*x[9]*x[11],
x[3]*x[8] - x[8]*x[10] - 9*x[9]*x[12],
-x[6]^2 + x[3]*x[9],
x[3]*x[10] - x[10]^2 - 9*x[11]*x[12],
x[3]*x[11] - x[10]*x[11] - 9*x[12]^2,
-x[6]*x[8] + x[3]*x[12],
x[4]*x[6] - x[7]*x[10] - 9*x[9]*x[11],
x[4]*x[8] - x[10]^2 - 9*x[11]*x[12],
x[4]*x[9] - x[10]*x[11] - 9*x[12]^2,
-x[7]^2 + x[4]*x[11],
-x[7]*x[8] + x[4]*x[12],
x[5]^2 - x[7]*x[10] - 9*x[9]*x[11],
x[5]*x[6] - x[8]*x[10] - 9*x[9]*x[12],
x[5]*x[7] - x[10]^2 - 9*x[11]*x[12],
x[5]*x[8] - x[10]*x[11] - 9*x[12]^2,
-x[6]*x[8] + x[5]*x[9],
-x[7]^2 + x[5]*x[10],
-x[7]*x[8] + x[5]*x[11],
-x[8]^2 + x[5]*x[12],
x[6]*x[7] - x[10]*x[11] - 9*x[12]^2,
-x[7]*x[8] + x[6]*x[10],
-x[8]^2 + x[6]*x[11],
-x[8]*x[9] + x[6]*x[12],
-x[8]^2 + x[7]*x[9],
-x[8]*x[10] + x[7]*x[11],
-x[9]*x[10] + x[7]*x[12],
-x[9]*x[10] + x[8]*x[11],
-x[9]*x[11] + x[8]*x[12],
-x[11]^2 + x[10]*x[12]]);

//We can't immediately find a point, so we intersect with hyperplanes.
pt := C1!Points(C1 meet Scheme(AmbientSpace(C1),x[3]))[1];
E := EllipticCurve(C1,pt);
Rank(E); //1
