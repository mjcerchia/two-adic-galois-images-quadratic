32.96.5.n.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.n.2.
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve by intersecting with a hyperplane. From this, we are able to construct an elliptic curve, 
which we find to be rank 1. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y^2 - z*w, 2*x^2 - y*z - 2*y*t, z^2 + 2*z*t + 8*w^2 + 2*t^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are three genus one quotients by an involution
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
P<[x]> := ProjectiveSpace(Rationals(),7);
C1 := Curve(P,[2*x[1]^2 + 2*x[1]*x[3] + x[3]^2 - 4*x[5]*x[8] - 2*x[7]*x[8],
4*x[1]*x[2] + 2*x[5]*x[8] + x[7]*x[8],
8*x[1]*x[4] - 2*x[5]*x[7] - 5*x[7]^2 + 32*x[8]^2,
8*x[1]*x[5] + 2*x[3]*x[5] - 5*x[3]*x[7] + 48*x[4]*x[8] + 24*x[6]*x[8],
4*x[1]*x[6] + 2*x[5]*x[7] + x[7]^2,
2*x[3]*x[5] + 4*x[1]*x[7] + 3*x[3]*x[7] - 16*x[4]*x[8] - 8*x[6]*x[8],
4*x[4]*x[5] + 4*x[5]*x[6] + x[6]*x[7] + 4*x[1]*x[8],
32*x[2]^2 - 4*x[6]^2 - 2*x[5]*x[8] - x[7]*x[8],
4*x[2]*x[3] - 2*x[5]*x[8] - 3*x[7]*x[8],
2*x[2]*x[4] + x[2]*x[6] - x[8]^2,
x[2]*x[5] - x[4]*x[8],
x[2]*x[7] - x[6]*x[8],
-4*x[4]*x[5] - 4*x[5]*x[6] - 5*x[6]*x[7] + 32*x[2]*x[8],
8*x[3]*x[4] - 2*x[5]*x[7] + 5*x[7]^2 - 32*x[8]^2,
4*x[3]*x[6] - 2*x[5]*x[7] - 3*x[7]^2,
-4*x[4]*x[5] - 8*x[5]*x[6] - 3*x[6]*x[7] + 4*x[3]*x[8],
4*x[4]^2 - x[6]^2 - 2*x[5]*x[8] + x[7]*x[8],
2*x[4]*x[6] + x[6]^2 - x[7]*x[8],
-x[5]*x[6] + x[4]*x[7],
4*x[5]^2 + 4*x[5]*x[7] + 5*x[7]^2 - 32*x[8]^2]);

//We can't immediately find a point, so we intersect with a hyperplane
Points(C1 meet Scheme(AmbientSpace(C1),x[3]));

//This produces a point
pt := C1![-2,1/2,0,3/2,3,-1,-2,1]; 
E := EllipticCurve(C1,pt);
Rank(E); //1
