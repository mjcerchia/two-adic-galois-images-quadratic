
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.be.2.
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve by intersecting with a hyperplane. From this, we are able to construct an elliptic curve, 
which we find to have rank 1. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[-x*t + y^2, x*t + y^2 - z*w, 16*x^2 - 4*z^2 - w^2 + 2*t^2]);

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
C1 := Curve(P,[162*x[1]^2 - 2592*x[2]*x[3] - 5184*x[3]^2 - 576*x[2]*x[6] - 2336*x[6]^2 - 
    972*x[5]*x[7] - 27*x[7]^2 + 828*x[5]*x[8] - 990*x[7]*x[8],
162*x[1]*x[2] + 1944*x[2]*x[3] + 2592*x[3]^2 - 216*x[2]*x[6] - 128*x[6]^2 + 
    675*x[5]*x[7] + 81*x[7]^2 - 900*x[5]*x[8] + 216*x[7]*x[8] + 408*x[8]^2,
648*x[1]*x[3] + 648*x[2]*x[3] + 2592*x[3]^2 + 144*x[2]*x[6] + 448*x[6]^2 + 
    27*x[5]*x[7] - 81*x[7]^2 - 252*x[5]*x[8] + 801*x[7]*x[8] + 672*x[8]^2,
12*x[1]*x[4] + 32*x[6]^2 - 3*x[5]*x[7] + 2*x[5]*x[8] + 12*x[8]^2,
3*x[1]*x[5] + 48*x[3]*x[5] - 64*x[4]*x[5] - 44*x[5]*x[6] + 12*x[3]*x[7] - 
    10*x[4]*x[7] + 8*x[6]*x[7],
4*x[1]*x[6] - 16*x[6]^2 - x[5]*x[8],
-24*x[3]*x[5] + 32*x[4]*x[5] + 16*x[5]*x[6] + 3*x[1]*x[7] - 12*x[6]*x[7],
3*x[4]*x[5] + 2*x[5]*x[6] + 3*x[1]*x[8] - 12*x[6]*x[8],
162*x[2]^2 - 1296*x[2]*x[3] - 2592*x[3]^2 - 576*x[2]*x[6] - 896*x[6]^2 - 
    432*x[5]*x[7] + 27*x[7]^2 + 576*x[5]*x[8] - 594*x[7]*x[8] - 288*x[8]^2,
72*x[2]*x[4] + 48*x[2]*x[6] + 18*x[5]*x[7] + 9*x[7]^2 - 24*x[5]*x[8] - 
    24*x[7]*x[8] + 16*x[8]^2,
9*x[2]*x[5] - 108*x[3]*x[5] + 132*x[4]*x[5] + 56*x[5]*x[6] - 36*x[3]*x[7] + 
    24*x[4]*x[7] - 4*x[6]*x[7] + 48*x[6]*x[8],
72*x[3]*x[5] - 84*x[4]*x[5] - 40*x[5]*x[6] + 9*x[2]*x[7] + 36*x[3]*x[7] - 
    36*x[4]*x[7] - 24*x[6]*x[7],
-18*x[4]*x[5] + 12*x[5]*x[6] - 9*x[4]*x[7] + 18*x[6]*x[7] + 18*x[2]*x[8] - 
    16*x[6]*x[8],
864*x[3]*x[4] - 128*x[6]^2 - 27*x[7]^2 + 306*x[7]*x[8] - 240*x[8]^2,
288*x[3]*x[6] + 64*x[6]^2 - 9*x[7]*x[8] + 96*x[8]^2,
9*x[4]*x[7] - 90*x[6]*x[7] + 72*x[3]*x[8] + 16*x[6]*x[8],
36*x[4]^2 - 16*x[6]^2 + 9*x[7]*x[8] - 12*x[8]^2,
12*x[4]*x[6] + 8*x[6]^2 + 3*x[8]^2,
-3*x[6]*x[7] + 3*x[4]*x[8] + 2*x[6]*x[8],
2*x[5]^2 + 4*x[5]*x[7] + x[7]^2 - 4*x[7]*x[8] - 8*x[8]^2]);

//We can't immediately find a point, so we intersect with a hyperplane
pt := C1!Points(C1 meet Scheme(AmbientSpace(C1),x[5]))[1];

E := EllipticCurve(C1,pt);
Rank(E); //1
