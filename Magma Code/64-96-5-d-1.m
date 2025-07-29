
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 64.96.5.d.1.
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve by intersecting with a hyperplane. From this, we are able to construct an elliptic curve, 
which we find to have rank 1. 
******************************************************************************/
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[x*z - y*w, -4*x*w + y*z + t^2, 4*x^2 + y^2 + z*w]);

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
C1 := Curve(P,[384*x[1]^2 + 4608*x[1]*x[5] + 15360*x[5]^2 + 5376*x[5]*x[6] - 96*x[6]^2 + 
    96*x[4]*x[7] + 1600*x[7]^2 - 54*x[4]*x[8] + 24*x[7]*x[8] - 3*x[8]^2,
12*x[1]*x[2] + 384*x[1]*x[5] - 288*x[5]*x[6] - 72*x[6]^2 + 24*x[4]*x[7] + 
    304*x[7]^2 - 3*x[4]*x[8] + 22*x[7]*x[8],
4*x[1]*x[3] + 16*x[6]^2 - 4*x[7]*x[8] - x[8]^2,
144*x[1]*x[4] - 33*x[3]*x[4] + 960*x[4]*x[5] + 264*x[4]*x[6] + 640*x[5]*x[7] + 
    704*x[6]*x[7] - 288*x[5]*x[8] - 76*x[6]*x[8],
24*x[1]*x[6] + 48*x[5]*x[6] + 12*x[6]^2 - 6*x[4]*x[7] - 40*x[7]^2 - 3*x[7]*x[8],
3*x[3]*x[4] - 192*x[4]*x[5] - 96*x[4]*x[6] + 192*x[1]*x[7] + 1024*x[5]*x[7] + 
    320*x[6]*x[7] - 16*x[6]*x[8],
-3*x[3]*x[4] - 12*x[4]*x[6] + 64*x[6]*x[7] + 12*x[1]*x[8] + 4*x[6]*x[8],
3*x[2]^2 - 3072*x[5]^2 - 2304*x[5]*x[6] - 432*x[6]^2 + 192*x[7]^2 - 3*x[4]*x[8] 
    + 28*x[7]*x[8] + 3*x[8]^2,
3*x[2]*x[3] + 384*x[5]*x[6] + 288*x[6]^2 + 48*x[4]*x[7] + 128*x[7]^2 + 
    12*x[4]*x[8] + 8*x[7]*x[8] - 9*x[8]^2,
18*x[2]*x[4] - 9*x[3]*x[4] + 576*x[4]*x[5] + 264*x[4]*x[6] + 224*x[6]*x[7] - 
    144*x[5]*x[8] - 24*x[6]*x[8],
64*x[2]*x[5] + 2048*x[5]^2 + 768*x[5]*x[6] - 64*x[7]^2 + 16*x[7]*x[8] - x[8]^2,
x[2]*x[6] + 32*x[5]*x[6] + 12*x[6]^2 - x[7]*x[8],
-3*x[4]*x[6] + 3*x[2]*x[7] + 96*x[5]*x[7] + 40*x[6]*x[7],
-3*x[3]*x[4] + 3*x[2]*x[8] + 96*x[5]*x[8] + 40*x[6]*x[8],
3*x[3]^2 - 192*x[6]^2 - 24*x[4]*x[8] - 64*x[7]*x[8] + 12*x[8]^2,
48*x[3]*x[5] + 384*x[5]*x[6] + 48*x[4]*x[7] + 128*x[7]^2 - 6*x[4]*x[8] - 
    40*x[7]*x[8] + 3*x[8]^2,
3*x[3]*x[6] - 96*x[5]*x[6] - 24*x[6]^2 - 12*x[4]*x[7] - 32*x[7]^2 + 6*x[7]*x[8],
x[3]*x[7] - x[6]*x[8],
64*x[6]*x[7] + x[3]*x[8] - 64*x[5]*x[8] - 16*x[6]*x[8],
72*x[4]^2 + 96*x[4]*x[7] + 320*x[7]^2 - 36*x[4]*x[8] + 48*x[7]*x[8] - 9*x[8]^2]);

//We can't immediately find a point, so we intersect with a hyperplane
pt := C1!Points(C1 meet Scheme(AmbientSpace(C1),x[4]))[1];

E := EllipticCurve(C1,pt);
Rank(E); //1
