/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.bo.1.
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,r,s,v>:=ProjectiveSpace(Rationals(),8);

C := Curve(P,[x^2 + u*v, x^2 - t*r, y*u - w*r, y*t + w*v, y*r - z*u, y*v + z\
*t, y^2 - z*w, z^2 - v^2 + r^2, y*z + t*v + u*r, y^2 - t^2 + u^2, t^2 + u^2 - \
v*r + s^2, y^2 + z*w - z*s - v*r + r^2, y^2 + z*w + z*s + v^2 + v*r, 2*w^2 - w\
*s + t*u + u^2, 2*w^2 + w*s + t^2 - t*u, y*t + y*u - z*v - w*v + w*r - v*s + r\
*s, y*t - y*u + z*r - w*v - w*r + v*s + r*s, x^2 + 2*y*w - y*s + u*r, x^2 - y*\
z - 2*y*s + t*r, y*r + 2*w*t - 2*w*u - t*s + u*s, y*v - 2*w*t - 2*w*u - t*s - \
u*s]);
 
C3 := Curve(Reduction(C,3)); 
S := AutomorphismGroup(C3); 

auts := [];
Stemp := Automorphisms(C3);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotient by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C3,[g]);
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
    Curve over GF(3) defined by
    x[3]*x[5] + 2*x[10]*x[12] + x[1]*x[16],
    x[1]*x[2] + x[3]*x[6] + 2*x[12]^2,
    x[3]*x[7] + 2*x[1]*x[12] + 2*x[10]*x[12],
    x[3]*x[8] + 2*x[1]*x[13] + 2*x[12]*x[14] + x[1]*x[15],
    x[3]*x[9] + 2*x[1]*x[14] + 2*x[12]*x[15] + x[15]*x[16],
    x[1]*x[6] + x[2]*x[6] + x[3]*x[10] + 2*x[7]*x[12],
    x[1]*x[4] + x[1]*x[9] + 2*x[9]*x[10] + x[3]*x[11],
    x[1]*x[7] + 2*x[7]*x[10] + x[3]*x[12],
    x[1]*x[8] + x[2]*x[8] + 2*x[9]*x[12] + x[3]*x[13],
    x[1]*x[9] + 2*x[9]*x[10] + x[3]*x[14],
    x[2]*x[8] + 2*x[9]*x[12] + x[3]*x[15],
    2*x[1]*x[5] + x[3]*x[16],
    2*x[1]^2 + x[4]^2 + 2*x[1]*x[10],
    x[4]*x[5] + 2*x[1]*x[13] + 2*x[12]*x[14] + x[1]*x[15],
    x[4]*x[6] + 2*x[1]*x[14] + 2*x[12]*x[15] + x[15]*x[16],
    x[4]*x[7] + 2*x[12]*x[14] + 2*x[1]*x[15],
    x[4]*x[8] + 2*x[1]*x[12] + 2*x[1]*x[16],
    x[1]*x[2] + x[4]*x[9] + 2*x[1]*x[10],
    x[1]*x[9] + x[4]*x[10] + 2*x[9]*x[10],
    2*x[1]*x[3] + x[4]*x[11],
    x[2]*x[8] + x[4]*x[12] + 2*x[9]*x[12],
    2*x[1]*x[5] + x[4]*x[13],
    2*x[1]*x[6] + x[4]*x[14],
    2*x[1]*x[7] + x[4]*x[15],
    2*x[1]*x[8] + x[4]*x[16],
    x[1]*x[2] + x[5]^2 + 2*x[12]^2,
    x[5]*x[6] + 2*x[1]*x[12] + 2*x[10]*x[12],
    2*x[1]*x[2] + x[5]*x[7] + 2*x[1]*x[10] + 2*x[12]^2,
    x[5]*x[8] + 2*x[1]*x[14] + 2*x[12]*x[15] + x[15]*x[16],
    x[5]*x[9] + 2*x[12]*x[14] + 2*x[1]*x[15],
    x[1]*x[7] + x[5]*x[10] + 2*x[7]*x[10],
    x[1]*x[8] + x[2]*x[8] + x[5]*x[11] + 2*x[9]*x[12],
    x[2]*x[6] + x[5]*x[12] + 2*x[7]*x[12],
    x[1]*x[9] + 2*x[9]*x[10] + x[5]*x[13],
    x[2]*x[8] + 2*x[9]*x[12] + x[5]*x[14],
    2*x[9]*x[10] + x[5]*x[15],
    2*x[1]*x[6] + x[5]*x[16],
    2*x[1]*x[2] + x[6]^2 + 2*x[1]*x[10] + 2*x[12]^2,
    x[6]*x[7] + 2*x[10]*x[12] + x[14]*x[15],
    x[6]*x[8] + 2*x[12]*x[14] + 2*x[1]*x[15],
    x[6]*x[9] + 2*x[12]*x[15] + 2*x[15]*x[16],
    x[2]*x[6] + x[6]*x[10] + 2*x[7]*x[12],
    x[1]*x[9] + 2*x[9]*x[10] + x[6]*x[11],
    2*x[7]*x[10] + x[6]*x[12],
    x[2]*x[8] + 2*x[9]*x[12] + x[6]*x[13],
    2*x[9]*x[10] + x[6]*x[14],
    2*x[9]*x[12] + x[6]*x[15],
    2*x[1]*x[7] + x[6]*x[16],
    x[1]*x[2] + x[2]^2 + x[7]^2 + x[1]*x[10] + 2*x[12]^2,
    x[7]*x[8] + 2*x[12]*x[15] + 2*x[15]*x[16],
    x[7]*x[9] + 2*x[12]*x[14] + x[2]*x[15],
    x[3]*x[4] + 2*x[1]*x[11] + x[1]*x[14] + 2*x[12]*x[15] + x[15]*x[16],
    x[2]*x[8] + x[7]*x[11] + 2*x[9]*x[12],
    x[1]^2 + x[1]*x[2] + x[3]^2 + x[1]*x[10] + 2*x[12]^2,
    2*x[9]*x[10] + x[7]*x[13],
    2*x[9]*x[12] + x[7]*x[14],
    2*x[2]*x[9] + 2*x[9]*x[10] + x[7]*x[15],
    2*x[2]*x[6] + x[7]*x[16],
    x[1]*x[2] + x[8]^2 + 2*x[1]*x[10],
    x[8]*x[9] + x[1]*x[12] + x[14]*x[15],
    x[2]*x[8] + x[8]*x[10] + 2*x[9]*x[12],
    2*x[1]*x[5] + x[8]*x[11],
    2*x[9]*x[10] + x[8]*x[12],
    2*x[1]*x[6] + x[8]*x[13],
    2*x[1]*x[7] + x[8]*x[14],
    2*x[2]*x[6] + x[8]*x[15],
    2*x[1]*x[9] + x[8]*x[16],
    2*x[1]*x[2] + x[2]^2 + x[9]^2 + 2*x[1]*x[10],
    x[1]*x[12] + 2*x[14]*x[15] + x[2]*x[16],
    2*x[1]*x[6] + x[9]*x[11],
    x[2]*x[14] + 2*x[15]*x[16],
    2*x[1]*x[7] + x[9]*x[13],
    2*x[2]*x[6] + x[9]*x[14],
    2*x[2]*x[7] + x[9]*x[15],
    2*x[2]*x[8] + x[9]*x[16],
    x[1]*x[2] + x[1]*x[10] + x[10]^2 + 2*x[12]^2,
    x[10]*x[11] + x[1]*x[14] + 2*x[12]*x[15] + x[15]*x[16],
    x[2]*x[13] + 2*x[1]*x[15],
    x[10]*x[13] + 2*x[12]*x[14] + x[1]*x[15],
    x[10]*x[14] + 2*x[12]*x[15] + x[15]*x[16],
    2*x[12]*x[14] + x[10]*x[15],
    2*x[1]*x[12] + x[10]*x[16],
    x[1]^2 + 2*x[1]*x[10] + x[11]^2,
    x[11]*x[12] + 2*x[12]*x[14] + x[1]*x[15],
    2*x[1]*x[12] + x[11]*x[13] + x[1]*x[16],
    2*x[1]*x[10] + x[11]*x[14],
    2*x[1]*x[12] + x[11]*x[15],
    2*x[1]*x[13] + x[11]*x[16],
    x[2]*x[12] + 2*x[14]*x[15],
    x[12]*x[13] + 2*x[12]*x[15] + x[15]*x[16],
    x[2]*x[11] + 2*x[1]*x[14],
    2*x[1]*x[2] + 2*x[1]*x[10] + x[2]*x[10],
    2*x[1]*x[2] + 2*x[1]*x[10] + x[12]*x[16],
    2*x[1]*x[10] + x[13]^2,
    2*x[1]*x[12] + x[13]*x[14],
    2*x[1]*x[2] + 2*x[1]*x[10] + x[13]*x[15],
    2*x[1]*x[14] + x[13]*x[16],
    2*x[1]*x[2] + 2*x[1]*x[10] + x[14]^2,
    x[2]*x[5] + 2*x[1]*x[7],
    2*x[1]*x[15] + x[14]*x[16],
    2*x[1]*x[2] + 2*x[2]^2 + 2*x[1]*x[10] + x[15]^2,
    x[2]*x[4] + 2*x[1]*x[9],
    2*x[1]*x[2] + x[16]^2,
    x[2]*x[3] + 2*x[1]*x[6]
]
*/
#EllipticCurve(Curve(l[1])); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
