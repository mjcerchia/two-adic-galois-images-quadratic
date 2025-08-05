32.192.9.i.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.i.2.
We compute the automorphism group over F5 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,u,v,r,w,s,t>:=ProjectiveSpace(Rationals(),8);
C := Curve(P,[y^2 - t*r, y^2 - y*z + t*s, y^2 + y*w + t*v, z*v - w*s, y^2 - \
u*r + v*s - s^2, y*z + u*s - r*s - s^2, y*w + u*v + v*s + r*s, y*w + u*v - v*r\
 - v*s, y^2 - t*s - u*s + v^2 + v*r, t^2 - t*s - u^2 - u*v, y*v + y*r + w*r, y\
*t + y*u - w*v + w*s, y*u - w*t + w*u + w*s, y*t - y*u - z*s - w*s, z*v + z*r \
+ w*r + w*s, y*t - z*t - z*u - w*s, y*r + y*s - z*r, z*w + w^2 + t^2 + t*u, y*\
z + y*w - z^2 - z*w + t^2 - t*u, y^2 + z*w - w^2 - t*u - u^2 + v*s + r*s, 2*x^\
2 + y*u - y*r - z*s + w*t - w*u + w*v + w*s]);

C5 := Curve(Reduction(C,5)); 
S := AutomorphismGroup(C5); 

auts := [];
Stemp := Automorphisms(C5);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotients by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C5,[g]);
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
\*
[
    Curve over GF(5) defined by
    x[3]*x[5] + 4*x[1]*x[12] + 2*x[1]*x[13] + x[1]*x[14] + 3*x[1]*x[15],
    x[3]*x[6] + 4*x[1]*x[13] + 2*x[1]*x[14] + x[1]*x[15] + 3*x[1]*x[16],
    3*x[1]*x[2] + x[3]*x[7] + 4*x[1]*x[14] + 2*x[1]*x[15] + x[1]*x[16],
    x[1]*x[2] + x[3]*x[8] + 3*x[2]*x[10] + 4*x[1]*x[15] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[3]*x[9] + x[2]*x[10] + 3*x[2]*x[11] + 4*x[1]*x[16],
    4*x[1]*x[4] + x[3]*x[10],
    4*x[1]*x[5] + x[3]*x[11],
    4*x[1]*x[6] + x[3]*x[12],
    4*x[1]*x[7] + x[3]*x[13],
    4*x[1]*x[8] + x[3]*x[14],
    4*x[1]*x[9] + x[3]*x[15],
    4*x[9]*x[10] + x[3]*x[16],
    x[4]^2 + 4*x[1]*x[12] + 2*x[1]*x[13] + x[1]*x[14] + 3*x[1]*x[15],
    x[4]*x[5] + 4*x[1]*x[13] + 2*x[1]*x[14] + x[1]*x[15] + 3*x[1]*x[16],
    3*x[1]*x[2] + x[4]*x[6] + 4*x[1]*x[14] + 2*x[1]*x[15] + x[1]*x[16],
    x[1]*x[2] + x[4]*x[7] + 3*x[2]*x[10] + 4*x[1]*x[15] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[4]*x[8] + x[2]*x[10] + 3*x[2]*x[11] + 4*x[1]*x[16],
    4*x[1]*x[2] + x[4]*x[9] + 2*x[2]*x[10] + x[2]*x[11] + 3*x[2]*x[12],
    4*x[1]*x[5] + x[4]*x[10],
    4*x[1]*x[6] + x[4]*x[11],
    4*x[1]*x[7] + x[4]*x[12],
    4*x[1]*x[8] + x[4]*x[13],
    4*x[1]*x[9] + x[4]*x[14],
    4*x[9]*x[10] + x[4]*x[15],
    4*x[2]*x[3] + x[4]*x[16],
    3*x[1]*x[2] + x[5]^2 + 4*x[1]*x[14] + 2*x[1]*x[15] + x[1]*x[16],
    x[1]*x[2] + x[5]*x[6] + 3*x[2]*x[10] + 4*x[1]*x[15] + 2*x[1]*x[16],
    2*x[1]*x[2] + x[5]*x[7] + x[2]*x[10] + 3*x[2]*x[11] + 4*x[1]*x[16],
    4*x[1]*x[2] + x[5]*x[8] + 2*x[2]*x[10] + x[2]*x[11] + 3*x[2]*x[12],
    x[5]*x[9] + 4*x[2]*x[10] + 2*x[2]*x[11] + x[2]*x[12] + 3*x[2]*x[13],
    4*x[1]*x[6] + x[5]*x[10],
    4*x[1]*x[7] + x[5]*x[11],
    4*x[1]*x[8] + x[5]*x[12],
    4*x[1]*x[9] + x[5]*x[13],
    4*x[9]*x[10] + x[5]*x[14],
    4*x[2]*x[3] + x[5]*x[15],
    4*x[2]*x[4] + x[5]*x[16],
    2*x[1]*x[2] + x[6]^2 + x[2]*x[10] + 3*x[2]*x[11] + 4*x[1]*x[16],
    4*x[1]*x[2] + x[6]*x[7] + 2*x[2]*x[10] + x[2]*x[11] + 3*x[2]*x[12],
    x[6]*x[8] + 4*x[2]*x[10] + 2*x[2]*x[11] + x[2]*x[12] + 3*x[2]*x[13],
    x[6]*x[9] + 4*x[2]*x[11] + 2*x[2]*x[12] + x[2]*x[13] + 3*x[2]*x[14],
    4*x[1]*x[7] + x[6]*x[10],
    4*x[1]*x[8] + x[6]*x[11],
    4*x[1]*x[9] + x[6]*x[12],
    4*x[9]*x[10] + x[6]*x[13],
    4*x[2]*x[3] + x[6]*x[14],
    4*x[2]*x[4] + x[6]*x[15],
    4*x[2]*x[5] + x[6]*x[16],
    x[7]^2 + 4*x[2]*x[10] + 2*x[2]*x[11] + x[2]*x[12] + 3*x[2]*x[13],
    x[7]*x[8] + 4*x[2]*x[11] + 2*x[2]*x[12] + x[2]*x[13] + 3*x[2]*x[14],
    x[7]*x[9] + 4*x[2]*x[12] + 2*x[2]*x[13] + x[2]*x[14] + 3*x[2]*x[15],
    4*x[1]*x[8] + x[7]*x[10],
    4*x[1]*x[9] + x[7]*x[11],
    4*x[9]*x[10] + x[7]*x[12],
    4*x[2]*x[3] + x[7]*x[13],
    4*x[2]*x[4] + x[7]*x[14],
    4*x[2]*x[5] + x[7]*x[15],
    4*x[2]*x[6] + x[7]*x[16],
    x[8]^2 + 4*x[2]*x[12] + 2*x[2]*x[13] + x[2]*x[14] + 3*x[2]*x[15],
    x[8]*x[9] + 4*x[2]*x[13] + 2*x[2]*x[14] + x[2]*x[15] + 3*x[2]*x[16],
    4*x[1]*x[9] + x[8]*x[10],
    4*x[9]*x[10] + x[8]*x[11],
    4*x[2]*x[3] + x[8]*x[12],
    4*x[2]*x[4] + x[8]*x[13],
    4*x[2]*x[5] + x[8]*x[14],
    4*x[2]*x[6] + x[8]*x[15],
    4*x[2]*x[7] + x[8]*x[16],
    3*x[2]^2 + x[9]^2 + 4*x[2]*x[14] + 2*x[2]*x[15] + x[2]*x[16],
    x[3]*x[4] + 4*x[1]*x[11] + 2*x[1]*x[12] + x[1]*x[13] + 3*x[1]*x[14],
    4*x[2]*x[3] + x[9]*x[11],
    4*x[2]*x[4] + x[9]*x[12],
    4*x[2]*x[5] + x[9]*x[13],
    4*x[2]*x[6] + x[9]*x[14],
    4*x[2]*x[7] + x[9]*x[15],
    4*x[2]*x[8] + x[9]*x[16],
    x[10]^2 + 4*x[1]*x[11],
    x[10]*x[11] + 4*x[1]*x[12],
    x[10]*x[12] + 4*x[1]*x[13],
    x[10]*x[13] + 4*x[1]*x[14],
    x[10]*x[14] + 4*x[1]*x[15],
    x[10]*x[15] + 4*x[1]*x[16],
    4*x[1]*x[2] + x[10]*x[16],
    x[11]^2 + 4*x[1]*x[13],
    x[11]*x[12] + 4*x[1]*x[14],
    x[11]*x[13] + 4*x[1]*x[15],
    x[11]*x[14] + 4*x[1]*x[16],
    4*x[1]*x[2] + x[11]*x[15],
    4*x[2]*x[10] + x[11]*x[16],
    x[12]^2 + 4*x[1]*x[15],
    x[12]*x[13] + 4*x[1]*x[16],
    4*x[1]*x[2] + x[12]*x[14],
    4*x[2]*x[10] + x[12]*x[15],
    4*x[2]*x[11] + x[12]*x[16],
    4*x[1]*x[2] + x[13]^2,
    4*x[2]*x[10] + x[13]*x[14],
    4*x[2]*x[11] + x[13]*x[15],
    4*x[2]*x[12] + x[13]*x[16],
    4*x[2]*x[11] + x[14]^2,
    4*x[2]*x[12] + x[14]*x[15],
    4*x[2]*x[13] + x[14]*x[16],
    4*x[2]*x[13] + x[15]^2,
    4*x[2]*x[14] + x[15]*x[16],
    4*x[2]*x[15] + x[16]^2,
    x[3]^2 + 4*x[1]*x[10] + 2*x[1]*x[11] + x[1]*x[12] + 3*x[1]*x[13]
]
*/
#EllipticCurve(Curve(l[1])); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
