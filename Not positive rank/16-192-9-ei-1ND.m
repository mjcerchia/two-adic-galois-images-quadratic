16.192.9.ei.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.9.ei.1
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F5, has 4 points.
The single rank one elliptic curve factor of the Jacobian has 10 points over F5.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,v,r,s> := ProjectiveSpace(Rationals(),8);
C:= Curve(P,[x*v - u*r, x*v - t*s - u*s, z*v - t*s + v*s, w*v - t*r + v*r, x*t - x*v - w*u, z*r - w*s, z*t + z*u - w*u, x*r - x*s - w*s - r*s, z*s + w*s + u*v - r*s - s^2, z*s - w*r - t*v + r^2 - s^2, x*z - x*w + z*w + w*s, z^2 + z*w + t*u - r*s - s^2, z^2 - w^2 - t^2 + r^2 - s^2, 2*x*u - z*t + z*v - v*s, 2*x^2 - z*w - r*s, x*z + x*w - x*r - x*s + t*u + u^2, x*z + x*w + x*r + x*s + z*s - w*r + t*u + u^2 - v^2, x*z - x*w - x*r + x*s + z^2 - z*w - t*u - t*v - u^2 + v^2 + r^2 - r*s, x*z - x*w - x*r + x*s - z*w + z*s + w^2 - w*r + t*u + u^2 - v^2 + r^2 - r*s, x*t + 2*y^2, x*t + x*v + z*t - z*v - w*t + w*u + w*v - v*r + v*s]);
C13 := Curve(Reduction(C,13)); 
S := AutomorphismGroup(C13); 

S := AutomorphismGroup(C13); 
auts := [];
Stemp := Automorphisms(C13);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is one genus one quotients by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C13,[g]);
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
    Curve over Rational Field defined by
    2*x[1]^2 - x[3]^2 + x[11]^2 + 2*x[15]^2,
    x[1]*x[2] - x[15]^2,
    x[1]*x[4] - x[11]^2 - 2*x[15]^2,
    x[1]*x[5] - x[3]*x[14],
    x[1]*x[6] - x[11]*x[12] - 2*x[15]*x[16],
    x[1]*x[7] - x[12]^2 - 2*x[16]^2,
    x[1]*x[8] - x[5]*x[14],
    x[1]*x[9] - x[10]*x[11],
    x[1]*x[10] - x[8]*x[14],
    x[1]*x[11] - x[10]*x[14],
    x[1]*x[12] - x[11]*x[14],
    x[1]*x[13] - x[12]*x[14],
    -x[3]*x[5] + x[11]*x[12] + 2*x[1]*x[14] + 2*x[15]*x[16],
    -x[14]^2 + x[1]*x[15],
    -x[14]*x[15] + x[1]*x[16],
    2*x[2]^2 - x[11]^2 + x[13]^2,
    x[2]*x[3] - x[10]*x[14],
    x[2]*x[4] - x[11]^2,
    x[2]*x[5] - x[11]*x[14],
    x[2]*x[6] - x[11]*x[12],
    x[2]*x[7] - x[12]^2,
    x[2]*x[8] - x[12]*x[14],
    x[2]*x[9] - x[12]*x[13],
    x[2]*x[10] - x[13]*x[14],
    x[2]*x[11] - x[13]*x[15],
    x[2]*x[12] - x[13]*x[16],
    x[2]*x[14] - x[15]*x[16],
    x[2]*x[15] - x[16]^2,
    -x[10]*x[11] + x[12]*x[13] + 2*x[2]*x[16],
    x[3]*x[4] - x[9]*x[12] - 2*x[10]*x[14] - 2*x[13]*x[15],
    x[3]*x[6] - x[9]*x[13] - 2*x[11]*x[14] - 2*x[13]*x[16],
    x[3]*x[7] - x[9]*x[10] - 2*x[12]*x[14],
    x[3]*x[8] - x[12]^2 - 2*x[14]^2 - 2*x[16]^2,
    x[3]*x[9] - x[9]*x[11] - 2*x[13]*x[14],
    x[3]*x[10] - x[10]*x[11] - 2*x[14]*x[15],
    x[3]*x[11] - x[11]^2 - 2*x[15]^2,
    x[3]*x[12] - x[11]*x[12] - 2*x[15]*x[16],
    -x[12]^2 + x[3]*x[13] - 2*x[16]^2,
    -x[5]*x[14] + x[3]*x[15],
    -x[8]*x[14] + x[3]*x[16],
    x[4]^2 - x[7]^2 - 2*x[11]^2,
    x[4]*x[5] - x[9]*x[13] - 2*x[11]*x[14] - 2*x[13]*x[16],
    x[4]*x[6] - x[7]*x[9] - 2*x[11]*x[12],
    x[4]*x[7] - x[9]^2 - 2*x[12]^2,
    x[4]*x[8] - x[9]*x[10] - 2*x[12]*x[14],
    -x[6]*x[7] + x[4]*x[9],
    x[4]*x[10] - x[9]*x[11] - 2*x[13]*x[14],
    x[4]*x[11] - x[9]*x[12] - 2*x[13]*x[15],
    x[4]*x[12] - x[9]*x[13] - 2*x[13]*x[16],
    -x[9]*x[10] + x[4]*x[13],
    -x[11]*x[12] + x[4]*x[14] - 2*x[15]*x[16],
    -x[12]^2 + x[4]*x[15] - 2*x[16]^2,
    -x[10]*x[11] + x[4]*x[16],
    x[5]^2 - x[12]^2 - 2*x[14]^2 - 2*x[16]^2,
    x[5]*x[6] - x[9]*x[10] - 2*x[12]*x[14],
    x[5]*x[7] - x[9]*x[11] - 2*x[13]*x[14],
    x[5]*x[8] - x[10]*x[11] - 2*x[14]*x[15],
    x[5]*x[9] - x[9]*x[12] - 2*x[13]*x[15],
    x[5]*x[10] - x[11]^2 - 2*x[15]^2,
    x[5]*x[11] - x[11]*x[12] - 2*x[15]*x[16],
    x[5]*x[12] - x[12]^2 - 2*x[16]^2,
    -x[10]*x[11] + x[5]*x[13],
    -x[8]*x[14] + x[5]*x[15],
    -x[10]*x[14] + x[5]*x[16],
    x[6]^2 - x[9]^2 - 2*x[12]^2,
    x[6]*x[8] - x[9]*x[11] - 2*x[13]*x[14],
    -x[7]^2 + x[6]*x[9],
    x[6]*x[10] - x[9]*x[12] - 2*x[13]*x[15],
    x[6]*x[11] - x[9]*x[13] - 2*x[13]*x[16],
    -x[9]*x[10] + x[6]*x[12],
    -x[9]*x[11] + x[6]*x[13],
    -x[12]^2 + x[6]*x[14] - 2*x[16]^2,
    -x[10]*x[11] + x[6]*x[15],
    -x[11]^2 + x[6]*x[16],
    x[7]*x[8] - x[9]*x[12] - 2*x[13]*x[15],
    x[7]*x[10] - x[9]*x[13] - 2*x[13]*x[16],
    -x[9]*x[10] + x[7]*x[11],
    -x[9]*x[11] + x[7]*x[12],
    -x[9]*x[12] + x[7]*x[13],
    -x[10]*x[11] + x[7]*x[14],
    -x[11]^2 + x[7]*x[15],
    -x[11]*x[12] + x[7]*x[16],
    x[8]^2 - x[11]^2 - 2*x[15]^2,
    x[8]*x[9] - x[9]*x[13] - 2*x[13]*x[16],
    x[8]*x[10] - x[11]*x[12] - 2*x[15]*x[16],
    x[8]*x[11] - x[12]^2 - 2*x[16]^2,
    -x[10]*x[11] + x[8]*x[12],
    -x[11]^2 + x[8]*x[13],
    -x[10]*x[14] + x[8]*x[15],
    -x[11]*x[14] + x[8]*x[16],
    -x[11]^2 + x[9]*x[14],
    -x[11]*x[12] + x[9]*x[15],
    -x[12]^2 + x[9]*x[16],
    x[10]^2 - x[12]^2 - 2*x[16]^2,
    -x[11]^2 + x[10]*x[12],
    -x[11]*x[12] + x[10]*x[13],
    -x[11]*x[14] + x[10]*x[15],
    -x[12]*x[14] + x[10]*x[16],
    -x[12]^2 + x[11]*x[13],
    -x[12]*x[14] + x[11]*x[15],
    -x[13]*x[14] + x[11]*x[16],
    -x[13]*x[14] + x[12]*x[15],
    -x[13]*x[15] + x[12]*x[16],
    -x[15]^2 + x[14]*x[16]
]
*/

#EllipticCurve(Curve(Reduction(l[1],5))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10
