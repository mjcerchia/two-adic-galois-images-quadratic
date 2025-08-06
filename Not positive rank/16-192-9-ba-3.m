16.192.9.ba.3
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.9.ba.3.
We compute the automorphism group over Q and find there is one genus one quotient by an involution. Over F3, it has 4 points.
But the single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,r,u,s,v> := ProjectiveSpace(Rationals(),8);
C := Curve(P,[w*u - w*s - u^2 + r*s, y^2 + y*z + y*u + y*s + z^2 + u^2 + s^2\
, y*r - y*s - w*u + w*s - u^2 + r*s, y^2 - y*w - y*u - y*r + w^2 + u^2 + r^2, \
y*z + y*u + y*r + z^2 + t^2 + v^2 - r^2, y*u + y*s - w*s + t*r - v*s + s^2, y*\
t + y*r + w*s - t*r - v*r - r^2, y*w + w*s + t*r - v*s + s^2, 2*x^2 + t^2 - u^\
2 + v^2 + v*r - v*s, y^2 + y*z + w*s - t*s - u*r - u*s - v*r - r^2, y*w - w^2 \
+ w*s + t*s - u^2 - u*r - u*s + v*r, y*u + y*r - w*u + t*u - u*v + u*s, y^2 - \
y*u - y*r + w*u + t*u + u*v - u*s, y*w - w*t + t*u + t*r, y*u + y*r + z*r - t*\
s - v*r - r^2, y*u + y*s + z*u - t*s + u^2 - v*r - r^2 - r*s, y^2 + y*z + z*t \
+ t*u - t*s, y*u + y*r + z*w - w*u - t*s - v*r - r^2, y*v - t^2 + u^2 - v^2 - \
v*r + v*s, w*u + w*v + w*r + t*u - v*r - r^2, y*t + y*s - z*v + z*s + w*u - w*\
s + t*u + t*s]);


S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There are four genus one quotients by an involution
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
\*
[
    Curve over Rational Field defined by
    x[3]*x[5] + x[1]*x[11] + 6*x[1]*x[13] + x[1]*x[15],
    x[3]*x[6] + x[1]*x[12] + 6*x[1]*x[14] + x[1]*x[16],
    x[1]*x[2] + x[3]*x[7] + x[1]*x[13] + 6*x[1]*x[15],
    x[3]*x[8] + x[2]*x[10] + x[1]*x[14] + 6*x[1]*x[16],
    6*x[1]*x[2] + x[3]*x[9] + x[2]*x[11] + x[1]*x[15],
    -x[1]*x[4] + x[3]*x[10],
    -x[1]*x[5] + x[3]*x[11],
    -x[1]*x[6] + x[3]*x[12],
    -x[1]*x[7] + x[3]*x[13],
    -x[1]*x[8] + x[3]*x[14],
    -x[1]*x[9] + x[3]*x[15],
    -x[9]*x[10] + x[3]*x[16],
    x[4]^2 + x[1]*x[11] + 6*x[1]*x[13] + x[1]*x[15],
    x[4]*x[5] + x[1]*x[12] + 6*x[1]*x[14] + x[1]*x[16],
    x[1]*x[2] + x[4]*x[6] + x[1]*x[13] + 6*x[1]*x[15],
    x[4]*x[7] + x[2]*x[10] + x[1]*x[14] + 6*x[1]*x[16],
    6*x[1]*x[2] + x[4]*x[8] + x[2]*x[11] + x[1]*x[15],
    x[4]*x[9] + 6*x[2]*x[10] + x[2]*x[12] + x[1]*x[16],
    -x[1]*x[5] + x[4]*x[10],
    -x[1]*x[6] + x[4]*x[11],
    -x[1]*x[7] + x[4]*x[12],
    -x[1]*x[8] + x[4]*x[13],
    -x[1]*x[9] + x[4]*x[14],
    -x[9]*x[10] + x[4]*x[15],
    -x[2]*x[3] + x[4]*x[16],
    x[1]*x[2] + x[5]^2 + x[1]*x[13] + 6*x[1]*x[15],
    x[5]*x[6] + x[2]*x[10] + x[1]*x[14] + 6*x[1]*x[16],
    6*x[1]*x[2] + x[5]*x[7] + x[2]*x[11] + x[1]*x[15],
    x[5]*x[8] + 6*x[2]*x[10] + x[2]*x[12] + x[1]*x[16],
    x[1]*x[2] + x[5]*x[9] + 6*x[2]*x[11] + x[2]*x[13],
    -x[1]*x[6] + x[5]*x[10],
    -x[1]*x[7] + x[5]*x[11],
    -x[1]*x[8] + x[5]*x[12],
    -x[1]*x[9] + x[5]*x[13],
    -x[9]*x[10] + x[5]*x[14],
    -x[2]*x[3] + x[5]*x[15],
    -x[2]*x[4] + x[5]*x[16],
    6*x[1]*x[2] + x[6]^2 + x[2]*x[11] + x[1]*x[15],
    x[6]*x[7] + 6*x[2]*x[10] + x[2]*x[12] + x[1]*x[16],
    x[1]*x[2] + x[6]*x[8] + 6*x[2]*x[11] + x[2]*x[13],
    x[6]*x[9] + x[2]*x[10] + 6*x[2]*x[12] + x[2]*x[14],
    -x[1]*x[7] + x[6]*x[10],
    -x[1]*x[8] + x[6]*x[11],
    -x[1]*x[9] + x[6]*x[12],
    -x[9]*x[10] + x[6]*x[13],
    -x[2]*x[3] + x[6]*x[14],
    -x[2]*x[4] + x[6]*x[15],
    -x[2]*x[5] + x[6]*x[16],
    x[1]*x[2] + x[7]^2 + 6*x[2]*x[11] + x[2]*x[13],
    x[7]*x[8] + x[2]*x[10] + 6*x[2]*x[12] + x[2]*x[14],
    x[7]*x[9] + x[2]*x[11] + 6*x[2]*x[13] + x[2]*x[15],
    -x[1]*x[8] + x[7]*x[10],
    -x[1]*x[9] + x[7]*x[11],
    -x[9]*x[10] + x[7]*x[12],
    -x[2]*x[3] + x[7]*x[13],
    -x[2]*x[4] + x[7]*x[14],
    -x[2]*x[5] + x[7]*x[15],
    -x[2]*x[6] + x[7]*x[16],
    x[8]^2 + x[2]*x[11] + 6*x[2]*x[13] + x[2]*x[15],
    x[8]*x[9] + x[2]*x[12] + 6*x[2]*x[14] + x[2]*x[16],
    -x[1]*x[9] + x[8]*x[10],
    -x[9]*x[10] + x[8]*x[11],
    -x[2]*x[3] + x[8]*x[12],
    -x[2]*x[4] + x[8]*x[13],
    -x[2]*x[5] + x[8]*x[14],
    -x[2]*x[6] + x[8]*x[15],
    -x[2]*x[7] + x[8]*x[16],
    x[2]^2 + x[9]^2 + x[2]*x[13] + 6*x[2]*x[15],
    x[3]*x[4] + x[1]*x[10] + 6*x[1]*x[12] + x[1]*x[14],
    -x[2]*x[3] + x[9]*x[11],
    -x[2]*x[4] + x[9]*x[12],
    -x[2]*x[5] + x[9]*x[13],
    -x[2]*x[6] + x[9]*x[14],
    -x[2]*x[7] + x[9]*x[15],
    -x[2]*x[8] + x[9]*x[16],
    x[10]^2 - x[1]*x[11],
    x[10]*x[11] - x[1]*x[12],
    x[10]*x[12] - x[1]*x[13],
    x[10]*x[13] - x[1]*x[14],
    x[10]*x[14] - x[1]*x[15],
    x[10]*x[15] - x[1]*x[16],
    -x[1]*x[2] + x[10]*x[16],
    x[11]^2 - x[1]*x[13],
    x[11]*x[12] - x[1]*x[14],
    x[11]*x[13] - x[1]*x[15],
    x[11]*x[14] - x[1]*x[16],
    -x[1]*x[2] + x[11]*x[15],
    -x[2]*x[10] + x[11]*x[16],
    x[12]^2 - x[1]*x[15],
    x[12]*x[13] - x[1]*x[16],
    -x[1]*x[2] + x[12]*x[14],
    -x[2]*x[10] + x[12]*x[15],
    -x[2]*x[11] + x[12]*x[16],
    -x[1]*x[2] + x[13]^2,
    -x[2]*x[10] + x[13]*x[14],
    -x[2]*x[11] + x[13]*x[15],
    -x[2]*x[12] + x[13]*x[16],
    -x[2]*x[11] + x[14]^2,
    -x[2]*x[12] + x[14]*x[15],
    -x[2]*x[13] + x[14]*x[16],
    -x[2]*x[13] + x[15]^2,
    -x[2]*x[14] + x[15]*x[16],
    -x[2]*x[15] + x[16]^2,
    x[1]^2 + x[3]^2 + 6*x[1]*x[11] + x[1]*x[13]
]
*/

#EllipticCurve(Curve(Reduction(l[1],3))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
