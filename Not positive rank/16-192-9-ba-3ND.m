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

*/

#EllipticCurve(Curve(Reduction(l[1],3))); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
