32.192.9.bs.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.bs.2.
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,r,s,v>:=ProjectiveSpace(Rationals(),8);

C := Curve(P,[w*r - u*v, t*r + u*s, x^2 - y*v - u*s, z*r + t*v, z*r - w*s, z\
*u + w*t, x^2 + y*v + t*s - u*s, y*r + y*s + z*s, y*t - y*u + z*t, y*z + y*w +\
 z^2, y*z - y*w + z^2 + t^2, 2*v^2 + r^2 + r*s, 2*w*v + u*r + u*s, 2*w^2 - t*u\
 + u^2, y*w + 2*w^2 + t*u - u^2 - v^2, 2*z*w + t^2 - t*u, 2*z*v - t*s + u*s, y\
*t - 2*z*u + 2*w*t + v*s, y*u + 4*w*u - v*r, 2*y^2 + y*z - z^2 + 2*z*w + t*u +\
 2*u^2 - v^2 - r^2, 2*y^2 - y*z + y*w - z^2 - 2*z*w + 2*t^2 + t*u - s^2]);
 
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
\*

*/
#EllipticCurve(Curve(l[1])); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6


