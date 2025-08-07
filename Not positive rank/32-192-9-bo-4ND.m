32.192.9.bo.4
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.bo.4.
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,r,s,v>:=ProjectiveSpace(Rationals(),8);

C := Curve(P,[u*r + v*s, y^2 + v^2 - v*r, y^2 - y*w + u*v + v*s, t*s + u*v +\
 v*s, y^2 + t*r, y*z + t*v, y*z - z*w + t*u, y*z - z^2 - t^2, y*u + z*s, y*u -\
 y*v + w*v, y*v - z*r, y*t + y*v - z*v, y*t + y*u - z*u - w*t, y*r + y*s - w*r\
, y*w - w^2 - t*s - u^2 - u*s, x^2 + y*t - y*u + w*v - w*r - w*s, x^2 - y*t - \
y*u - y*s + w*v + w*s, y*z + w^2 - t*v + t*s + u^2 - u*v - u*r - u*s + v^2, z^\
2 - w^2 - t^2 + t*v + u^2 + u*s, y^2 - t*r + r^2 + 2*r*s + 2*s^2, y*v - 2*z*t \
- 2*w*u]);
 
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

