32.192.9.l.3
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.l.3.
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,r,s,v>:=ProjectiveSpace(Rationals(),8);

C := Curve(P,[y^2 - t*r, y^2 - y*w + t*s, y*z - t*v, v^2 + v*r + r^2 + r*s +\
 s^2, y^2 + u*r + v^2 - v*s - r*s, y^2 - u*r + v*r + v*s - s^2, y^2 + y*z + y*\
w + u*v - u*s, y*w + t^2 + u^2 - u*s, y*z - t*s - u*v - u*s + v^2 + s^2, y*z +\
 z^2 + w^2 - t*s, y*v - z*r, y*z + y*w + z^2 + z*w + t^2 - t*u, z*w + w^2 - t^\
2 - t*u, x^2 + y*u + z*s, y*v + z*s - w*v, y*v + y*r + z*v + w*s, y*r + y*s - \
w*r, y*t + y*u - w*v - w*r - w*s, y*t - y*u + y*s + w*v - w*s, y*s + z*t + z*u\
 - z*v + z*s + w*t - w*u, z*t - z*u + z*s - w*t - w*u + w*s]);
 
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
