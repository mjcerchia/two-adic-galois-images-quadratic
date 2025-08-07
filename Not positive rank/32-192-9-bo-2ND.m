32.192.9.bo.2
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.bo.2.
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
\*

*/
#EllipticCurve(Curve(l[1])); //4

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
#EllipticCurve(Curve(Reduction(E,3))); //6
