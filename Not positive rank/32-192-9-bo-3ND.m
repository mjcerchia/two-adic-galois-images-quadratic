32.192.9.bo.3
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.192.9.bo.3.
We compute the automorphism group over F3 and find there is one genus one quotient by an involution. It has 4 points.
The single rank one elliptic curve factor of the Jacobian has 6 points over F3.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,r,s,v>:=ProjectiveSpace(Rationals(),8);

C := Curve(P,[x*v - u*s, v^2 + r^2 - r*s, x*v + y*v - u*r, y*r + u*v, w*r + \
t*s, x*v - w*r - t*r, x*r + t*v, x^2 - w*u, x*y + y^2 + u^2, x^2 + x*y + t*u, \
x^2 + w*t + t^2, x*s - w*v, x*r - x*s - y*s, y*s + z^2 - t*v, x*u - y*t, x*w +\
 x*t + y*w, x^2 - x*y + w^2 - w*t + w*u - t^2 - t*u + s^2, 2*x*y - y^2 - w^2 +\
 t^2 - t*u - u^2 + r*s - s^2, x^2 - x*y - y^2 + 3*u^2 + v^2, x*w + 2*x*u + 2*y\
*t + v*s, y*w + 4*y*u + v*r - v*s]);
 
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
