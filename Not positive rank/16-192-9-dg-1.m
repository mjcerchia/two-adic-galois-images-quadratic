/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.192.9.dg.1
We compute the automorphism group over Q and find there are no genus one quotient by an involution. 
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t,u,v,r,s> := ProjectiveSpace(Rationals(),8);
C := Curve(P,[x*y - z^2, x^2 + x*y + u*r, x*y + x*v - t*u, x*y - x*v + w*r, x^2 - x*y - w*t, x*t + y*t + y*r + v*r, x*w + y*w - y*u + u*v, x*r + y*t - y*r - t*v, x*w + x*u + 2*y*r - r*s, x*y + z^2 - w*u - t*r + v*s, x^2 - 2*y^2 + v^2, x*t - x*r + 2*y*u + u*s, x*y - x*s + y*s + w^2 - w*t - w*u - v^2, x*y + x*s + y*s - t*r - u*r + v^2 + r^2, x*y - x*s - y*s - w*u - u^2 - u*r + v^2, x^2 - x*y + w*t + u^2 + v*s - r^2 - s^2, x^2 + x*y + w^2 - t^2 - u*r - v*s - s^2, x*u - y*w - y*u - w*v, x*s + 2*y*v + w*u - t*r, x*w - x*u + 2*y*t - t*s, x*t + x*r - 2*y*w - w*s]);

S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;


//There is no genus one quotient by an involution
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
