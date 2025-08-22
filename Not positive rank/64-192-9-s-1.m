//64.192.9.s.1

P<x,y,z,w,r,s,t,u,v>:=ProjectiveSpace(Rationals(),8);
C:=Curve(P,[
r*t - s*u,
x*u + z*t,
-x*u + y*r,
x*r + z*s,
-x*t + y*s,
-x*t + z*u + w*u,
-x*s + z*r + w*r,
x^2 + z^2 + z*w,
2*x*y + s*u,
-2*y*z + r*u,
2*y*w + r*u + s*t,
y*z - y*w + r*u - s*t + v^2,
x^2 - z^2 - z*w + r*s,
-2*x*z + r^2,
-x*z + x*w + y^2 - r^2,
x*z - x*w + y^2 + r^2 + t*u,
2*x*w + r^2 + s^2,
x*r + y*t - 2*z*s + w*s,
-x^2 + z^2 - z*w - 2*w^2 + 2*r*s + t^2,
x*s - y*u + 2*z*r - 2*w*r,
-x^2 - 3*z^2 + 3*z*w - 2*w^2 + r*s + t^2 + u^2]);
Cp:=Curve(Reduction(C,3));

G:=AutomorphismGroup(Cp);
S:=Automorphisms(Cp);

for s in S do

s1:=G!s;
 
Order(s1);
if Order(s1) eq 2 then

AG := AutomorphismGroup(Cp,[s1]);
CG,prj := CurveQuotient(AG);
Genus(CG);
if Genus(CG) eq 1 then


EG:=EllipticCurve(CG);
#Points(EG);
end if;

end if;
print ".......";
end for;

// There is exactly one genus 1 factor, it has 4 points mod 3

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);
print "Points of Jacobian factor ";
#EllipticCurve(Curve(Reduction(E,3))); // this has 6 points
