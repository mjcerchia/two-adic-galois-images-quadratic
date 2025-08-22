//27.108.7.g.1

P<x,y,z,w,t,u,v> := ProjectiveSpace(Rationals(), 6);

C := Curve(P,[x*y - x*z + w*t, 2*y*w - y*u - y*v + z*w - z*u, x*t + x*u + 2*x*v + y^2 - y*z - z^2, -3*x*w - x*t - x*u + x*v + y*z, y*w + y*t + 2*y*u - z*w - z*t - 3*z*u + z*v, -x*t - 4*x*u + x*v + y^2, 3*w^2 - w*t + 2*w*u - 2*w*v - t*u - t*v - 2*u^2 - u*v, 2*x*y - 2*x*z + 3*w^2 - 3*w*t - 2*w*u - w*v - t*v + u^2 - v^2, 5*x*y + 4*x*z + 3*w^2 - w*t + 3*w*u + t^2 + 3*t*u + 2*u^2 + 2*u*v - v^2, 9*x^2 - y*w + y*v + z*w + z*t + 2*z*u]);

Cp:=Curve(Reduction(C,7));

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

#Points(CG);
end if;

end if;
print ".......";
end for;

//243.2.a.a
E := EllipticCurve([0, 0, 1, 0, -1]);
print "Points of Jacobian factor ";
#EllipticCurve(Curve(Reduction(E,7)));

//243.2.a.b
E := EllipticCurve([0, 0, 1, 0, -61]);
print "Points of Jacobian factor ";
#EllipticCurve(Curve(Reduction(E,7)));


