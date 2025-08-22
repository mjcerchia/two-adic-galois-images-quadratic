//16.96.5.v.2

P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[2*y^2 - w^2 + w*t, 2*z^2 - w*t - t^2, 2*x^2 + y*z]);

//Cp:=Curve(Reduction(C,5));

G:=AutomorphismGroup(C);
S:=Automorphisms(C);

for s in S do

s1:=G!s;
 
Order(s1);
if Order(s1) eq 2 then

AG := AutomorphismGroup(C,[s1]);
CG,prj := CurveQuotient(AG);

if Genus(CG) eq 1 then

CG;
Cp:=Curve(Reduction(CG,3));
E:=EllipticCurve(Cp);
#Points(E);
end if;

end if;
print ".......";
end for;


//128.2.a.a
E := EllipticCurve([0, 1, 0, -9, 7]);
print "Points of Jacobian factor ";
#EllipticCurve(Curve(Reduction(E,3)));
