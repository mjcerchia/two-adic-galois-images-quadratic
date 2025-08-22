//16.96.5.dx.1

P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[-x*w + y*z, 4*x^2 + y^2 - 4*z^2 - w^2 - w*t, 4*x^2 - y^2 + 4*z^2 + 2*w^2 + 3*w*t - t^2]);

Cp:=Curve(Reduction(C,3));

G:=AutomorphismGroup(Cp);
S:=Automorphisms(Cp);

for s in S do

s1:=G!s;
 
Order(s1);
if Order(s1) eq 2 then

AG := AutomorphismGroup(Cp,[s1]);
CG,prj := CurveQuotient(AG);

if Genus(CG) eq 1 then

CG;

end if;

end if;
print ".......";
end for;
