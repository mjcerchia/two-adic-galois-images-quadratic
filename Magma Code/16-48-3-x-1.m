16.48.3.x.1

/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.48.3.x.1. 
We find there are five genus one quotients by an involution. 
We are able to find a point on the second curve. From this, we are able to construct an elliptic curve, 
which we find to be rank 1. 
******************************************************************************/
P<x,y,z> := ProjectiveSpace(Rationals(),2);
C := Curve(P,[2*x^3*y + x^3*z + 3*x^2*y*z + 2*x*y^3 + 3*x*y^2*z - 2*x*z^3 + y^3*z - 2*y*z^3 - z^4]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are three genus one quotients by an involution
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

//One of these quotients has the following model:
P<[x]> := ProjectiveSpace(Rationals(),3);
C1 := Curve(P,[872*x[1]^2 - 1744*x[1]*x[2] + 872*x[2]^2 - 126420*x[2]*x[3] - 457905*x[3]^2 - 
    10816*x[1]*x[4] + 32768*x[2]*x[4] + 246960*x[3]*x[4] - 13184*x[4]^2,
872*x[1]*x[3] + 2544*x[2]*x[3] + 24031*x[3]^2 - 128*x[1]*x[4] - 768*x[2]*x[4] - 
    17056*x[3]*x[4] + 2816*x[4]^2]);

//We search for a point
  rationalPoints := function(D : Bound := 1)
    return {@D![t :t in tup]
              : tup in CartesianPower([-Bound..Bound],Dimension(AmbientSpace(D)\
)+1)
              | not {i : i in tup} eq {0}
                and {Evaluate(eqns,[t : t in tup]) : eqns in
DefiningEquations(D)} eq {0}
            @};
  end function;
  rationalPoints(C1:Bound := 1);
pt :=   rationalPoints(C1:Bound := 1)[1];

E:=EllipticCurve(C1,pt);
Rank(E);
