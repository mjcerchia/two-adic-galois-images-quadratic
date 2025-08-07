/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.96.4.f.1.
We find there are two genus one quotients by an involution, but both are rank 0.
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w> := ProjectiveSpace(Rationals(),3);
C := Curve(P,[2*x^2 + 8*y^2 - z^2 - w^2, 2*x^2*z + 2*x^2*w - z^3 - 2*z^2*w + z*w^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are two genus one quotients by an involution
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
/*
  [
    Curve over Rational Field defined by
    2*x[1]^2 + 4*x[2]*x[4] - 10*x[4]^2 - x[5]^2 + x[5]*x[6] + 2*x[6]^2,
    x[1]*x[2] - x[4]^2,
    x[1]*x[3] - x[5]^2,
    2*x[1]*x[4] - 2*x[2]*x[4] + 4*x[4]^2 - x[5]*x[6] - x[6]^2,
    2*x[1]*x[5] - x[3]*x[5] + 4*x[4]*x[5] - x[3]*x[6] - 2*x[4]*x[6],
    -x[4]*x[5] + x[1]*x[6],
    x[2]*x[3] - x[6]^2,
    x[2]*x[5] - x[4]*x[6],
    x[3]*x[4] - x[5]*x[6],
    Curve over Rational Field defined by
    x[1]*x[2] - x[4]^2,
    x[1]*x[3] - x[5]^2,
    x[1]*x[4] - x[2]*x[4] - 2*x[4]^2 + 8*x[5]^2 + 8*x[5]*x[6],
    -x[4]*x[5] + x[1]*x[6],
    x[2]^2 + 2*x[2]*x[4] - x[4]^2 - 8*x[5]*x[6] - 8*x[6]^2,
    x[2]*x[3] - x[6]^2,
    x[2]*x[5] - x[4]*x[6],
    -8*x[3]*x[5] - x[4]*x[5] + x[2]*x[6] - 8*x[3]*x[6] + 2*x[4]*x[6],
    x[3]*x[4] - x[5]*x[6]
]
*/

C1 := l[1];
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
Rank(E); // 0

C1 := l[2];
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
Rank(E); //0
