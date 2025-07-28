/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.a.1. 
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve, from which we are able to construct an elliptic curve, which
we find to have rank 1. 
******************************************************************************/


P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y^2 - z*w, 2*x^2 - y*z + 2*y*t, z*t + 4*w^2 - t^2]);

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

//The second of these quotients has the following model 
P<[x]> := ProjectiveSpace(Rationals(),7);
C1 := Curve(P,[240*x[1]*x[2] - 3*x[4]^2 - 6*x[6]^2 - 120*x[5]*x[7] - 52*x[5]*x[8],
15*x[1]*x[3] - 30*x[7]^2 + 2*x[7]*x[8],
15*x[1]*x[4] - 30*x[7]*x[8] + 2*x[8]^2,
45*x[1]*x[5] + x[5]*x[6] - 45*x[3]*x[7] - 12*x[4]*x[7],
-x[4]^2 + 40*x[1]*x[6] - 2*x[6]^2 - 20*x[5]*x[8],
-45*x[3]*x[5] - 78*x[6]*x[7] + 180*x[1]*x[8] + 4*x[6]*x[8],
8*x[2]^2 - x[4]^2 - 2*x[6]^2 - 4*x[5]*x[8],
x[2]*x[3] - 2*x[7]*x[8],
x[2]*x[4] - 2*x[8]^2,
3*x[2]*x[5] - x[5]*x[6] - 3*x[4]*x[7],
-x[4]^2 + 4*x[2]*x[6] - 2*x[6]^2,
-x[3]*x[5] + 4*x[2]*x[7] - 2*x[6]*x[7],
-3*x[6]*x[7] + 3*x[2]*x[8] - x[6]*x[8],
9*x[3]^2 - x[4]^2 - 18*x[5]*x[7] + 6*x[5]*x[8],
3*x[3]*x[4] + x[4]^2 - 6*x[5]*x[8],
3*x[3]*x[6] + x[4]*x[6] - 6*x[8]^2,
-x[4]*x[7] + x[3]*x[8],
3*x[4]*x[5] - 12*x[6]*x[7] + 2*x[6]*x[8],
-x[5]*x[6] + x[4]*x[8],
3*x[5]^2 - 12*x[7]*x[8] + 2*x[8]^2]);

//The following produces a rational point
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

pt := C1!pt;
E := EllipticCurve(C1,pt);
Rank(E); //1
