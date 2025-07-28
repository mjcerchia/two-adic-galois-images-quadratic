/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 332.96.5.a.2.
We find there are three genus one quotients by an involution. 
We are able to find a point on the second curve, from which we are able to construct an elliptic curve, which
we find to have rank 1. 
******************************************************************************/

//Model obtained from the LMFDB
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[x^2 + z*w, x*z - 2*x*t - y^2, z*t + 4*w^2 - t^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//We find that there are 3 genus one quotients by involutions
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

//One of the quotients is given by
P<[x]> := ProjectiveSpace(Rationals(),7);
C1 := Curve(P,[11520*x[1]*x[2] + 180*x[4]^2 + 138*x[4]*x[5] + 35*x[5]^2 - 288*x[8]^2,
45*x[1]*x[3] - 45*x[4]*x[6] - 12*x[5]*x[6] + x[5]*x[7],
15*x[1]*x[4] - 15*x[6]^2 + x[6]*x[7],
15*x[1]*x[5] - 15*x[6]*x[7] + x[7]^2,
45*x[3]*x[4] + 2880*x[1]*x[7] - 1248*x[6]*x[8] + 64*x[7]*x[8],
30*x[4]*x[5] + 13*x[5]^2 + 1920*x[1]*x[8] - 96*x[8]^2,
384*x[2]^2 + 6*x[4]*x[5] + 5*x[5]^2 - 96*x[8]^2,
3*x[2]*x[3] - 3*x[5]*x[6] - x[5]*x[7],
x[2]*x[4] - x[6]*x[7],
x[2]*x[5] - x[7]^2,
x[3]*x[4] + 64*x[2]*x[6] - 32*x[6]*x[8],
3*x[2]*x[7] - 3*x[6]*x[8] - x[7]*x[8],
x[5]^2 + 64*x[2]*x[8] - 32*x[8]^2,
3*x[3]^2 + 192*x[6]*x[7] - 32*x[7]^2,
3*x[3]*x[5] + 192*x[6]*x[8] - 32*x[7]*x[8],
-3*x[4]^2 - x[4]*x[5] + 3*x[3]*x[6],
-3*x[4]*x[5] - x[5]^2 + 3*x[3]*x[7],
-x[5]*x[7] + x[3]*x[8],
-x[5]*x[6] + x[4]*x[7],
-3*x[7]^2 + 3*x[4]*x[8] + x[5]*x[8]]);


//We can find a rational point on C1
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

E := EllipticCurve(C1,pt);
Rank(E); //1
