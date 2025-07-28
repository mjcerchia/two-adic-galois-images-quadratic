/****************************************************************************** 
 Here is a summary of the argument.

Let C be the modular curve with lmfdb label  16.96.5.l.2. 
We find there are three genus one quotients by an involution. 
Two of these quotients have a different number of points mod 7 than the single rank 1 elliptic curve 
factor of Jac(C). 
We find a point on the remaining curve, from which we are able to construct an elliptic curve, which
we find to have rank 1. 
******************************************************************************/
//Model obtained from the LMFDB
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
C := Curve(P,[y*w + y*t - z*w + z*t, 2*x^2 + y*t + z*t, 4*y*z + w^2 + t^2]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

/*
There are 3 genus one quotients, all in l. 
[
    Curve over Rational Field defined by
    x[1]^2 - x[4]^2 - 4*x[6]^2 - 4*x[7]^2,
    x[1]*x[2] - x[6]^2,
    x[1]*x[3] - x[7]^2,
    x[1]*x[5] - x[4]*x[6],
    x[1]*x[6] - x[4]*x[7],
    -x[4]*x[6] + x[1]*x[7] - 4*x[3]*x[7] - 4*x[6]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    4*x[2]^2 + x[5]^2 - x[6]^2 + 4*x[8]^2,
    x[2]*x[3] - x[8]^2,
    x[2]*x[4] - x[5]*x[6],
    x[2]*x[6] - x[5]*x[8],
    x[2]*x[7] - x[6]*x[8],
    x[5]*x[6] - x[6]*x[7] + 4*x[2]*x[8] + 4*x[3]*x[8],
    4*x[3]^2 + x[6]^2 - x[7]^2 + 4*x[8]^2,
    x[3]*x[4] - x[6]*x[7],
    x[3]*x[5] - x[6]*x[8],
    x[3]*x[6] - x[7]*x[8],
    x[4]*x[5] - x[4]*x[7] + 4*x[5]*x[8] + 4*x[7]*x[8],
    -x[6]^2 + x[4]*x[8],
    -x[6]^2 + x[5]*x[7],
    Curve over Rational Field defined by
    x[1]^2 - x[3]^2 - x[7]^2 + x[7]*x[8],
    4*x[1]*x[2] + x[7]*x[8],
    4*x[1]*x[4] + x[5]*x[7],
    x[1]*x[5] - x[3]*x[7],
    4*x[1]*x[6] + x[7]^2,
    -x[3]*x[5] + x[1]*x[7] + 4*x[6]*x[7] - 4*x[6]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    16*x[2]^2 + 4*x[4]^2 + x[7]*x[8] + 4*x[8]^2,
    4*x[2]*x[3] + x[5]*x[8],
    x[2]*x[5] - x[4]*x[8],
    4*x[2]*x[6] + x[8]^2,
    x[2]*x[7] - x[6]*x[8],
    x[4]*x[5] - x[6]*x[7] + 4*x[2]*x[8] - 4*x[6]*x[8],
    4*x[3]*x[4] + x[7]^2 + 4*x[7]*x[8] - 4*x[8]^2,
    4*x[3]*x[6] + x[5]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    4*x[4]*x[6] + x[5]*x[8],
    -x[5]*x[6] + x[4]*x[7],
    x[5]^2 - x[7]^2 - 4*x[7]*x[8] + 4*x[8]^2,
    4*x[6]^2 + x[7]*x[8],
    Curve over Rational Field defined by
    x[1]^2 - x[3]^2 + x[7]^2 + x[7]*x[8],
    4*x[1]*x[2] + x[7]*x[8],
    4*x[1]*x[4] + x[5]*x[7],
    x[1]*x[5] - x[3]*x[7],
    4*x[1]*x[6] + x[7]^2,
    -x[3]*x[5] + x[1]*x[7] - 4*x[6]*x[7] - 4*x[6]*x[8],
    -x[6]*x[7] + x[1]*x[8],
    16*x[2]^2 + 4*x[4]^2 + x[7]*x[8] - 4*x[8]^2,
    4*x[2]*x[3] + x[5]*x[8],
    x[2]*x[5] - x[4]*x[8],
    4*x[2]*x[6] + x[8]^2,
    x[2]*x[7] - x[6]*x[8],
    x[4]*x[5] - x[6]*x[7] + 4*x[2]*x[8] + 4*x[6]*x[8],
    4*x[3]*x[4] + x[7]^2 - 4*x[7]*x[8] - 4*x[8]^2,
    4*x[3]*x[6] + x[5]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    4*x[4]*x[6] + x[5]*x[8],
    -x[5]*x[6] + x[4]*x[7],
    x[5]^2 - x[7]^2 + 4*x[7]*x[8] + 4*x[8]^2,
    4*x[6]^2 + x[7]*x[8]
]
*/

//The rank one factor of the Jacobian corresponding to 128.2.a.a is 
//"Elliptic Curve defined by y^2 = x^3 + x^2 + x + 1 over Rational Field"

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3+x^2+x+1);

//The first and second genus one quotient in l have a different number of pts than E mod 7
#EllipticCurve(Curve(Reduction(l[1],7))); // 8
#EllipticCurve(Curve(Reduction(l[2],7))); // 4
#EllipticCurve(Curve(Reduction(E,7))); // 12

C1 := l[3];

//Search for a rational point on C1
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
Rank(E); // 1



