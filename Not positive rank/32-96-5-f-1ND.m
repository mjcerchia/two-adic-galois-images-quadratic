32.96.5.f.1
/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 32.96.5.f.1.
We compute the automorphism group over Q and find there are three genus one quotient by an involution. Over F5, the first and second
have 4 and 2 points, respectively, while the single rank one elliptic curve factor of the Jacobian has 10 points over F5. This leaves us 
with the second curve. We construct a new model of this quotient curve using Riemann-Roch and find that it fails local solubility. 
NOT bielliptic. 
******************************************************************************/
P<x,y,z,w,t>:=ProjectiveSpace(Rationals(),4);
C := Curve(P,[2*x^2 + w*t, 2*y*z - w^2, 2*y^2 + 32*z^2 + t^2]);
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

l;
\*
[
    Curve over Rational Field defined by
    128*x[1]^2 + x[3]^2 + 2*x[7]^2,
    x[1]*x[2] - x[7]^2,
    x[1]*x[4] - x[3]*x[6],
    x[1]*x[5] - x[4]*x[6],
    x[3]*x[4] + 128*x[1]*x[6] + 2*x[7]*x[8],
    -x[6]^2 + x[1]*x[7],
    -x[6]*x[7] + x[1]*x[8],
    2*x[2]^2 + x[5]^2 + 128*x[7]^2,
    x[2]*x[3] - x[5]*x[7],
    x[2]*x[4] - x[5]*x[8],
    x[2]*x[6] - x[7]*x[8],
    x[2]*x[7] - x[8]^2,
    x[4]*x[5] + 128*x[6]*x[7] + 2*x[2]*x[8],
    x[3]*x[5] + 128*x[6]^2 + 2*x[8]^2,
    -x[4]*x[6] + x[3]*x[7],
    -x[5]*x[6] + x[3]*x[8],
    x[4]^2 + 128*x[6]^2 + 2*x[8]^2,
    -x[5]*x[6] + x[4]*x[7],
    -x[5]*x[7] + x[4]*x[8],
    -x[7]^2 + x[6]*x[8],
    Curve over Rational Field defined by
    2048*x[1]^2 + 16*x[3]^2 + x[5]*x[8],
    512*x[1]*x[2] - x[5]*x[8],
    16*x[1]*x[4] + 4*x[7]^2 - x[8]^2,
    4*x[1]*x[5] - 4*x[3]*x[7] + x[4]*x[8],
    512*x[1]*x[6] - x[5]*x[7],
    x[3]*x[5] + 128*x[1]*x[7],
    x[4]*x[5] + 128*x[1]*x[8],
    128*x[2]^2 - 512*x[6]^2 + x[5]*x[8],
    4*x[2]*x[3] + x[7]*x[8],
    4*x[2]*x[4] + x[8]^2,
    x[2]*x[5] - x[4]*x[8],
    x[2]*x[7] - x[6]*x[8],
    -x[4]*x[5] - 128*x[6]*x[7] + 32*x[2]*x[8],
    4*x[3]*x[4] + x[5]*x[7],
    4*x[3]*x[6] + x[7]^2,
    -x[5]*x[6] + x[3]*x[8],
    4*x[4]^2 + x[5]*x[8],
    4*x[4]*x[6] + x[7]*x[8],
    -x[5]*x[6] + x[4]*x[7],
    x[5]^2 + 128*x[7]^2 - 32*x[8]^2,
    Curve over Rational Field defined by
    2048*x[1]^2 + 16*x[3]^2 - x[5]*x[8],
    512*x[1]*x[2] - x[5]*x[8],
    16*x[1]*x[4] + 4*x[7]^2 + x[8]^2,
    4*x[1]*x[5] - 4*x[3]*x[7] - x[4]*x[8],
    512*x[1]*x[6] - x[5]*x[7],
    x[3]*x[5] + 128*x[1]*x[7],
    x[4]*x[5] + 128*x[1]*x[8],
    128*x[2]^2 + 512*x[6]^2 - x[5]*x[8],
    4*x[2]*x[3] + x[7]*x[8],
    4*x[2]*x[4] + x[8]^2,
    x[2]*x[5] - x[4]*x[8],
    x[2]*x[7] - x[6]*x[8],
    x[4]*x[5] + 128*x[6]*x[7] + 32*x[2]*x[8],
    4*x[3]*x[4] + x[5]*x[7],
    4*x[3]*x[6] + x[7]^2,
    -x[5]*x[6] + x[3]*x[8],
    4*x[4]^2 + x[5]*x[8],
    4*x[4]*x[6] + x[7]*x[8],
    -x[5]*x[6] + x[4]*x[7],
    x[5]^2 + 128*x[7]^2 + 32*x[8]^2
]
*/
#EllipticCurve(Curve(Reduction((l[1]),5))); //4
#EllipticCurve(Curve(Reduction((l[2]),5))); //2

Qx<x> := PolynomialRing(Rationals());
E := EllipticCurve(x^3 - 2*x);
#EllipticCurve(Curve(Reduction(E,5))); //10

C1 := l[3];
//Search for degree 2 divisors on C1
P<[X]> := AmbientSpace(C1);
deg := 2;
tups := {};
d := Dimension(P);
ideals := {@@};
tupsOld := tups;
coeffs := [ 0, 1 ];
    tups   :=
        [ P![ tup[i] : i in [1..#tup] ] :
          tup in CartesianPower(coeffs, d+1) |
          not tup eq < 0 : i in [1..d+1] >  
          and not P![tup[i] : i in [1..#tup]] in tupsOld
        ];
    for c in tups  do  
            L := Scheme(P,&+[c[i]*X[i] : i in [1..Dimension(AmbientSpace(C1)) + 1]]);
              irr := IrreducibleComponents(L meet C1);
              degs := [Degree(ReducedSubscheme(i)) : i in irr];
              if  deg in degs then
        for cpt in irr do  
         if  Degree(ReducedSubscheme(cpt)) eq deg
             and not Basis(Ideal(ReducedSubscheme(cpt))) in ideals then      
               ideals := ideals join {@Basis(Ideal(ReducedSubscheme(cpt)))@}; 
        end if;
    end for;
    end if;
    end for;
    ideals;
