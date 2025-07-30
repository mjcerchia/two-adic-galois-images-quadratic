Is this actually bielliptic? There's a rank 1 quotient of a singular model 

/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 37.114.4.a.1. 
We find there are four genus one quotients by an involution. 

******************************************************************************/
P<x,y,z> := ProjectiveSpace(Rationals(),2);
C := Curve(P,[778*x^6 - 2294*x^5*y + 547*x^5*z + 814*x^4*y^2 + 3367*x^4*y*z - 3513*x^4*z^2 - 962*x^3*y^3 - 14023*x^3*y^2*z + 18907*x^3*y*z^2 - 8863*x^3*z^3 - 259*x^2*y^3*z - 24235*x^2*y^2*z^2 + 31413*x^2*y*z^3 - 8564*x^2*z^4 - 11063*x*y^3*z^2 - 24013*x*y^2*z^3 + 19832*x*y*z^4 - 4944*x*z^5 - 3145*y^3*z^3 - 4292*y^2*z^4 + 5328*y*z^5 - 1664*z^6]);

S := AutomorphismGroup(C); 
auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
#auts eq #S;

//There are four genus one quotients by an involution
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
C1:=l[3];

//This is quite complicated, so we search for divisors to construct a simpler model.
P<[X]> := AmbientSpace(C1);
deg := 3;
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

//Using a degree three divisor found above, we build a new model using Riemann-Roch.
gens := ideals[1];
P7 := AmbientSpace(C1);
R  := CoordinateRing(P7);
I := ideal< R | gens >;
Pl := Place(C1, I);
D  := Divisor(Pl);
a,b := RiemannRochSpace(D);
x1 := b(a.1); 
x2 := b(a.2);
x3:=b(a.3); 
F := BaseField(C1);
P := PolynomialRing( F, [1,1,1], "grevlexw", [1,1,1] );
vars := ["x", "y","z"]; 
P := ProjectiveSpace(P);
AssignNames(~P, vars);
phi := map<C1 -> P | [x1,x2,x3]>;
C2:=Image(phi); 

E := EllipticCurve(C2);
    
