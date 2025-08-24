/****************************************************************************** 
Here is a summary of the argument.

Let C be the modular curve with lmfdb label 16.96.5.bq.1.
Using a singular model we find a rank 1 quotient. This is positive rank bielliptic. 
******************************************************************************/

P<x,y,z>:=ProjectiveSpace(Rationals(),2);

C := Curve(P,[x^4*y^4 - 4*x^4*y^2*z^2 + 4*x^4*z^4 + 2*x^2*y^6 - 44*x^2*y^4*z^2 + 116*x^2*y^2*z^4 - 72*x^2*z^6 + y^8 + 24*y^6*z^2 + 116*y^4*z^4 - 336*y^2*z^6 + 196*z^8]);
S := AutomorphismGroup(C); 

auts := [];
Stemp := Automorphisms(C);
for s in Stemp do
auts := Append(auts, S!s);
end for;
assert #auts eq #S;

//There are three genus one quotients by an involution
l := []; //list of genus 1 quotients by involutions
m:= []; //in case Magma complains that genus 1 curves and elliptic curves can't be in the same list.
for g in auts do
if Order(g) eq 2 then
AG := AutomorphismGroup(C,[g]);
CG,prj := CurveQuotient(AG);
if Genus(CG) eq 1 then
try
l := Append(l,CG);
catch e
m := Append(m,CG);
end try;
end if;

end if;

end for;

C1:=l[1];

// Constructing a simpler model

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
Pt:=PointSearch(C2,100)[1];
E:=EllipticCurve(C2,Pt);
Rank(E); // 1 true
