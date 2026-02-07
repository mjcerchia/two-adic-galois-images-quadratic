
//We observe an automorphism for 16.96.3.ds.1 that gives a rank 1 elliptic curve as a quotient. This is positive rank bielliptic.

P<x,y,z,w,t,u> := ProjectiveSpace(Rationals(),5);
C := Curve(P,[y^2 + 2*y*u + w*t + t^2, 
2*y^2 + w^2 + t^2,
-2*x*u + z^2, 
4*x*y + y*u + w*t,
4*x*w + 4*x*t - y*w - w*u + t*u, 
4*x*w - 4*x*t - y*t + w*u + t*u, 
32*x^2 - y^2 + 2*u^2]);

/*********************************Optional code to verify that models taken from LMFDB match with that taken from Zywina's repo

for tuple in data211 do;
          if tuple[1] eq "16.96.3.ds.1" then
                      level:=Split(tuple[1],".")[1];
                      level:=StringToInteger(level);
                      GL2:=GL(2,Integers(level));
                      G:=sub<GL2|tuple[4]>;
                      Gt:=sub<GL2|[Transpose(GL2!g):g in Generators(G)]>;
                      X:=CreateModularCurveRec(Gt);
                      XG:=FindModelOfXG(X);
                      D := Curve(ProjectiveSpace(Rationals(), Rank(Parent((XG`psi)[1]))-1),XG`psi);
           end if;
end for;

P<x,y,z,w,t,u> := ProjectiveSpace(Rationals(),5);
C := Curve(P,[y^2 + 2*y*u + w*t + t^2, 
2*y^2 + w^2 + t^2,
-2*x*u + z^2, 
4*x*y + y*u + w*t,
4*x*w + 4*x*t - y*w - w*u + t*u, 
4*x*w - 4*x*t - y*t + w*u + t*u, 
32*x^2 - y^2 + 2*u^2]);

assert IsIsomorphic(C,D);

**************************************************************/


s:=iso<C->C|[u/4,y,z,-t,-w,4*x],[1/4*u,y,z,-t,-w,4*x]>;

AG := AutomorphismGroup(C,[s]);
CG,prj := CurveQuotient(AG);

Pts:=PointSearch(CG,10000);

E:=EllipticCurve(CG,Pts[1]);
Rank(E); // 1 true
