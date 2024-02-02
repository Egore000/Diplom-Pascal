unit TwoBody;

interface
uses service in '..\Tools\Service\service.pas',
    constants in '..\..\params\CONSTANT\constants.pas',
    config in '..\..\params\config\config.pas';

function f(x, M: extended): extended;
function df(x: extended): extended;
procedure TwoPoints(t, n, M0, w, Omega, i, a:extended; var X, V: mas);
procedure CoordsToElements(Coords, Velocities: mas; mu: extended; 
                          var a, e, i, Omega, w, M: extended);

implementation

function f(x, M: extended): extended;
// Уравнение Кеплера
begin
    f := x - ecc*sin(x) - M;
end; {f}

function df(x: extended): extended;
// Производная от уравнения Келпера
begin
    df := 1 - ecc*cos(x);
end; {df}

procedure TwoPoints(t, n, M0, w, Omega, i, a: extended; var X, V: mas);
// Задача двух тел
var M, E, E0, v0, u, dif, sum: extended;
    parametr, alpha0, beta0, gamma0, alpha, beta, gamma: extended;
    i1,i2,i3,i4: integer;
    Z1, Z2, Xi, prod: matrix;
    Orbit: mas;
begin
    Omega := Omega * toRad;
    M0 := M0 * toRad;
    i := i * toRad;
    w := w * toRad;

    for i1 := 1 to 3 do
    begin  
        X[i1] := 0;
        for i2 := 1 to 3 do
        begin
            Z1[i1, i2] := 0;
            Z2[i1, i2] := 0;
            Xi[i1, i2] := 0;
            prod[i1, i2] := 0;
        end;
    end;
    
    M := n*(t - t0) + M0;
    
    E0 := M;
    dif := 1;
    
    while (abs(dif) > eps) do
    begin
        E := E0 - f(E0,M)/df(E0);
        dif := E - E0;
        E0 := E;
    end;
    
    v0 := ArcTg((sqrt(1 - sqr(ecc)) * sin(E)), (cos(E) - ecc));
    
    Orbit[1] := a * (cos(E) - ecc);
    Orbit[2] := a * sqrt(1 - sqr(ecc)) * sin(E);
    Orbit[3] := 0;
    
    Z1[1,1] := cos(Omega);  Z2[1,1] := cos(w);
    Z1[1,2] := -sin(Omega); Z2[1,2] := -sin(w);
    Z1[2,1] := sin(Omega);  Z2[2,1] := sin(w);
    Z1[2,2] := cos(Omega);  Z2[2,2] := cos(w);
    Z1[3,3] := 1;             Z2[3,3] := 1; 
        
    Xi[1,1] := 1;
    Xi[2,2] := cos(i);
    Xi[2,3] := -sin(i);
    Xi[3,2] := sin(i);
    Xi[3,3] := cos(i);
    
    for i1 := 1 to 3 do
    for i2 := 1 to 3 do
    for i3 := 1 to 3 do
    for i4 := 1 to 3 do
        prod[i1, i4] := prod[i1, i4] + Z1[i1, i2]*Xi[i2, i3]*Z2[i3, i4];      
              
    for i1 := 1 to 3 do
    begin
        sum := 0;
        for i2 := 1 to 3 do  
            sum := sum + prod[i1,i2] * orbit[i2];
        X[i1] := sum;  
    end;
  
    parametr := a*(1 - sqr(ecc));
    
    u := v0 + w;
    
    alpha0 := cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i);
    alpha := -sin(u)*cos(Omega) - cos(u)*sin(Omega)*cos(i);
    
    beta0 := cos(u)*sin(Omega) + sin(u)*cos(Omega)*cos(i);
    beta := -sin(u)*sin(Omega) + cos(u)*cos(Omega)*cos(i);
    
    gamma0 := sin(u)*sin(i);
    gamma := cos(u)*sin(i);
    
    V[1] := sqrt(mu/parametr)*(ecc*sin(v0)*alpha0 + (1 + ecc*cos(v0))*alpha);
    V[2] := sqrt(mu/parametr)*(ecc*sin(v0)*beta0 + (1 + ecc*cos(v0))*beta);
    V[3] := sqrt(mu/parametr)*(ecc*sin(v0)*gamma0 + (1 + ecc*cos(v0))*gamma);
end;



procedure CoordsToElements(Coords, Velocities: mas; mu: extended;
                           var a, e, i, Omega, w, M: extended);
var x, y, z, Vx, Vy, Vz: extended;
    r, V2, h, c1, c2, c3, l1, l2, l3, E0: extended;
    c, l: extended;
begin
    x := Coords[1];
    y := Coords[2];
    z := Coords[3];
    Vx := Velocities[1];
    Vy := Velocities[2];
    Vz := Velocities[3];

    r := sqrt(x*x + y*y + z*z);
    V2 := Vx*Vx + Vy*Vy + Vz*Vz;
    
    h := V2/2 - mu/r;
    
    c1 := y*Vz - z*Vy;
    c2 := z*Vx - x*Vz;
    c3 := x*Vy - y*Vx;
    
    l1 := -mu*x/r + Vy*c3 - Vz*c2;
    l2 := -mu*y/r + Vz*c1 - Vx*c3;
    l3 := -mu*z/r + Vx*c2 - Vy*c1;
    
    c := sqrt(c1*c1 + c2*c2 + c3*c3);
    l := sqrt(l1*l1 + l2*l2 + l3*l3);
    
    a := -mu/(2*h);
    e := l/mu;
    i := ArcTg(sqrt(1 - (c3*c3/(c*c))), c3/c);

    if (i = 0) then i := i + 1e-12;

    Omega := Arctg2(c1/(c*sin(i)), -c2/(c*sin(i)));
    w := Arctg2(l3/(l*sin(i)), l1*cos(Omega)/l + l2*sin(Omega)/l);
    
    E0 := Arctg2((x*Vx + y*Vy + z*Vz)/(e * sqrt(mu * a)), (1 - r/a)/e);

    M := E0 - e*sin(E0);
end;  

begin {Main}
end.