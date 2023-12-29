unit service;

interface

procedure Reduce(A:extended;var B:extended);
function Sign(n:extended):Shortint;
function ArcTg(x,y: extended):extended;
function ArcSin(x: extended):extended;
function date_jd(year,month:integer;day:extended):extended;
function sid2000(jd:extended):extended; {v radianah}

implementation

procedure Reduce(A:extended;var B:extended);
const PI2 = 2*pi;
begin
    B := A - trunc(A/PI2)*PI2;
    if (A < 0) then B := B + PI2;
end; //Reduce

function Sign(n:extended):Shortint;
begin
    if (n > 0) then Sign := 1
    else 
        if (n < 0) then Sign := -1
        else Sign := 0;
end; //Sign

function ArcTg(x, y: extended):extended;
var a: extended;
begin
    if (abs(y) < 1e-18) then ArcTg := sign(x)*0.5*Pi
    else
    begin
        a := ArcTan(x/y);
        if (y < 0) then a := a+Pi;
        if (a < 0) then a := a+2*Pi;
        ArcTg := a;
    end;
end; //ArcTg

function ArcSin(x: extended):extended;
{ arcsin in [Pi/2,-Pi/2] }
begin
    if (x < 1.0) and (x > -1.0) then
        Arcsin := ArcTan(x/sqrt(1 - sqr(x)))
    else
        if ((abs(x) - 1.0) < 1e-18) then ArcSin := Sign(x)*Pi/2
        else
        begin 
            Writeln('Error: Argument of ArcSin >1.0 or <-1.0', x); 
            Halt; 
        end;
end; {Arcsin}

function date_jd(year, month: integer;
                day: extended): extended;
label 1;
const han=100;
var m1,d,date,jd: extended;
    i,me,ja,jb: integer;
begin
    date := year + month/100 + day/1e4;
    i := trunc(date);
    m1 := (date - i)*han;
    me := trunc(m1);
    d := (m1 - me)*han;
    if (me > 2) then goto 1;
    i := i - 1;
    me := me + 12;
1:  jd := trunc(365.25*i) + trunc(30.6001*(me+1)) + d + 1720994.5;
    if (date<1582.1015) then
    begin
        date_jd := jd;
        exit;
    end;
    ja := trunc(i/100.0);
    jb := 2 - ja + trunc(ja/4);
    jd := jd+jb;
    date_jd := jd;
end; //date_jd()

function sid2000(jd: extended): extended; {v radianah}
const jd2000 = 2451545;
      jdyear = 36525;
var m,mm,d,t,s,sr: extended;
begin
    m := frac(jd) - 0.5;
    d := jd - m - jd2000;

    t := (d + m)/jdyear;
    mm := m*86400;
    s := (24110.54841+mm+236.555367908*(d+m)+(0.093104*t-6.21E-6*sqr(t))*t)/86400*2*pi;
    sid2000 := s;
end; // sid2000()

begin
end.