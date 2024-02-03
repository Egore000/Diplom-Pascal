unit service;

interface
uses SysUtils,
    readfond in '..\ReadFond\readfond.pas',
    constants in '..\..\params\CONSTANT\constants.pas',
    config in '..\..\params\config\config.pas';

type
    matrix = array[1..3,1..3] of extended; // Матрицы поворота в задаче двух тел (модуль TwoBody.pas)
    mas = array[1..3] of extended; // Массив скоростей или координат
    CLS = array[1..5] of integer; // Массив для классификации компонент резонансного аргумента
    arr = array[1..5] of extended; // Массив резонансных аргументов или частот
    angle_data = array[res_start..res_end, 1..2000] of extended; // Матрица с полным набором резонансных углов 
    time_data = array[1..2000] of extended; // Вектор с моментами времени
    NETWORK = array[res_start..res_end, 1..rows + 1, 1..cols + 1] of integer; // Сетка разбиения данных для классификации
    FLAGS = array[res_start..res_end, 1..libration_rows] of integer;

procedure Create_File(var f: text; path: string);
procedure WriteToFile(var f: text; time: extended; angles, freq: arr);
procedure WriteClassification(var f: text; folder, number: integer; a0, i0, megno: extended; classes, classes2, classes3: CLS);
procedure WriteHeader(var f: text; min, max: integer);
procedure fond405(jd:extended; var xm_,xs_,vm_,vs_:mas);  
function sid2000(jd:extended):extended; {v radianah}
function date_jd(year,month:integer;day:extended):extended;
procedure Reduce(A:extended;var B:extended);
function Sign(n:extended):Shortint;
function ArcTg(x,y: extended):extended;
function Arctg2(x, y: extended):extended;
function ArcSin(x: extended):extended;
procedure perehod(HH: extended; xx: mas; var lambda, phi: extended);
procedure OutNET(net: NETWORK);
procedure FillZero(var net, net2, net3: NETWORK;
                    var flag, flag2, flag3: FLAGS;
                    var t: time_data;
                    var phi, phi2, phi3: angle_data;
                    var dot_phi, dot_phi2, dot_phi3: angle_data);

implementation

procedure Create_File(var f: text;
                    path: string);
// Создание файла и запись в него заголовка
var i: integer;
begin
    assign(f, path);
    rewrite(f);
    write(f, 't', delimiter);
    if ORBITAL then
        for i := res_start to res_end do write(f, 'F', i, delimiter);
    if SECONDARY then
        for i := res_start to res_end do write(f, 'dF', i, delimiter);
    writeln(f);
end; {Create_File}


procedure WriteToFile(var f: text;
                      time: extended;
                      angles, freq: arr);
// Запись данных о резонансах в файл f
var i: integer;
begin
    write(f, time/(86400 * 365), delimiter);
    for i := res_start to res_end do write(f, angles[i] * toDeg, delimiter);
    for i := res_start to res_end do write(f, freq[i], delimiter);
    writeln(f);
end; {WriteToFile}


procedure WriteClassification(var f: text;
                            folder, number: integer;
                            a0, i0, megno: extended; 
                            classes, classes2, classes3: CLS);
// Запись в файл с классификацией
var i: integer;
begin
    if INITIAL_ELEMENTS then
        write(f, folder, delimiter, 
                number, delimiter, 
                a0, delimiter,
                i0, delimiter,
                megno, delimiter)
    else
        write(f, folder, delimiter, number, delimiter);
    
    
    if ORBITAL then
        for i := res_start to res_end do write(f, classes[i], delimiter);
    if SECONDARY then
    begin
        for i := res_start to res_end do write(f, classes2[i], delimiter);
        for i := res_start to res_end do write(f, classes3[i], delimiter);
    end;
    writeln(f);
end; {WriteClassification}



procedure WriteHeader(var f: text;
                    min, max: integer);
var i: integer;
begin
    if INITIAL_ELEMENTS then
        write(f, 'folder', delimiter, 
                'file', delimiter,
                'a, km', delimiter,
                'i, grad', delimiter,
                'MEGNO', delimiter)
    else        
        write(f, 'folder', delimiter, 
                'file', delimiter);

    if ORBITAL then
        for i := min to max do write(f, 'F', i, delimiter);
    if SECONDARY then
    begin
        for i := min to max do write(f, 'dF' + inttostr(i) + '(+)', delimiter);
        for i := min to max do write(f, 'dF' + inttostr(i) + '(-)', delimiter);
    end;
    writeln(f);
end; {WriteHeader}



procedure fond405(jd:extended; var xm_,xs_,vm_,vs_:mas);  
var x_planet:masc;
begin
    read405(0, 1, jd, x_planet);    
    // процедура выдает координаты и скорости Луны и Солнца в геоцентрической 
    // экваториальной — км и км/с

    // переход в геоцентрическую экваториальную — в км

    // координаты Луны и Солнца
    xm_[1] := (x_planet[55] - x_planet[13]) * a_e;  //
    xm_[2] := (x_planet[56] - x_planet[14]) * a_e;  //
    xm_[3] := (x_planet[57] - x_planet[15]) * a_e;  //   {xle=xl-xe}
        
    xs_[1] := -x_planet[13] * a_e;  //
    xs_[2] := -x_planet[14] * a_e; //
    xs_[3]:= - x_planet[15] * a_e; //    {xs=-xe }

    //скорости Луны и Солнца
    vm_[1] := (x_planet[58] - x_planet[16]) * a_e/86400;
    vm_[2] := (x_planet[59] - x_planet[17]) * a_e/86400;
    vm_[3] := (x_planet[60] - x_planet[18])*a_e/86400;     {vle=vl-ve}

    vs_[1] := -x_planet[16] * a_e/86400;
    vs_[2] := -x_planet[17] * a_e/86400;
    vs_[3] := -x_planet[18] * a_e/86400;     {vs=-ve}
end; {fond405}


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
end; {sid2000}


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
end; {date_jd()}


procedure Reduce(A:extended;var B:extended);
const PI2 = 2*pi;
begin
    B := A - trunc(A/PI2)*PI2;
    if (A < 0) then B := B + PI2;
end; {Reduce}


function Sign(n:extended):Shortint;
begin
    if (n > 0) then Sign := 1
    else 
        if (n < 0) then Sign := -1
        else Sign := 0;
end; {Sign}


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
end; {ArcTg}


function Arctg2(x, y: extended):extended;
var a: extended;
begin
    if (abs(y)<1e-18) then
        if (x>0) then  
            Arctg2 := sign(x)*0.5*pi
        else
            Arctg2 := -sign(x)*0.5*pi
    else
    begin
        a := ArcTan(x/y);
        if (y>0) then
        
        else if (x>=0) then 
                a := a + pi
            else
                a := a - pi;
    
        Arctg2 := a;
    end;
end; {Arctg2}


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


procedure perehod(HH: extended; 
                  xx: mas; 
                  var lambda, phi: extended);
// Переход во вращающуюся СК 
var A: array[1..3,1..3] of extended;
    i, j: integer;
    s, r, argum: extended;
    yy: mas;
begin

    A[1,1] := cos(HH);     A[1,2] := sin(HH);    A[1,3] := 0;
    A[2,1] := -sin(HH);    A[2,2] := cos(HH);    A[2,3] := 0;
    A[3,1] := 0;           A[3,2] := 0;          A[3,3] := 1;

    for i:=1 to 3 do
    begin
        s := 0;
        for j:=1 to 3 do
            s := s + A[i,j]*xx[j];
        yy[i] := s;
    end;
    r := sqrt(yy[1]*yy[1] + yy[2]*yy[2] + yy[3]*yy[3]);

    lambda := arctg(yy[2], yy[1]);

    argum := yy[3]/r;
    phi := arcsin(argum);
end; {perehod}


procedure OutNET(net: NETWORK);
// Вывод матриц
var 
  num, row, col: integer;
begin
    for num := 1 to 5 do
    begin
        writeln('num = ', num);
        for row := 1 to rows do
        begin
            for col := 1 to cols do
            begin
                write(net[num, row, col], #9);
            end;
            writeln;
        end;
        writeln;
    end;
end; {OutNET}


procedure FillZero(var net, net2, net3: NETWORK;
                    var flag, flag2, flag3: FLAGS;
                    var t: time_data;
                    var phi, phi2, phi3: angle_data;
                    var dot_phi, dot_phi2, dot_phi3: angle_data);
var num, row, col: integer;
begin
    for num := res_start to res_end do 
    begin
        for row := 1 to rows do
            for col := 1 to cols do
            begin
            net[num, row, col] := 0;
            net2[num, row, col] := 0;
            net3[num, row, col] := 0;
            end;
        
        for row := 1 to 2000 do
        begin
            t[row] := 0;

            phi[num, row] := 0;
            phi2[num, row] := 0;
            phi3[num, row] := 0;

            dot_phi[num, row] := 0;
            dot_phi2[num, row] := 0;
            dot_phi3[num, row] := 0;
        end;

        for row := 1 to libration_rows do
        begin
            flag[num, row] := 0;
            flag2[num, row] := 0;
            flag3[num, row] := 0;
        end;
    end;
end; {FillZero}

begin

end.