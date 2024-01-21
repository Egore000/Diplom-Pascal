unit service;

interface
uses readfond;

const a_e = 149597870.691; //[km]
      pi = 3.1415926535897932;
      mu = 3.986004418e+5; // Гравитационная постоянная для спутника
      Gmu = 1.32712442099e+11; // Гравитационная постоянная для Солнца
      toRad = pi/180; // Перевод в радианы
      toDeg = 180/pi; // Перевод в градусы
      
      rows = 6; // Количество строк в сетке
      cols = 60; // Количество столбцов в сетке
      row_step = 360 / rows; // Шаг по строкам
      col_step = 100 / cols; // Шаг по столбцам

      eps = 1e-12; // Точность вычисления аномалии в задаче двух тел
      t0 = 0; // Начальная эпоха

      ecc = 1e-3; // Эксцентриситет орбиты

type
    matrix = array[1..3,1..3] of extended; // Матрицы поворота в задаче двух тел (модуль TwoBody.pas)
    mas = array[1..3] of extended; // Массив скоростей или координат
    CLS = array[1..5] of integer; // Массив для классификации компонент резонансного аргумента
    arr = array[1..5] of extended; // Массив резонансных аргументов или частот
    angle_data = array[1..5, 1..2000] of extended; // Матрица с полным набором резонансных углов 
    time_data = array[1..2000] of extended; // Вектор с моментами времени
    NETWORK = array[1..5, 1..rows + 1, 1..cols + 1] of integer; // Сетка разбиения данных для классификации

procedure Create_File(var f: text; path: string);
procedure WriteToFile(var f: text; time: extended; angles, freq: arr);
procedure WriteClassification(var f: text; folder, number: integer; classes, classes2, classes3: CLS);
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

implementation

procedure Create_File(var f: text;
                    path: string);
// Создание файла и запись в него заголовка
begin
    assign(f, path);
    rewrite(f);
    writeln(f, 't', #9, 'F1', #9, 'F2', #9, 'F3', #9, 'F4', #9, 'F5', #9, 'dF1', #9, 'dF2', #9, 'dF3', #9, 'dF4', #9, 'dF5');
end;


procedure WriteToFile(var f: text;
                      time: extended;
                      angles, freq: arr);
// Запись данных о резонансах в файл f
begin
    writeln(f, time/(86400 * 365), ' ',
              angles[1] * toDeg, ' ',
              angles[2] * toDeg, ' ',
              angles[3] * toDeg, ' ',
              angles[4] * toDeg, ' ',
              angles[5] * toDeg, ' ',
              freq[1], ' ', 
              freq[2], ' ',
              freq[3], ' ',
              freq[4], ' ',
              freq[5]);
end;


procedure WriteClassification(var f: text;
                            folder, number: integer;
                            classes, classes2, classes3: CLS);
// Запись в файл с классификацией
const delimiter = ';';
begin
    writeln(f, folder, delimiter,
            number, delimiter,
            classes[1], delimiter,
            classes[2], delimiter,
            classes[3], delimiter,
            classes[4], delimiter,
            classes[5], delimiter,
            classes2[1], delimiter,
            classes2[2], delimiter,
            classes2[3], delimiter,
            classes2[4], delimiter,
            classes2[5], delimiter,
            classes3[1], delimiter,
            classes3[2], delimiter,
            classes3[3], delimiter,
            classes3[4], delimiter,
            classes3[5]);
end;


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
end;


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
end; //perehod


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
end;

begin
end.