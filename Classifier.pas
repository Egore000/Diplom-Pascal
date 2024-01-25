unit Classifier;

interface
uses SysUtils,
    service;

procedure Classification(net: NETWORK; 
                        t: time_data;
                        phi, dot_phi: angle_data;
                        num: integer;
                        var classes: CLS);

implementation

procedure Classification(net: NETWORK;
                        t: time_data;
                        phi, dot_phi: angle_data;
                        num: integer;
                        var classes: CLS);
// Классификация резонанса
// Параметры:
// net - сетка графика
// t - массив времени
// phi, dot_phi - массивы углов и частот
// num - номер исследуемого объекта
// classes - выходной массив классов резонанса

// 0 - циркуляция
// 1 - либрация
// 2 - смешанный тип
var
  res, i, j: integer;
  zero_cols_counter, zero_rows_counter, zero_counter, count: integer;
  cols_count, rows_count: integer;
  inc_count, dec_count, perehod_count, class_: integer;

begin
  // Цикл по компонентам резонанса
  for res := res_start to res_end do
  begin  
    class_ := 0;
    zero_counter := 0;
    count := 0;
    perehod_count := 0;
    rows_count := 0;
    cols_count := 0;

    // Цикл по строчкам сетки
    for i := 1 to rows do
    begin
      zero_rows_counter := 0; {Счётчик нулевых ячеек в строке}
    
      // Цикл по столбцам
      for j := 1 to cols do
      begin
        count := count + net[res, i, j];

        if (net[res, i, j] = 0) then 
        begin
          inc(zero_counter); 
          inc(zero_rows_counter);
        end;
      end;

      if (zero_rows_counter = cols) then inc(rows_count);
    end;

    for j:= 1 to cols do
    begin
      zero_cols_counter := 0;

      for i := 1 to rows do
        if (net[res, i, j] = 0) then inc(zero_cols_counter);

      if (zero_cols_counter > 1) then inc(cols_count);
    end;

    if (zero_counter = 0) then class_ := 0;
    if (cols_count <> 0) then class_ := 2;
    if (rows_count <> 0) then class_ := 1;

    inc_count := 0;
    dec_count := 0;
    for i := 1 to count-1 do
    begin
      // Подсчёт убывющих точек
      if (phi[res, i] > phi[res, i+1]) or 
      (phi[res, i] < phi[res, i+1]/2) then inc(dec_count);
      // ((phi[res, i] < 10 * toRad) and (phi[res, i+1] > 350 * toRad)) then inc(dec_count);


      // Подсчёт возрастающих точек
      if (phi[res, i] < phi[res, i+1]) or 
      (phi[res, i]/2 > phi[res, i+1]) then inc(inc_count);
      // ((phi[res, i] > 350 * toRad) and (phi[res, i+1] < 10 * toRad)) then inc(inc_count);         
    
      // Подсчёт переходов частоты через 0
      if (dot_phi[res, i] * dot_phi[res, i+1]) < 0 then inc(perehod_count);
    end;

    // writeln('[COUNT]    ', count);
    // writeln('[INCREASE]   ', inc_count);
    // writeln('[DECREASE]   ', dec_count);

    // writeln('[ZERO FREQUENCE TRANSITION]   ', perehod_count);
    // if (perehod_count > count * coef) then class_ := 2;

    if (inc_count = count-1) then class_ := 0;
    if (dec_count = count-1) then class_ := 0;

    // writeln('[RESONANCE]', #9, res, #9, '[CLASS]', #9, class_);
    classes[res] := class_;
  end;
end;

begin {Main}
end.
