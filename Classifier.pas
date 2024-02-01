unit Classifier;

interface
uses SysUtils,
    service;

procedure Classification(net: NETWORK; 
                        flag: FLAGS;
                        t: time_data;
                        phi, dot_phi: angle_data;
                        var classes: CLS);

implementation

procedure Classification(net: NETWORK;
                        flag: FLAGS;
                        t: time_data;
                        phi, dot_phi: angle_data;
                        var classes: CLS);
// Классификация резонанса
// Параметры:
// net - сетка графика
// flag - сетка полос для выявления либрации
// t - массив времени
// phi, dot_phi - массивы углов и частот
// classes - выходной массив классов резонанса

// 0 - циркуляция
// 1 - смешанный тип
// 2 - либрация
var
    res, i, j: integer;
    zero_counter, count, libration: integer;
    perehod_count: integer;
    increase, decrease: array[1..5] of integer;

begin
    for res := res_start to res_end do
    begin
        increase[res] := 0;
        decrease[res] := 0;
    end;
    

    // Цикл по компонентам резонанса
    for res := res_start to res_end do
    begin  
        zero_counter := 0;
        count := 0;
        perehod_count := 0;

        // Цикл по строчкам сетки
        for i := 1 to rows do
            // Цикл по столбцам
            for j := 1 to cols do
            begin
                count := count + net[res, i, j];

                if (net[res, i, j] = 0) then inc(zero_counter); 
            end; {for j}


        // if (zero_counter < 3) then classes[res] := 0;
        // if (cols_count <> 0) then classes[res] := 1;
        // if (rows_count <> 0) then classes[res] := 2;

        for i := 1 to count-1 do
        begin
            // Подсчёт убывющих точек
            if ((decrease[res] = 0) and (phi[res, i] > phi[res, i+1])) or 
            ((phi[res, i]/2 < phi[res, i+1]) and (phi[res, i] < 100)) then 
            else
                decrease[res] := 1;

            // Подсчёт возрастающих точек
            if ((increase[res] = 0) and (phi[res, i] < phi[res, i+1])) or 
            ((phi[res, i]/2 > phi[res, i+1]) and (phi[res, i] > 260)) then
            else
                increase[res] := 1;   
    
            // Подсчёт переходов частоты через 0
            if (dot_phi[res, i] * dot_phi[res, i+1]) < 0 then inc(perehod_count);
        end; {for i}
        // if (perehod_count > count * coef) then class_ := 2;
    end; {for res}



    for res := res_start to res_end do
    begin
        libration := 0;

        for i := 1 to libration_rows do
            if (flag[res, i] = 0) then inc(libration);

        if (libration > 0) and (increase[res] <> 0) and (decrease[res] <> 0) then
            classes[res] := 2
        else
            if (increase[res] = 0) or (decrease[res] = 0) or (zero_counter < 3) then
                classes[res] := 0
            else
                classes[res] := 1;

        if DEBUG then
        begin
            writeln('[RESONANCE]', #9, res, #9, '[CLASS]', #9, classes[res]);
            writeln('[COUNT]  ', count);
            writeln('[ZEROS]  ',  zero_counter);
            writeln('[ZERO FREQUENCE TRANSITION]   ', perehod_count);
            writeln;
        end; {if}
    end; {for res}


    if DEBUG then
    begin
        write('[INC]  ');
        for i := res_start to res_end do
            write(increase[i], #9);
        writeln;

        write('[DEC]  ');
        for i := res_start to res_end do
            write(decrease[i], #9);
        writeln;
    end;
end;

begin {Main}

end.
