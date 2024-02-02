uses SysUtils, 
    Classifier in 'MODULES\Classifier\Classifier.pas',
    ResonanceUnit in 'MODULES\Resonance\Resonance.pas',
    readfond in 'MODULES\Tools\ReadFond\readfond.pas',     
    TwoBody in 'MODULES\TwoBody\TwoBody.pas', 
    service in 'MODULES\Tools\Service\service.pas',
    constants in 'params\CONSTANT\constants.pas',
    config in 'params\config\config.pas';

const 
      // Пути к файлам и директориям
      TARGER_FOLDER = 'Без светового давления';
      // TARGER_FOLDER = 'Со световым давлением';

      PATH_DATA = '..\Исходные данные\' + TARGER_FOLDER + '\'; // Путь к папке с исходными данными
      PATH_OUT = 'Initial_coords\ecc.DAT';
      
var coords, velocities: mas; // Массивы координат и скоростей
    
    a, e, max_e, i, Omega, w, M, megno, mean_megno: extended; 
    tm, time, day: extended;
    year, month, num, number, x: integer;

    data, eccentr: text; // Файлы
                                                // data - файл с исходными данными
                                                // eccentr - файл для записи эксцентриситетов
    ss: string[7]; // Служебная строка

    folder: integer; // Папка с исходными файлами
    file_num: string;


begin {Main}
    {$WARNINGS-}
    assign(eccentr, PATH_OUT);
    rewrite(eccentr);

    {Цикл по папкам}
    for folder := start_folder to finish_folder do
    { Цикл по файлам в папке folder }
        for number := start to finish do
        begin
            if (number < 10) then file_num := '000' + inttostr(number);
            if (number >= 10) and (number < 100) then file_num := '00' + inttostr(number);
            if (number >= 100) and (number < 1000) then file_num := '0' + inttostr(number);
            if (number >= 1000) then file_num := inttostr(number);

            if FileExists(PATH_DATA + inttostr(folder) + '\EPH_' + file_num + '.DAT') then
            begin
                assign(data, PATH_DATA + inttostr(folder) + '\EPH_' + file_num + '.DAT');
                reset(data);
                writeln('[FILE]', #9, number);
            end
            else
            begin
                writeln('Finished!');
                halt;
            end;

            max_e := 0;
            while not eof(data) do
            begin
                {Считывание данных из файла}
                readln(data, tm, time, ss, year, month, day);
                readln(data, x, coords[1], coords[2], coords[3], megno);
                readln(data, velocities[1], velocities[2], velocities[3], mean_megno);

                {Расчёт элментов орбиты}
                CoordsToElements(coords, velocities, mu, a, e, i, Omega, w, M);

                if (e > max_e) then max_e := e;  
            end;  
            writeln(eccentr, max_e);

            close(data);
        end;
    close(eccentr);
end.