unit config;

interface
const 
    ORBITAL = true; // Исследование орбитального резонанса
    SECONDARY = false; // Исследование вторичных резонансов

    // Порядок резонанса (u:v)
    u = 1;
    v = 2;

    DEBUG = false; // Дебаггинг

    start_folder = 2; // Начальная папка
    finish_folder = 2; // Конечная папка

    start = 1000; // Начальный файл
    finish = 1100; // Конечный файл

    res_start = 1; // Начальная компонента резонанса
    res_end = 5; // Конечная компонента резонанса
    
    rows = 36; // Количество строк в сетке
    cols = 20; // Количество столбцов в сетке
    row_step = 360 / rows; // Шаг по строкам
    col_step = 100 / cols; // Шаг по столбцам

    libration_rows = 24; // Количество строк в разбиении для определения либрации
    libration_step = 360 / libration_rows; // Шаг по строкам при определении либрации

    delimiter = #9; // Разделитель в выходных файлах
implementation
begin
end.