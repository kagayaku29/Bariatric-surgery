# Установка необходимых пакетов, если они еще не установлены
install.packages("vegan")

# Загрузка библиотеки vegan для анализа микробиоты
library(vegan)

# Загрузка данных из файла CSV
data <- read.csv("ilr_table.csv")

# Предположим, что первый столбец - это метки групп
group_labels <- data[, 1]

# Предположим, что остальные столбцы содержат данные
# Преобразование данных в матрицу расстояний
distance_matrix <- vegdist(data[, -1])

# Выполнение анализа перманова
result <- adonis(distance_matrix ~ group_labels)

# Открываем файл для записи
sink("adonis_output.txt")

# Выводим результаты анализа
print(result)

# Закрываем файл
sink()





