# Установка параметров
options <- list(encoding = "UTF-8")

# Загрузка данных
data <- read.csv("ilr_table.csv")

# Выполнение теста Уилкоксона для каждого столбца
wilcox_test_results <- lapply(data, function(x) wilcox.test(x, alternative = "two.sided"))

# Фильтрация признаков
filtered_features <- names(which(sapply(wilcox_test_results, function(x) x$p.value > 0.5)))

# Подготовка данных для Random Forest
filtered_data <- data[, filtered_features]

# Загрузка библиотеки Random Forest
library(randomForest)

# Открываем файл для записи
sink("random_forest_output.txt")

# Обучение модели Random Forest
model <- randomForest(filtered_data, ntree = 100)

# Оценка модели
print(model)

# Вывод важности признаков
print(importance(model))

# Закрываем файл
sink()

