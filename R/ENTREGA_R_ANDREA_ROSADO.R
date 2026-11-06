##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

# EJERCICIO 1: Cargar y explorar los datos

datos <- read.table("archivo-entrada.txt",
                    header = TRUE,
                    sep = "\t",
                    check.names = FALSE)

head(datos)
tail(datos)
dim(datos)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

# EJERCICIO 2: Medias y desviaciones estándar

# Seleccionar las tres réplicas de cada grupo
control <- datos[, 2:4]
condicion <- datos[, 5:7]

# Calcular las medias
datos$Media_Control <- apply(control, 1, mean)
datos$Media_Condicion <- apply(condicion, 1, mean)

# Calcular las desviaciones estándar
datos$DE_Control <- apply(control, 1, sd)
datos$DE_Condicion <- apply(condicion, 1, sd)

# Comprobar los resultados
head(datos)
dim(datos)

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

# EJERCICIO 3: Histogramas

# Dos gráficos en una sola imagen
par(mfrow = c(1, 2))

# Histograma del control
hist(log2(datos$Media_Control + 1),
     col = "skyblue",
     main = "Control",
     xlab = "log2(TPM + 1)",
     ylab = "Número de genes")

# Histograma de la condición
hist(log2(datos$Media_Condicion + 1),
     col = "salmon",
     main = "Condición",
     xlab = "log2(TPM + 1)",
     ylab = "Número de genes")

# Volver a la disposición habitual
par(mfrow = c(1, 1))

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

# EJERCICIO 4: Boxplots

# Comparar las medias de control y condición
boxplot(log2(datos$Media_Control + 1),
        log2(datos$Media_Condicion + 1),
        names = c("Control", "Condición"),
        col = c("skyblue", "salmon"),
        main = "Expresión génica: Control y Condición",
        ylab = "log2(TPM + 1)")

# Añadir una leyenda en la esquina superior izquierda
legend("topleft",
       legend = c("Control", "Condición"),
       fill = c("skyblue", "salmon"))

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

# EJERCICIO 5: Normalidad y comparación estadística

# Comprobar la normalidad del control
shapiro.test(log2(datos$Media_Control + 1))

# Comprobar la normalidad de la condición
shapiro.test(log2(datos$Media_Condicion + 1))

# Comparar control y condición para cada gen
# Guardar los valores p en una nueva columna

datos$p_value <- rep(NA_real_, nrow(datos))

for (i in 1:nrow(datos)) {
  
  # Seleccionar las tres réplicas del gen
  x <- as.numeric(control[i, ])
  y <- as.numeric(condicion[i, ])
  
  # Aplicar la prueba t si es posible
  prueba <- try(t.test(x, y), silent = TRUE)
  
  if (!inherits(prueba, "try-error")) {
    datos$p_value[i] <- prueba$p.value
  }
}

# Comprobar los resultados
head(datos)
dim(datos)
summary(datos$p_value)

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

# EJERCICIO 6: Genes diferencialmente expresados

# Seleccionar los genes con p-value menor de 0.05
genes_DE <- datos[!is.na(datos$p_value) &
                    datos$p_value < 0.05, ]

# Mostrar los genes seleccionados
genes_DE

# Contar los genes diferencialmente expresados
nrow(genes_DE)

# Genes sobreexpresados en la condición
sobreexpresados <- genes_DE[
  genes_DE$Media_Condicion > genes_DE$Media_Control, ]

nrow(sobreexpresados)

# Genes reprimidos en la condición
reprimidos <- genes_DE[
  genes_DE$Media_Condicion < genes_DE$Media_Control, ]

nrow(reprimidos)

# Interpretación de resultados:
# Se identifican 362 genes con p-value < 0.05.
# De ellos, 50 están sobreexpresados y 312 reprimidos
# en la condición experimental respecto al control.
# Se utiliza el criterio de significación sin ajustar
# solicitado en el enunciado.
#
# La normalidad se evaluó sobre las medias de los genes,
# no sobre las réplicas individuales de cada gen.
# Con tres réplicas por grupo, la prueba t de Welch
# debe interpretarse con cautela.

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# EJERCICIO 7: Heatmap

# Seleccionar las seis columnas de expresión
matriz <- as.matrix(genes_DE[, 2:7])

# Asignar el nombre de cada gen a su fila
rownames(matriz) <- make.unique(genes_DE[, 1])

# Transformación logarítmica
matriz_log <- log2(matriz + 1)

# Representar el heatmap
heatmap(matriz_log,
        col = colorRampPalette(c("skyblue",
                                 "white",
                                 "salmon"))(100),
        scale = "row",
        Colv = NA,
        labRow = NA,
        labCol = c("Control 1", "Control 2", "Control 3",
                   "Condición 1", "Condición 2", "Condición 3"),
        main = "Genes diferencialmente expresados",
        margins = c(9, 5))

#TRABAJO FINALIZADO