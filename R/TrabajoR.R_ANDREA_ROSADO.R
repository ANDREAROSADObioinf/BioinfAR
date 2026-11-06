##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)
# EJERCICIO 1: Carga y exploración de datos
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

# Seleccionar las réplicas
control <- datos[, c("C2-TPM", "C1-TPM", "C3-TPM")]
condicion <- datos[, c("1W5-TPM", "1W6-TPM", "1W5-COV")]

# Calcular las medias
datos$Media_Control <- rowMeans(control)
datos$Media_Condicion <- rowMeans(condicion)

# Calcular las desviaciones estándar
datos$DE_Control <- apply(control, 1, sd)
datos$DE_Condicion <- apply(condicion, 1, sd)

# Comprobar los resultados
head(datos)
dim(datos)
### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).
# EJERCICIO 3: Histogramas de expresión génica

# Mostrar los dos histogramas en una sola imagen
par(mfrow = c(1, 2))

# Histograma del control
hist(log10(datos$Media_Control + 1),
     col = "skyblue",
     main = "Expresión génica - Control",
     xlab = "log10(TPM + 1)",
     ylab = "Número de genes")

# Histograma de la condición
hist(log10(datos$Media_Condicion + 1),
     col = "salmon",
     main = "Expresión génica - Condición",
     xlab = "log10(TPM + 1)",
     ylab = "Número de genes")

# Restablecer la disposición de gráficos
par(mfrow = c(1, 1))
### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)
# EJERCICIO 4: Boxplots de expresión génica

# Representar control y condición en una sola imagen
boxplot(log10(datos$Media_Control + 1),
        log10(datos$Media_Condicion + 1),
        names = c("Control", "Condición"),
        col = c("skyblue", "salmon"),
        main = "Comparación de expresión génica",
        ylab = "log10(TPM + 1)")

# Añadir leyenda en la esquina superior izquierda
legend("topleft",
       legend = c("Control", "Condición"),
       fill = c("skyblue", "salmon"),
       bty = "n")
### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.
# EJERCICIO 5: Normalidad y comparación estadística

# Comprobar normalidad de los datos
normalidad_control <- shapiro.test(log10(datos$Media_Control + 1))
normalidad_condicion <- shapiro.test(log10(datos$Media_Condicion + 1))

# Mostrar resultados
normalidad_control
normalidad_condicion
# Comparar control y condición para cada gen
datos$p_value <- sapply(1:nrow(datos), function(i) {
  
  valores_control <- as.numeric(control[i, ])
  valores_condicion <- as.numeric(condicion[i, ])
  
  # Si ambos grupos tienen varianza cero, no se puede aplicar la prueba t
  if (sd(valores_control) == 0 &&
      sd(valores_condicion) == 0) {
    return(NA_real_)
  }
  
  # Prueba t de Welch
  resultado <- tryCatch(
    t.test(valores_control, valores_condicion)$p.value,
    error = function(e) NA_real_
  )
  
  return(resultado)
})

# Comprobar los resultados
head(datos)
dim(datos)
summary(datos$p_value)
### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)
# EJERCICIO 6: Genes diferencialmente expresados

# Ajustar los valores p por comparaciones múltiples
datos$p_ajustado <- p.adjust(datos$p_value, method = "BH")

# Seleccionar genes significativos (FDR < 0.05)
genes_DE <- datos[!is.na(datos$p_ajustado) &
                    datos$p_ajustado < 0.05, ]

# Mostrar el listado de genes
print(genes_DE[, c("Gene Name", "Media_Control",
                   "Media_Condicion", "p_value",
                   "p_ajustado")])

# Contar cuántos genes significativos hay
cat("Número de genes diferencialmente expresados:",
    nrow(genes_DE), "\n")
### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
# EJERCICIO 7: Heatmap de genes diferencialmente expresados

# Seleccionar las seis réplicas
columnas <- c("C2-TPM", "C1-TPM", "C3-TPM",
              "1W5-TPM", "1W6-TPM", "1W5-COV")

# Extraer la expresión de los genes significativos
matriz <- as.matrix(genes_DE[, columnas])

# Identificar los genes
rownames(matriz) <- make.unique(genes_DE[["Gene Name"]])

# Transformación logarítmica
matriz_log <- log10(matriz + 1)

# Normalizar cada gen para comparar patrones de expresión
matriz_z <- t(scale(t(matriz_log)))

# Sustituir valores no definidos de genes sin variación
matriz_z[!is.finite(matriz_z)] <- 0

# Colores del heatmap
colores <- colorRampPalette(c("skyblue", "white", "salmon"))(100)

# Crear el mapa de calor
heatmap(matriz_z,
        Colv = NA,
        scale = "none",
        col = colores,
        labRow = NA,
        labCol = c("Control 1", "Control 2", "Control 3",
                   "Condición 1", "Condición 2", "Condición 3"),
        main = "Genes diferencialmente expresados",
        margins = c(9, 5))
# Interpretación:
# El heatmap muestra diferencias en los patrones de expresión
# génica entre el grupo control y la condición experimental.
# La mayoría de los genes seleccionados presentan una menor
# expresión relativa en la condición, mientras que un grupo
# más reducido presenta una mayor expresión relativa.