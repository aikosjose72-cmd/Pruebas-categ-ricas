# Práctica de Pruebas Categóricas en Biotecnología

## Caso 1: Inductores de expresión génica vs Nivel de expresión

Archivo de datos: `caso1_cat.csv` (columnas: `Muestra`, `Tratamiento`, `Expresion`)

Dos variables categóricas (Tratamiento: Control/Inductor_A/Inductor_B; Expresión: Alta/Baja) medidas en muestras independientes → **prueba de Chi-cuadrado de independencia**.

```r
datos <- read.csv("caso1_cat.csv")

# Convertir a factores las variables categóricas
datos$Tratamiento <- as.factor(datos$Tratamiento)
datos$Expresion <- as.factor(datos$Expresion)

# Tabla de contingencia (frecuencias observadas)
tabla <- table(datos$Tratamiento, datos$Expresion)
tabla

# Prueba de Chi-cuadrado de independencia
prueba <- chisq.test(tabla, correct = FALSE)
prueba

# Verificación del supuesto: todas las frecuencias esperadas deben ser >= 5
prueba$expected

# Tamaño del efecto (Cramér's V)
n <- sum(tabla)
k <- min(dim(tabla))
cramersV <- sqrt(prueba$statistic / (n * (k - 1)))
cramersV

# Gráfico de barras agrupadas
barplot(tabla, beside = TRUE, legend.text = TRUE,
        col = c("steelblue", "orange", "seagreen"),
        xlab = "Nivel de expresión", ylab = "Frecuencia",
        main = "Expresión génica según tratamiento")

# Mosaico (proporciones relativas por grupo)
mosaicplot(tabla, main = "Tratamiento vs Expresión", color = TRUE)
```

### Pasos aplicados

1. **Cargar datos** y convertir `Tratamiento` y `Expresion` a factores.
2. **Tabla de contingencia** 3×2 (Control/Inductor_A/Inductor_B × Alta/Baja).
3. **Verificar supuesto** de Chi-cuadrado: todas las frecuencias esperadas ≥ 5 (se cumple: mínimo 14.67).
4. **Ejecutar chisq.test** (sin corrección de Yates, ya que no es tabla 2×2).
5. **Calcular Cramér's V** como medida de tamaño del efecto.
6. **Graficar** con barplot agrupado y mosaico para visualizar la asociación.

### Interpretación

Las frecuencias esperadas mínimas (14.67) cumplen la regla de ≥5 por celda, validando el uso de Chi-cuadrado. El resultado fue **χ² = 14.93, df = 2, p = 0.000573**. Como p < 0.05, se **rechaza la hipótesis nula de independencia**: existe una **asociación estadísticamente significativa entre el tipo de inductor y el nivel de expresión génica**.

El tamaño del efecto (**Cramér's V = 0.353**) indica una asociación de **magnitud moderada**. Al observar la tabla de contingencia, el grupo Control presenta la menor proporción de expresión Alta (16/40 = 40%), mientras que Inductor_B alcanza la mayor proporción (32/40 = 80%), sugiriendo que **ambos inductores incrementan la expresión génica respecto al control, siendo Inductor_B el más efectivo**.

---

## Caso 2: Método de transformación bacteriana vs Éxito de transformación

Archivo de datos: `caso2_cat.csv` (columnas: `Muestra`, `Metodo`, `Resultado`)

Tabla 2×2 con **muestra pequeña (n = 14)**: se espera que alguna frecuencia esperada sea < 5, lo que viola el supuesto de Chi-cuadrado → se usa la **prueba exacta de Fisher**.

```r
datos <- read.csv("caso2_cat.csv")

datos$Metodo <- as.factor(datos$Metodo)
datos$Resultado <- as.factor(datos$Resultado)

# Tabla de contingencia
tabla <- table(datos$Metodo, datos$Resultado)
tabla

# Verificación del supuesto de Chi-cuadrado (frecuencias esperadas)
chisq.test(tabla, correct = FALSE)$expected

# Como hay frecuencias esperadas < 5, se aplica la prueba exacta de Fisher
prueba <- fisher.test(tabla)
prueba

# Razón de momios (odds ratio) estimada por Fisher
prueba$estimate

# Gráfico de barras agrupadas
barplot(tabla, beside = TRUE, legend.text = TRUE,
        col = c("steelblue", "orange"),
        xlab = "Resultado de la transformación", ylab = "Frecuencia",
        main = "Método de transformación vs Resultado")
```

### Pasos aplicados

1. **Cargar datos** y convertir `Metodo` y `Resultado` a factores.
2. **Tabla de contingencia** 2×2 (Electroporación/Choque_térmico × Éxito/Fallo).
3. **Verificar supuesto**: las frecuencias esperadas son todas 3.5, **menores a 5**, por lo que no se cumple el requisito de Chi-cuadrado.
4. **Aplicar fisher.test**, prueba exacta adecuada para tablas 2×2 con muestras pequeñas o frecuencias esperadas bajas.
5. **Graficar** con barplot agrupado.

### Interpretación

Dado que las frecuencias esperadas (3.5 en cada celda) son menores a 5, se utilizó la **prueba exacta de Fisher** en lugar de Chi-cuadrado. El resultado fue **OR = 0.0278, p = 0.0291**. Como p < 0.05, se **rechaza la hipótesis nula**: existe una **asociación estadísticamente significativa entre el método de transformación y el éxito de la misma**.

La tabla de contingencia muestra que **Electroporación tuvo 6 de 7 éxitos**, mientras que **Choque térmico solo tuvo 1 de 7 éxitos**. El odds ratio (0.0278, calculado en dirección Choque_térmico/Electroporación) confirma que las probabilidades de éxito son considerablemente menores con choque térmico. En conclusión, **la electroporación es un método significativamente más eficiente para la transformación bacteriana** en este experimento.

---

## Caso 3: Contaminación antes y después de una intervención (mismos biorreactores)

Archivo de datos: `caso3_cat.csv` (columnas: `Biorreactor`, `Antes`, `Despues`)

Las mismas unidades (biorreactores) se midieron dos veces (**antes** y **después**) → datos **pareados/dependientes** con variable binaria → **prueba de McNemar**.

```r
datos <- read.csv("caso3_cat.csv")

datos$Antes <- as.factor(datos$Antes)
datos$Despues <- as.factor(datos$Despues)

# Tabla de contingencia pareada (mismas unidades medidas dos veces)
tabla <- table(datos$Antes, datos$Despues)
tabla

# Prueba de McNemar (con corrección de continuidad, adecuada para 2x2 pareado)
prueba <- mcnemar.test(tabla, correct = TRUE)
prueba

# Casos discordantes (los únicos informativos para McNemar)
discordantes_mejoraron <- tabla["Si", "No"]  # contaminado -> no contaminado
discordantes_empeoraron <- tabla["No", "Si"] # no contaminado -> contaminado
discordantes_mejoraron
discordantes_empeoraron

# Gráfico de barras comparando antes y después
barplot(rbind(table(datos$Antes), table(datos$Despues)),
        beside = TRUE, legend.text = c("Antes", "Despues"),
        col = c("firebrick", "seagreen"),
        xlab = "Estado de contaminación", ylab = "Frecuencia",
        main = "Contaminación antes vs después de la intervención")
```

### Pasos aplicados

1. **Cargar datos** y convertir `Antes` y `Despues` a factores.
2. **Tabla de contingencia pareada** 2×2, ya que cada biorreactor aporta una medición en ambos momentos (no son muestras independientes).
3. **Aplicar mcnemar.test con corrección de continuidad**, la prueba correcta para datos categóricos dicotómicos medidos en las mismas unidades (antes/después).
4. **Identificar los pares discordantes** (los que cambiaron de estado), que son los únicos que aportan información a la prueba.
5. **Graficar** la frecuencia de contaminación antes vs después.

### Interpretación

De los 40 biorreactores, 12 pasaron de **contaminados a no contaminados** y solo 3 pasaron de **no contaminados a contaminados** (los 25 restantes no cambiaron de estado y no aportan a la prueba). La prueba de McNemar (con corrección de continuidad) arrojó **χ² = 4.267, p = 0.0389**.

Como p < 0.05, se **rechaza la hipótesis nula de que no hubo cambio en la proporción de contaminación**: existe evidencia estadísticamente significativa de que **la intervención redujo la tasa de contaminación** en los biorreactores. La proporción de reactores que mejoraron (pasaron de contaminados a limpios) es significativamente mayor que la de los que empeoraron, lo que respalda la efectividad de la intervención aplicada.

---

## Caso 4: Genotipo bacteriano vs Respuesta a antibiótico

Archivo de datos: `caso4_cat.csv` (columnas: `Muestra`, `Genotipo`, `Respuesta`)

Dos variables categóricas con 3 niveles cada una (Genotipo: WT/MutA/MutB; Respuesta: Sensible/Intermedia/Resistente) en muestras independientes → **prueba de Chi-cuadrado de independencia**.

```r
datos <- read.csv("caso4_cat.csv")

datos$Genotipo <- as.factor(datos$Genotipo)
datos$Respuesta <- as.factor(datos$Respuesta)

# Tabla de contingencia
tabla <- table(datos$Genotipo, datos$Respuesta)
tabla

# Prueba de Chi-cuadrado de independencia
prueba <- chisq.test(tabla, correct = FALSE)
prueba

# Verificación del supuesto: frecuencias esperadas >= 5
prueba$expected

# Tamaño del efecto (Cramér's V)
n <- sum(tabla)
k <- min(dim(tabla))
cramersV <- sqrt(prueba$statistic / (n * (k - 1)))
cramersV

# Gráfico de barras agrupadas
barplot(tabla, beside = TRUE, legend.text = TRUE,
        col = c("steelblue", "orange", "seagreen"),
        xlab = "Respuesta al antibiótico", ylab = "Frecuencia",
        main = "Genotipo vs Respuesta a antibiótico")

# Mosaico
mosaicplot(tabla, main = "Genotipo vs Respuesta", color = TRUE)
```

### Pasos aplicados

1. **Cargar datos** y convertir `Genotipo` y `Respuesta` a factores.
2. **Tabla de contingencia** 3×3.
3. **Verificar supuesto** de Chi-cuadrado: todas las frecuencias esperadas ≥ 5 (mínimo 12).
4. **Ejecutar chisq.test** sin corrección (no es tabla 2×2).
5. **Calcular Cramér's V**.
6. **Graficar** con barplot agrupado y mosaico.

### Interpretación

Todas las frecuencias esperadas superan el umbral de 5 (mínimo 12), validando el uso de Chi-cuadrado. El resultado fue **χ² = 23.89, df = 4, p = 0.000084**. Como p < 0.05, se **rechaza la hipótesis nula de independencia**: existe una **asociación estadísticamente significativa entre el genotipo bacteriano y la respuesta al antibiótico**.

El tamaño del efecto (**Cramér's V = 0.316**) indica una asociación de **magnitud moderada**. Observando la tabla, el genotipo **WT** concentra la mayor proporción de respuesta **Sensible** (24/40 = 60%), mientras que **MutB** presenta la mayor proporción de **Resistente** (23/40 = 57.5%). Esto sugiere que **las mutaciones evaluadas (especialmente MutB) están asociadas con un incremento en la resistencia al antibiótico** respecto a la cepa silvestre.

---

## Caso 5: Pretratamiento enzimático vs Regeneración de explantes

Archivo de datos: `caso5_cat.csv` (columnas: `Explante`, `Pretratamiento`, `Regeneracion`)

Tabla 2×2 (Pretratamiento: Con_enzima/Sin_enzima; Regeneración: Regenera/No_regenera) en muestras independientes, con **n pequeño (24)** → se verifica el supuesto de Chi-cuadrado y, al ser tabla 2×2, se aplica **corrección de continuidad de Yates**.

```r
datos <- read.csv("caso5_cat.csv")

datos$Pretratamiento <- as.factor(datos$Pretratamiento)
datos$Regeneracion <- as.factor(datos$Regeneracion)

# Tabla de contingencia
tabla <- table(datos$Pretratamiento, datos$Regeneracion)
tabla

# Verificación del supuesto: frecuencias esperadas >= 5
chisq.test(tabla, correct = FALSE)$expected

# Al cumplirse el supuesto (todas >= 5) pero ser tabla 2x2,
# se aplica Chi-cuadrado CON corrección de continuidad de Yates
prueba <- chisq.test(tabla, correct = TRUE)
prueba

# Prueba exacta de Fisher como confirmación (recomendable con n pequeño)
fisher.test(tabla)

# Gráfico de barras agrupadas
barplot(tabla, beside = TRUE, legend.text = TRUE,
        col = c("steelblue", "orange"),
        xlab = "Regeneración del explante", ylab = "Frecuencia",
        main = "Pretratamiento enzimático vs Regeneración")
```

### Pasos aplicados

1. **Cargar datos** y convertir `Pretratamiento` y `Regeneracion` a factores.
2. **Tabla de contingencia** 2×2.
3. **Verificar supuesto**: frecuencias esperadas de 5.5 y 6.5, **todas ≥ 5**, por lo que Chi-cuadrado es válido.
4. **Aplicar chisq.test con corrección de Yates**, obligatoria en tablas 2×2 para evitar sobreestimar la significancia.
5. **Confirmar con fisher.test**, dado el tamaño de muestra pequeño (n = 24), como respaldo adicional.
6. **Graficar** con barplot agrupado.

### Interpretación

Las frecuencias esperadas (5.5 y 6.5) cumplen el criterio de ≥5, por lo que se usó Chi-cuadrado con corrección de Yates. El resultado fue **χ² ≈ 0, p = 1.0** (confirmado con Fisher exacto: p = 1.0). Como p > 0.05, **no se rechaza la hipótesis nula**: **no existe evidencia de asociación entre el pretratamiento enzimático y la capacidad de regeneración del explante**.

En la tabla de contingencia, la proporción de regeneración es prácticamente idéntica entre ambos grupos (Con_enzima: 7/12 = 58.3%; Sin_enzima: 6/12 = 50%), una diferencia mínima y atribuible al azar. Se concluye que, **con el tamaño de muestra evaluado, el pretratamiento enzimático no mejora significativamente la tasa de regeneración** de los explantes.

---

## Caso 6 (Bonus): Segregación fenotípica observada vs proporción esperada 9:3:3:1

Archivo de datos: `caso6_cat.csv` (columnas: `Fenotipo`, `Observado`, `Proporcion_esperada`)

Una sola variable categórica con frecuencias observadas que se comparan contra una **proporción teórica conocida** (9:3:3:1, herencia mendeliana con dos genes) → **prueba de Chi-cuadrado de bondad de ajuste**.

```r
datos <- read.csv("caso6_cat.csv")

datos$Fenotipo <- as.factor(datos$Fenotipo)

# Frecuencias observadas y proporciones esperadas
observado <- datos$Observado
proporcion_esperada <- datos$Proporcion_esperada
names(observado) <- datos$Fenotipo

observado
proporcion_esperada
sum(proporcion_esperada)  # debe ser 1

# Frecuencias esperadas bajo H0 (proporción 9:3:3:1)
n <- sum(observado)
esperado <- n * proporcion_esperada
esperado

# Verificación del supuesto: frecuencias esperadas >= 5
esperado >= 5

# Prueba de Chi-cuadrado de bondad de ajuste
prueba <- chisq.test(x = observado, p = proporcion_esperada)
prueba

# Gráfico comparativo observado vs esperado
barplot(rbind(observado, esperado), beside = TRUE,
        names.arg = datos$Fenotipo,
        legend.text = c("Observado", "Esperado"),
        col = c("steelblue", "gray70"),
        xlab = "Fenotipo", ylab = "Frecuencia",
        main = "Segregación fenotípica: observado vs esperado 9:3:3:1")
```

### Pasos aplicados

1. **Cargar datos** ya resumidos en frecuencias por fenotipo.
2. **Definir la hipótesis nula**: los datos siguen la proporción mendeliana 9:3:3:1.
3. **Calcular las frecuencias esperadas** multiplicando n total por cada proporción.
4. **Verificar supuesto**: todas las frecuencias esperadas ≥ 5 (mínimo 10).
5. **Ejecutar chisq.test de bondad de ajuste** (`x` = observados, `p` = proporciones esperadas).
6. **Graficar** observado vs esperado por fenotipo.

### Interpretación

Con n = 160 individuos, las frecuencias esperadas bajo el modelo 9:3:3:1 son 90, 30, 30 y 10 respectivamente, todas ≥ 5, validando el uso de Chi-cuadrado. El resultado fue **χ² = 0.511, df = 3, p = 0.9164**.

Como p > 0.05, **no se rechaza la hipótesis nula**: los datos observados (92, 29, 31, 8) son **consistentes con la proporción mendeliana esperada de 9:3:3:1**. Esto indica que la segregación fenotípica observada en la progenie **no se desvía significativamente** de lo esperado bajo un modelo de herencia de dos genes independientes con dominancia completa, apoyando el patrón de herencia mendeliana clásico para este cruce.

---

## Tabla resumen

| Caso | Variables | Prueba usada | Resultado | Conclusión |
|---|---|---|---|---|
| 1 | Tratamiento vs Expresión | Chi-cuadrado independencia (3×2) | χ²=14.93, p=0.00057 | Asociación significativa |
| 2 | Método vs Resultado | Fisher exacto (esperados <5) | p=0.0291 | Electroporación más efectiva |
| 3 | Antes vs Después (pareado) | McNemar | χ²=4.27, p=0.0389 | La intervención redujo la contaminación |
| 4 | Genotipo vs Respuesta | Chi-cuadrado independencia (3×3) | χ²=23.89, p<0.001 | Asociación significativa (WT sensible, MutB resistente) |
| 5 | Pretratamiento vs Regeneración | Chi-cuadrado con Yates (2×2) | p=1.0 | Sin asociación |
| 6 (bonus) | Fenotipo observado vs esperado | Chi-cuadrado bondad de ajuste | χ²=0.51, p=0.916 | Consistente con proporción 9:3:3:1 |
