# Informe técnico: eliminación gaussiana con mínimo costo operacional y control de mantisa

## 1. Identificación

| Campo | Detalle |
|---|---|
| Integrantes | Cristobal Espinosa, Steven Cifuentes |
| Asignatura | Computación Numérica (INF-321), Ingeniería Civil Informática, UCM |
| Docentes | Juan Manuel Cabezas Contardo |
| Evaluación | Prueba Parcial "El pivote" |
| Sección | 2 |
| Fecha | 20 Noviembre de 2026 |
| Repositorio | https://github.com/Cristobal-Espinosa/cn-parcial-gauss12--Espinosa-y-Cifuentes-  |

## 2. Planteamiento

Se resuelve el sistema lineal denso **A x = b** con n = 12. La matriz no tiene ceros (144 elementos no nulos, 9 cifras significativas cada uno) y b se ingresa con sus 20 dígitos, sin truncar. Los archivos están en `datos/matriz_a.csv` y `datos/vector_b.csv`.

Las comprobaciones de los datos coinciden con los valores de referencia del enunciado: det(A) ≈ 1.3646×10⁹, κ₂(A) ≈ 3.08×10¹³ y un multiplicador máximo sin pivoteo de 2.75×10¹³.

Tres rasgos hacen costoso el sistema:

1. **Escalas de fila incompatibles.** El mayor valor absoluto de cada fila va desde 8.05×10⁻⁵ (fila 5) hasta 9.67×10⁵ (fila 6), un cociente de 1.2×10¹⁰. Al restar un múltiplo de una fila grande a una fila pequeña, el alineamiento de exponentes empuja fuera de los 53 bits de la mantisa los bits bajos del operando pequeño, y esa información no se recupera. Por eso el criterio de pivoteo debe tener en cuenta la escala de cada fila.
2. **Pivote inicial minúsculo.** a(1,1) = 1.3407×10⁻⁸, mientras que el mayor elemento de su columna es 3.690×10⁵ (fila 6). Sin pivoteo, el multiplicador máximo es 2.75×10¹³ y el mayor elemento de la matriz crece por un factor de 2641 durante la eliminación. Ningún pivote se anula, así que el sistema se resuelve sin pivoteo y las dos corridas son comparables.
3. **Casi dependencia lineal.** Ajustando la fila 7 por mínimos cuadrados con las filas 3 y 5 se obtiene fila 7 ≈ 4.18·fila 3 + 583.2·fila 5, con un residuo relativo de 2.3×10⁻³ en norma 2 (el enunciado la describe como una perturbación de orden 10⁻⁴; en ambos casos es una dependencia casi exacta). Esto lleva el valor singular mínimo a σ_min = 7.32×10⁻⁸, frente a σ_max = 2.26×10⁶, y por eso κ₂(A) = 3.08×10¹³ sin que la matriz sea singular.

## 3. Estrategia de eliminación

El algoritmo (`src/gauss.py`, función `gauss(A, b, pivoteo)`) es una eliminación hacia adelante por columnas, k = 1, …, 11, seguida de sustitución regresiva. La misma función hace las dos corridas que pide el enunciado: sin pivoteo (`pivoteo=False`) y con pivoteo parcial escalado (`pivoteo=True`).

**Orden de pasos y por qué no excede el mínimo teórico:**

- **Un multiplicador por fila.** En el paso k se calcula m(i,k) = a(i,k)/a(k,k) una sola vez por fila (66 divisiones en total) y se reutiliza para toda la fila, incluida la actualización de b. Recalcularlo por columna sería trabajo redundante.
- **Barrido desde la columna k.** El paso k no toca las columnas anteriores, ya anuladas. La subdiagonal se deja en 0.0 exacto.
- **Una sola factorización.** Los multiplicadores se guardan en L y el orden de filas en `perm`. Para resolver con otro lado derecho (el residuo del refinamiento) se usa `resolver(info, rhs)`, que reutiliza L y U sin factorizar de nuevo.
- **Sin inversa.** Nunca se calcula A⁻¹.

**Pivoteo parcial escalado** (Burden y Faires): al inicio se calcula una vez la escala s(i) de cada fila, igual a su mayor valor absoluto en la matriz original. Antes del paso k se elige la fila p ≥ k que maximiza |a(p,k)| / s(p) y se intercambia con la fila k, junto con b, s y `perm`. El criterio compara cada candidato contra el tamaño de su propia fila, y por eso no elige mal cuando las escalas difieren en diez órdenes de magnitud. El costo de buscar el pivote y de calcular las escalas se contabiliza aparte (sección 4).

## 4. Conteo de operaciones

**Convención:** una flop por cada operación aritmética elemental (suma, resta, multiplicación o división). El barrido del paso k va de la columna k en adelante e incluye la actualización de b. Los números salen de contadores dentro de `gauss.py`, no de una fórmula, y la tabla la escribe `src/conteo.py` en `resultados/conteo_flops.md`.

| Paso | Mult/Div | Suma/Resta | Total |
|---|---|---|---|
| 1 | 154 | 143 | 297 |
| 2 | 130 | 120 | 250 |
| 3 | 108 | 99 | 207 |
| 4 | 88 | 80 | 168 |
| 5 | 70 | 63 | 133 |
| 6 | 54 | 48 | 102 |
| 7 | 40 | 35 | 75 |
| 8 | 28 | 24 | 52 |
| 9 | 18 | 15 | 33 |
| 10 | 10 | 8 | 18 |
| 11 | 4 | 3 | 7 |
| Sustitución regresiva | 78 | 66 | 144 |
| **Total** | **782** | **704** | **1486** |

La tabla es idéntica en las dos corridas: el pivoteo no cambia la cantidad de flops, solo agrega operaciones de otro tipo.

**Comparación con la cota.** La cota de referencia del enunciado es ≈ 1486 flops para una eliminación densa bien implementada, con actualización de b y sustitución regresiva. El conteo da 1486, una diferencia de 0.00 %. El término dominante 2n³/3 vale 1152 para n = 12; el conteo exacto incluye además términos de orden n², la actualización de b y la sustitución. Si no se calculara la entrada de la columna k que termina en cero, el total sería 1354 (8.9 % menos, también dentro de la tolerancia de ±10 %); se mantuvo la convención del enunciado.

**Sobrecosto del pivoteo, contado aparte:**

| Concepto | Cantidad |
|---|---|
| Comparaciones en la búsqueda del pivote | 66 |
| Comparaciones para calcular las escalas de fila | 132 |
| Divisiones \|a(i,k)\| / s(i) | 77 |
| Intercambios de fila | 8 |

## 5. Análisis de mantisa y pivoteo

Las dos corridas se compararon contra una solución de referencia en precisión extendida (la solución refinada de la sección 8). Un `float` de 64 bits tiene una mantisa de 53 bits y una unidad de redondeo u = 2⁻⁵³ ≈ 1.11×10⁻¹⁶. Los bits perdidos en el resultado se definen como log₂(error relativo / u). Los datos salen de `src/mantisa.py` y `resultados/mantisa.txt`.

| Medida | Sin pivoteo | Pivoteo escalado |
|---|---|---|
| Multiplicador máximo | 2.753×10¹³ | 2.466×10⁹ |
| log₂ del multiplicador máximo (bits) | 44.6 | 31.2 |
| Factor de crecimiento | 2.641×10³ | 1.783 |
| log₂ del factor de crecimiento (bits) | 11.4 | 0.8 |
| Error relativo final (norma 2) | 2.38×10⁻¹⁰ | 2.702×10⁻¹³ |
| Error relativo máximo por incógnita | 7.506×10⁻⁸ | 8.536×10⁻¹¹ |
| Bits perdidos en el resultado (norma 2) | 21.0 | 11.2 |
| Bits perdidos en la peor incógnita (derivado) | 29.3 | 19.6 |

La fila marcada como derivada no la escribe `mantisa.txt`: se obtiene como log₂(error máximo por incógnita / u) con los valores de la fila correspondiente.

**Por qué el pivote a(1,1) es inaceptable.** |a(1,1)| = 1.341×10⁻⁸ y el mayor elemento de su columna vale 3.690×10⁵. El multiplicador que resulta de dividir por a(1,1) es del orden de 2.75×10¹³, equivalente a 44.6 bits. Cada error de redondeo de la fila pivote se amplifica por ese factor al propagarse a las demás filas, y esa amplificación consume la mayor parte de los 53 bits de la mantisa.

**Información que se destruye en el paso 1.** Si el término m·a(p,j) que se resta es 2ᵏ veces mayor que el dato a(i,j) que se quiere conservar, el resultado conserva solo 53 − k bits de ese dato:

| Pivote del paso 1 | k en b (bits) | Bits de b que quedan | k en A (bits) | Bits de A que quedan |
|---|---|---|---|---|
| a(1,1), sin pivoteo | 15.8 | 37.2 | 18.1 | 34.9 |
| Fila 8, pivoteo escalado | 3.2 | 49.8 | 6.1 | 46.9 |

**Por qué difieren las dos corridas.** Con a(1,1) como pivote, el primer paso ya deja 35 bits útiles de A y 37 de b, y esa pérdida se hereda en los 10 pasos siguientes, porque el error de cada operación se propaga a la siguiente. Con el pivote escalado el multiplicador máximo baja 13 bits (de 44.6 a 31.2), el factor de crecimiento casi desaparece (1.78 frente a 2641) y el primer paso conserva 47 bits de A y 50 de b. El resultado final pierde alrededor de 10 bits menos (21.0 frente a 11.2 en norma 2; 29.3 frente a 19.6 en la peor incógnita), lo que equivale a unos tres dígitos decimales: 7 a 10 dígitos correctos sin pivoteo y 10 a 13 con pivoteo. El error en norma 2 mejora unas 880 veces.

Un matiz: los 44.6 bits del multiplicador son una medida de la amplificación posible, no de la pérdida efectiva. El resultado sin pivoteo pierde 21 bits en norma 2 y no 44, porque los errores no se alinean todos en el peor sentido.

## 6. Condicionamiento

Se estimó κ₂(A) con `numpy.linalg.cond` y `numpy.linalg.svd`, usados solo como diagnóstico y no para resolver nada. El cálculo está en `src/verificacion.py` y los valores se escriben en la sección "Condicionamiento y digitos recuperables" de `resultados/residuo.txt`: σ_max = 2.258×10⁶, σ_min = 7.320×10⁻⁸ y **κ₂(A) = 3.085×10¹³**.

**Cota clásica.** Para un método estable, el error relativo satisface aproximadamente ‖δx‖/‖x‖ ≲ κ(A)·ε (Higham, 2002; Trefethen y Bau, 1997). Con ε = 2.22×10⁻¹⁶, la cota vale κ·ε = 6.85×10⁻³, es decir, solo unos **2 dígitos recuperables** (−log₁₀ de la cota = 2.2). El resultado real es muy superior: 10 a 13 dígitos con pivoteo escalado y 7 a 10 sin pivoteo.

**Por qué la cota resulta pesimista aquí.** Hay dos razones.

1. La cota supone el peor caso: una perturbación alineada con la dirección peor condicionada. Un error de redondeo concreto rara vez cumple eso.
2. Buena parte de κ₂(A) proviene de las escalas de fila, no de la dependencia casi lineal entre filas. Si se divide cada fila por su mayor valor absoluto, la matriz equilibrada tiene κ₂ = 8.13×10³, cuatro órdenes de magnitud menor. El pivoteo escalado se comporta aproximadamente como si trabajara sobre esa matriz equilibrada, de modo que el número de condición efectivo es mucho menor. Con él, κ_eq·ε ≈ 1.8×10⁻¹², alrededor de 11.7 dígitos, que concuerda con los 10 a 13 observados. Sin pivoteo, multiplicar por el factor de crecimiento (2641) da ≈ 4.8×10⁻⁹, un orden de magnitud que también es compatible con los 7 a 10 dígitos observados.

Esta lectura es una interpretación respaldada por los números anteriores, no una demostración. Los números efectivos (residuo y dígitos correctos) son los que deciden la calidad de la solución.

## 7. Resultados

Cada incógnita se reporta con los dígitos significativos que se afirman correctos, calculados como ⌊−log₁₀(error relativo)⌋ contra la solución refinada. Los valores están en `resultados/solucion.csv`.

| Incógnita | Valor (pivoteo escalado) | Dígitos correctos, pivoteo escalado | Dígitos correctos, sin pivoteo |
|---|---|---|---|
| x₁ | 0.7384915627000 | 13 | 10 |
| x₂ | −1.20593741680 | 12 | 9 |
| x₃ | 0.004387195200 | 10 | 7 |
| x₄ | −0.091726453800 | 11 | 8 |
| x₅ | 2.84610937250 | 12 | 9 |
| x₆ | −0.006281540900 | 10 | 7 |
| x₇ | 1.47382059160 | 12 | 9 |
| x₈ | −0.038492716400 | 11 | 8 |
| x₉ | 0.915264803700 | 12 | 9 |
| x₁₀ | −2.60713948520 | 12 | 9 |
| x₁₁ | 0.051873629400 | 11 | 8 |
| x₁₂ | −1.08493716250 | 12 | 9 |

Las 12 incógnitas tienen 10 o más dígitos correctos con pivoteo. La solución refinada coincide con números de 10 cifras decimales con una diferencia menor que 10⁻⁴⁸, lo que es coherente con que b se construyó con aritmética exacta a partir de una solución de ese tipo. Los ceros finales de la columna de valores son cifras afirmadas como correctas.

## 8. Verificación

**Residuo en precisión extendida.** r = b − A x se calculó con `mpmath` a 50 dígitos, leyendo A y b desde el texto de los CSV. Norma infinito del residuo:

| Corrida | ‖r‖∞ |
|---|---|
| Sin pivoteo | 1.77175×10⁻⁷ |
| Pivoteo escalado | 2.29088×10⁻¹⁰ |

Relativo a ‖b‖∞ = 7.81×10⁵, el residuo del pivoteo escalado es 2.9×10⁻¹⁶, del orden de la unidad de redondeo: la solución calculada resuelve exactamente un sistema muy cercano al original. Que el error en la solución (10⁻¹¹ a 10⁻¹⁰) sea mucho mayor que ese residuo es lo esperable con un κ tan alto: un residuo pequeño no garantiza un error pequeño.

**Refinamiento iterativo** (sobre la corrida con pivoteo escalado). En cada pasada se calcula el residuo en precisión extendida, se resuelve la corrección con `resolver` (reutilizando L y U) y se actualiza x en precisión extendida:

| Pasada | ‖r‖∞ |
|---|---|
| 0 | 2.29088×10⁻¹⁰ |
| 1 | 1.36624×10⁻²² |
| 2 | 3.92751×10⁻³⁵ |
| 3 | 1.4013×10⁻⁴⁵ |

Una sola pasada baja el residuo doce órdenes de magnitud, y las siguientes lo llevan al límite de los 50 dígitos de trabajo. La solución refinada coincide con `mpmath.lu_solve` (usado solo como comprobación independiente) con una diferencia menor que 10⁻⁴⁵, y es la referencia para contar los dígitos correctos.

**Control con biblioteca.** `numpy.linalg.solve` se usa solo como control. La mayor diferencia relativa entre nuestra solución con pivoteo escalado y la de NumPy fue 7.307×10⁻¹² en la ejecución de los autores. Ese valor depende de la biblioteca de álgebra lineal que use NumPy en cada plataforma, y es la única línea de `resultados/` que cambia al ejecutar el programa en otro entorno.

## 9. Conclusiones

- **Compromiso costo/exactitud.** El pivoteo parcial escalado no cambia las 1486 flops. Agrega 66 comparaciones de búsqueda, 132 para calcular las escalas, 77 divisiones y 8 intercambios, y a cambio reduce el error en norma 2 unas 880 veces (de 2.4×10⁻¹⁰ a 2.7×10⁻¹³), recupera unos 10 bits de mantisa y sube de 7–10 a 10–13 los dígitos correctos. Con ese retorno el pivoteo no se discute en este caso.
- **Refinamiento.** Cuesta una pasada de dos sustituciones triangulares (del orden de n² operaciones, sin refactorizar) más el cálculo del residuo en precisión extendida. Es barato frente a la eliminación, y lleva el residuo de 10⁻¹⁰ a 10⁻²².
- **Lectura del condicionamiento.** κ₂ = 3.1×10¹³ anuncia solo 2 dígitos recuperables y se obtuvieron 10 a 13. La diferencia se explica en buena parte por las escalas de fila (κ de la matriz equilibrada: 8.1×10³), que el pivoteo escalado neutraliza. Sin ese criterio se pierden unos tres dígitos.
- **Qué se haría distinto con n = 1000.** La eliminación costaría 2n³/3 ≈ 6.7×10⁸ flops, mientras que el pivoteo escalado agregaría del orden de 2×10⁶ operaciones (búsqueda, escalas y divisiones), cerca del 0.3 %: su costo relativo baja con n, así que se mantendría. En cambio, los ciclos de Python puro no son viables a ese tamaño. Se usaría una factorización LU por bloques con pivoteo parcial, como la que usa LAPACK en `dgesv` (Anderson et al., 1999), y el refinamiento iterativo se mantendría reutilizando la factorización.
- **Limitaciones.** Los dígitos correctos se miden contra la solución refinada del propio código, comprobada contra una referencia independiente de `mpmath`. El control con NumPy depende de la plataforma.

## 10. Reproducibilidad

- **Entorno de los autores:** Python 3.11, NumPy 2.2.6, mpmath 1.3.0 (Windows). El código requiere Python 3.10 o superior.
- **Instalación y ejecución**, desde la raíz del repositorio:

```
pip install -r requirements.txt
python src/main.py
```

- **Archivos que genera** en `resultados/`: `traza_eliminacion.txt` (matriz tras cada uno de los 11 pasos, en las dos corridas), `conteo_flops.md`, `residuo.txt` (residuo, refinamiento, dígitos correctos y condicionamiento), `mantisa.txt` (archivo de apoyo adicional) y `solucion.csv`.
- **Comprobación.** El programa se ejecutó desde un clon limpio de la rama `main` y reprodujo los mismos resultados. Una ejecución adicional en Linux (Python 3.12) dio archivos idénticos salvo la línea del control con NumPy.
- **Datos que no genera `main.py`.** Los datos descriptivos de la sección 2 (determinante, escalas de fila y ajuste de la fila 7) se calcularon aparte con NumPy y no están en `resultados/`.

### Declaración de uso de IA

Se usó **Claude (Anthropic)** como asistente, tal como permite el enunciado, para lo siguiente:

- Escribir el código inicial de `gauss.py`, `main.py`, `mantisa.py`, `conteo.py` y `verificacion.py`. Para el bloque de condicionamiento de `verificacion.py`, uno de los integrantes usó además otro asistente de IA, y el resultado se revisó y ajustó (por ejemplo, se cambió u por ε = 2⁻⁵² para coincidir con el enunciado).
- Calcular los datos de condicionamiento y de la sección 2.

Lo verificado por los autores: la ejecución del código en su computador y desde un clon limpio de `main`, y la lectura de los resultados generados. Con ayuda de la IA se verificaron, además, los datos de entrada contra los valores del enunciado, la solución contra una referencia independiente de `mpmath`, y el conteo de 1486 flops contra la cota.

### Aclaración puntual
- Los commits realizados por mister_noah pertenecen a Steven Cifuentes.
### Referencias

- Anderson, E. et al. (1999). *LAPACK Users' Guide*, 3.ª ed. SIAM.
- Burden, R. L. y Faires, J. D. *Análisis numérico*. Cengage Learning (pivoteo parcial escalado).
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2.ª ed. SIAM (cota κ·u, factor de crecimiento, escalado de filas).
- Trefethen, L. N. y Bau, D. (1997). *Numerical Linear Algebra*. SIAM (estabilidad de la eliminación gaussiana con pivoteo parcial).
- Enunciado de la Prueba Parcial "El pivote" (INF-321, UCM, 2026).
