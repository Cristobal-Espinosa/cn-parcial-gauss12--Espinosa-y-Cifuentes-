# Conteo de operaciones

Convencion: una flop por cada operacion aritmetica elemental. El barrido de cada paso va de la columna k en adelante e incluye la actualizacion de b. Los numeros salen de los contadores dentro de gauss.py.

## Sin pivoteo

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
| sustitucion | 78 | 66 | 144 |
| **Total** | 782 | 704 | 1486 |

## Con pivoteo parcial escalado

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
| sustitucion | 78 | 66 | 144 |
| **Total** | 782 | 704 | 1486 |

## Comparacion con la cota

| Corrida | Flops contadas | Cota minima | Diferencia |
|---|---|---|---|
| Sin pivoteo | 1486 | 1486 | 0.00 % |
| Pivoteo escalado | 1486 | 1486 | 0.00 % |

## Sobrecosto del pivoteo (aparte de las flops anteriores)

| Concepto | Cantidad |
|---|---|
| Comparaciones en la busqueda del pivote | 66 |
| Comparaciones para calcular las escalas de fila | 132 |
| Divisiones de escala (valor absoluto de a(i,k) entre s(i)) | 77 |
| Intercambios de fila | 8 |
