# Eliminación gaussiana con pivoteo escalado: sistema denso 12×12

Programa en Python que resuelve un sistema lineal denso A x = b de 12 incógnitas, casi singular (κ₂(A) ≈ 3×10¹³), con eliminación gaussiana propia. Compara la corrida sin pivoteo con la de pivoteo parcial escalado, cuenta las operaciones desde el código y verifica la solución en precisión extendida.

## Autoría y asignatura

Cristobal Espinosa y Steven Cifuentes · Computación Numérica (INF-321), Ingeniería Civil Informática, Universidad Católica del Maule · Prueba Parcial 3 "El pivote", noviembre de 2026

## Requisitos

- Python 3.10 o superior (probado con 3.11).
- Dependencias (numpy 2.2.6 y mpmath 1.3.0), desde la raíz del repositorio:

```
pip install -r requirements.txt
```

## Ejecución

Desde la raíz del repositorio:

```
python src/main.py
```

El programa lee los datos, hace las dos corridas (sin pivoteo y con pivoteo parcial escalado) y escribe en `resultados/`, sin intervención manual:

- `traza_eliminacion.txt`: la matriz después de cada uno de los 11 pasos.
- `conteo_flops.md`: operaciones por paso y sobrecosto del pivoteo.
- `residuo.txt`: residuo en precisión extendida, refinamiento iterativo, dígitos correctos y condicionamiento.
- `mantisa.txt`: bits perdidos de cada corrida.
- `solucion.csv`: las 12 incógnitas con los dígitos que se afirman correctos.

## Estructura de carpetas

```
├─ readme.md          esta página
├─ informe.md         informe técnico y análisis de error
├─ requirements.txt   dependencias con versión fija
├─ datos/             matriz_a.csv (A, 12×12) y vector_b.csv (b, 20 dígitos)
├─ src/               gauss.py, conteo.py, verificacion.py, mantisa.py y main.py (punto de entrada)
└─ resultados/        archivos que genera main.py
```

## Resultado resumido

Solución con pivoteo escalado, tras refinamiento, a 10 dígitos significativos:

| Incógnita | Valor | Incógnita | Valor |
|---|---|---|---|
| x₁ | 0.7384915627 | x₇ | 1.473820592 |
| x₂ | −1.205937417 | x₈ | −0.03849271640 |
| x₃ | 0.004387195200 | x₉ | 0.9152648037 |
| x₄ | −0.09172645380 | x₁₀ | −2.607139485 |
| x₅ | 2.846109372 | x₁₁ | 0.05187362940 |
| x₆ | −0.006281540900 | x₁₂ | −1.084937162 |

| Medida | Valor |
|---|---|
| Flops por corrida (sin y con pivoteo) | 1486 |
| Sobrecosto del pivoteo, aparte | 66 comparaciones de búsqueda, 132 para las escalas, 77 divisiones, 8 intercambios |
| Norma infinito del residuo, pivoteo escalado | 2.29×10⁻¹⁰ |
| Norma infinito del residuo tras 1 pasada de refinamiento | 1.37×10⁻²² |

## Informe

El análisis completo está en [informe.md](informe.md).

## Nota de uso académico

Trabajo realizado con fines académicos para la asignatura Computación Numérica (INF-321) de la Universidad Católica del Maule. No está pensado para otro uso. El uso de herramientas de IA se declara en `informe.md`.
