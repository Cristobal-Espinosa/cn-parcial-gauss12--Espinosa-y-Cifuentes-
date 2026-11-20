import math
import mpmath as mp
import numpy as np
from gauss import resolver

mp.mp.dps = 50


def leer_exactos(ruta_a, ruta_b):
    # los valores se leen como texto para no perder ninguna cifra
    A = []
    for linea in open(ruta_a):
        linea = linea.strip()
        if linea != "":
            fila = []
            for valor in linea.split(","):
                fila.append(mp.mpf(valor))
            A.append(fila)
    b = []
    for linea in open(ruta_b):
        linea = linea.strip()
        if linea != "":
            b.append(mp.mpf(linea))
    return A, b


def residuo(A, b, x):
    n = len(b)
    r = []
    for i in range(n):
        suma = mp.mpf(0)
        for j in range(n):
            suma = suma + A[i][j] * x[j]
        r.append(b[i] - suma)
    return r


def norma_inf(v):
    mayor = mp.mpf(0)
    for valor in v:
        if abs(valor) > mayor:
            mayor = abs(valor)
    return mayor


def refinar(info, A, b, x, pasadas):
    # x queda en precision extendida, la correccion se resuelve reutilizando L y U
    x = [mp.mpf(v) for v in x]
    normas = [norma_inf(residuo(A, b, x))]
    for p in range(pasadas):
        r = residuo(A, b, x)
        r_float = [float(v) for v in r]
        d = resolver(info, r_float)
        for i in range(len(x)):
            x[i] = x[i] + mp.mpf(d[i])
        normas.append(norma_inf(residuo(A, b, x)))
    return x, normas


def digitos_correctos(x, x_ref):
    # cantidad de cifras significativas que coinciden con la solucion de referencia
    digitos = []
    for i in range(len(x)):
        error = abs(mp.mpf(x[i]) - x_ref[i]) / abs(x_ref[i])
        if error == 0:
            digitos.append(16)
        else:
            digitos.append(int(math.floor(-mp.log10(error))))
    return digitos


def escribir_verificacion(ruta, A_float, b_float, sin, con, A, b):
    texto = "# Verificacion\n\n"
    texto = texto + "Residuo r = b - A x calculado con mpmath (50 digitos), A y b leidos del texto de los CSV.\n\n"

    texto = texto + "## Residuo de cada corrida (norma infinito)\n\n"
    r_sin = norma_inf(residuo(A, b, sin["x"]))
    r_con = norma_inf(residuo(A, b, con["x"]))
    texto = texto + "- Sin pivoteo: " + mp.nstr(r_sin, 6) + "\n"
    texto = texto + "- Pivoteo escalado: " + mp.nstr(r_con, 6) + "\n\n"

    # refinamiento sobre la corrida con pivoteo (la misma L y U, sin factorizar de nuevo)
    x1, normas1 = refinar(con, A, b, con["x"], 1)
    x5, normas5 = refinar(con, A, b, con["x"], 5)

    texto = texto + "## Refinamiento iterativo (pivoteo escalado)\n\n"
    texto = texto + "| Pasada | Norma infinito del residuo |\n"
    texto = texto + "|---|---|\n"
    i = 0
    for n in normas5:
        texto = texto + "| " + str(i) + " | " + mp.nstr(n, 6) + " |\n"
        i = i + 1
    texto = texto + "\n"
    if normas1[1] < normas1[0]:
        texto = texto + "Con una pasada el residuo baja.\n\n"
    else:
        texto = texto + "Con una pasada el residuo no mejora.\n\n"

    # digitos correctos contra la solucion refinada
    d_sin = digitos_correctos(sin["x"], x5)
    d_con = digitos_correctos(con["x"], x5)
    texto = texto + "## Digitos significativos correctos (contra la solucion refinada)\n\n"
    texto = texto + "| Incognita | Sin pivoteo | Pivoteo escalado |\n"
    texto = texto + "|---|---|---|\n"
    for i in range(len(d_sin)):
        texto = texto + "| x" + str(i + 1) + " | " + str(d_sin[i]) + " | " + str(d_con[i]) + " |\n"
    texto = texto + "\n"

    # condicionamiento: kappa estimado con numpy (solo estimacion)
    A_np = np.array(A_float)
    sigmas = np.linalg.svd(A_np, compute_uv=False)
    kappa = np.linalg.cond(A_np, 2)
    eps = 2.0 ** -52
    cota = kappa * eps
    if cota < 1:
        digitos_cota = int(math.floor(-math.log10(cota)))
    else:
        digitos_cota = 0

    # lo mismo con las filas equilibradas (cada fila dividida por su mayor valor absoluto)
    escalas = np.abs(A_np).max(axis=1)
    A_eq = A_np / escalas[:, None]
    kappa_eq = np.linalg.cond(A_eq, 2)
    cota_eq = kappa_eq * eps
    if cota_eq < 1:
        digitos_cota_eq = int(math.floor(-math.log10(cota_eq)))
    else:
        digitos_cota_eq = 0

    texto = texto + "## Condicionamiento y digitos recuperables\n\n"
    texto = texto + "- Mayor valor singular: %.3e\n" % sigmas[0]
    texto = texto + "- Menor valor singular: %.3e\n" % sigmas[-1]
    texto = texto + "- kappa2(A) estimado con numpy: %.3e\n" % kappa
    texto = texto + "- epsilon de maquina 2^-52: %.3e\n" % eps
    texto = texto + "- Cota del error relativo, kappa * epsilon: %.3e\n" % cota
    texto = texto + "- Digitos recuperables que indica la cota: " + str(digitos_cota) + "\n"
    texto = texto + "- kappa2 con filas equilibradas: %.3e\n" % kappa_eq
    texto = texto + "- Cota con filas equilibradas: %.3e (%d digitos)\n" % (cota_eq, digitos_cota_eq)
    texto = texto + "- Digitos correctos reales sin pivoteo (minimo): " + str(min(d_sin)) + "\n"
    texto = texto + "- Digitos correctos reales con pivoteo escalado (minimo): " + str(min(d_con)) + "\n\n"

    # control con la biblioteca: solo para comparar, no reemplaza la implementacion propia
    x_lib = np.linalg.solve(np.array(A_float), np.array(b_float))
    dif = 0.0
    for i in range(len(x_lib)):
        d = abs(con["x"][i] - x_lib[i]) / abs(x_lib[i])
        if d > dif:
            dif = d
    texto = texto + "## Control con numpy.linalg.solve\n\n"
    texto = texto + "Mayor diferencia relativa entre la solucion propia (pivoteo escalado) y numpy: %.3e\n" % dif

    archivo = open(ruta, "w")
    archivo.write(texto)
    archivo.close()

    return x5, d_con