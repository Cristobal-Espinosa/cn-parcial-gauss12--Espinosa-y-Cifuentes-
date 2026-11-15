import math
import mpmath as mp

U = 2.0 ** -53
BITS_MANTISA = 53


def norma2(v):
    suma = mp.mpf(0)
    for valor in v:
        suma = suma + valor * valor
    return mp.sqrt(suma)


def error_relativo(x, x_ref):
    dif = []
    for i in range(len(x)):
        dif.append(mp.mpf(x[i]) - x_ref[i])
    return norma2(dif) / norma2(x_ref)


def error_maximo(x, x_ref):
    mayor = mp.mpf(0)
    for i in range(len(x)):
        e = abs(mp.mpf(x[i]) - x_ref[i]) / abs(x_ref[i])
        if e > mayor:
            mayor = e
    return mayor


def bits_perdidos(error):
    razon = float(error / U)
    if razon <= 1:
        return 0.0
    return math.log2(razon)


def perdida_paso1(A, b, fila_pivote):
    n = len(A)
    p = fila_pivote
    mayor_b = 0.0
    mayor_a = 0.0
    for i in range(n):
        if i == p:
            continue
        m = A[i][0] / A[p][0]
        if b[i] != 0:
            r = abs(m * b[p]) / abs(b[i])
            if r > mayor_b:
                mayor_b = r
        for j in range(1, n):
            if A[i][j] != 0:
                r = abs(m * A[p][j]) / abs(A[i][j])
                if r > mayor_a:
                    mayor_a = r
    return math.log2(mayor_b), math.log2(mayor_a)


def escribir_mantisa(ruta, A, b, sin, con, x_ref):
    n = len(A)
    texto = "# Analisis de mantisa y pivoteo\n\n"
    texto = texto + "Mantisa de un float de 64 bits: " + str(BITS_MANTISA) + " bits. "
    texto = texto + "Unidad de redondeo u = 2^-53 = %.3e\n" % U
    texto = texto + "Error de referencia: solucion refinada en precision extendida (verificacion.py).\n\n"

    texto = texto + "## Comparacion de las dos corridas\n\n"
    texto = texto + "| Medida | Sin pivoteo | Pivoteo escalado |\n"
    texto = texto + "|---|---|---|\n"
    m1 = sin["mult_max"]
    m2 = con["mult_max"]
    texto = texto + "| Multiplicador maximo | %.3e | %.3e |\n" % (m1, m2)
    texto = texto + "| log2 del multiplicador maximo (bits) | %.1f | %.1f |\n" % (math.log2(m1), math.log2(m2))
    g1 = sin["crecimiento"]
    g2 = con["crecimiento"]
    texto = texto + "| Factor de crecimiento | %.3e | %.3e |\n" % (g1, g2)
    texto = texto + "| log2 del factor de crecimiento (bits) | %.1f | %.1f |\n" % (math.log2(g1), math.log2(g2))
    e1 = error_relativo(sin["x"], x_ref)
    e2 = error_relativo(con["x"], x_ref)
    texto = texto + "| Error relativo final (norma 2) | " + mp.nstr(e1, 4) + " | " + mp.nstr(e2, 4) + " |\n"
    f1 = error_maximo(sin["x"], x_ref)
    f2 = error_maximo(con["x"], x_ref)
    texto = texto + "| Error relativo maximo por incognita | " + mp.nstr(f1, 4) + " | " + mp.nstr(f2, 4) + " |\n"
    texto = texto + "| Bits perdidos en el resultado: log2(error / u) | %.1f | %.1f |\n" % (bits_perdidos(e1), bits_perdidos(e2))
    texto = texto + "| Bits de mantisa que se conservan | %.1f | %.1f |\n\n" % (BITS_MANTISA - bits_perdidos(e1), BITS_MANTISA - bits_perdidos(e2))

    mayor = 0.0
    fila_mayor = 0
    for i in range(n):
        if abs(A[i][0]) > mayor:
            mayor = abs(A[i][0])
            fila_mayor = i
    razon = mayor / abs(A[0][0])
    texto = texto + "## Por que el pivote a(1,1) es inaceptable\n\n"
    texto = texto + "- |a(1,1)| = %.3e\n" % abs(A[0][0])
    texto = texto + "- Mayor |a(i,1)| de la columna 1 = %.3e (fila %d)\n" % (mayor, fila_mayor + 1)
    texto = texto + "- Cociente = %.3e, que equivale a %.1f bits\n" % (razon, math.log2(razon))
    texto = texto + "- Dividir por a(1,1) multiplica por ese cociente cualquier error de redondeo de la fila pivote.\n\n"

    pb0, pa0 = perdida_paso1(A, b, 0)
    p = con["perm"][0]
    pb1, pa1 = perdida_paso1(A, b, p)
    texto = texto + "## Informacion que se destruye en el paso 1\n\n"
    texto = texto + "Se compara el termino m*a(p,j) que se resta contra el dato a(i,j) que se quiere conservar. "
    texto = texto + "Si el termino es 2^k veces mayor, al restarlo quedan solo " + str(BITS_MANTISA) + " - k bits del dato.\n\n"
    texto = texto + "| Pivote del paso 1 | k en b (bits) | Bits de b que quedan | k en A (bits) | Bits de A que quedan |\n"
    texto = texto + "|---|---|---|---|---|\n"
    texto = texto + "| a(1,1), sin pivoteo | %.1f | %.1f | %.1f | %.1f |\n" % (pb0, BITS_MANTISA - pb0, pa0, BITS_MANTISA - pa0)
    texto = texto + "| Fila %d, pivoteo escalado | %.1f | %.1f | %.1f | %.1f |\n" % (p + 1, pb1, BITS_MANTISA - pb1, pa1, BITS_MANTISA - pa1)

    archivo = open(ruta, "w")
    archivo.write(texto)
    archivo.close()
