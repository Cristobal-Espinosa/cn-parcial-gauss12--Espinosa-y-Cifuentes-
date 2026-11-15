import os
import mpmath as mp
from gauss import gauss
from conteo import escribir_conteo
from verificacion import leer_exactos, escribir_verificacion
from mantisa import escribir_mantisa

carpeta = os.path.dirname(os.path.abspath(__file__))
raiz = os.path.dirname(carpeta)
ruta_a = os.path.join(raiz, "datos", "matriz_a.csv")
ruta_b = os.path.join(raiz, "datos", "vector_b.csv")
ruta_resultados = os.path.join(raiz, "resultados")

A_exacta, b_exacto = leer_exactos(ruta_a, ruta_b)
A = []
for fila in A_exacta:
    nueva = []
    for valor in fila:
        nueva.append(float(valor))
    A.append(nueva)
b = []
for valor in b_exacto:
    b.append(float(valor))

print("Matriz leida:", len(A), "x", len(A[0]))
print("Vector b leido:", len(b), "valores")

sin = gauss(A, b, False)
con = gauss(A, b, True)

print()
print("Corrida sin pivoteo")
print("  multiplicador maximo:", sin["mult_max"])
print("  factor de crecimiento:", sin["crecimiento"])
print("  flops:", sin["contadores"]["mult_div"] + sin["contadores"]["sum_res"])

print()
print("Corrida con pivoteo parcial escalado")
print("  multiplicador maximo:", con["mult_max"])
print("  factor de crecimiento:", con["crecimiento"])
print("  flops:", con["contadores"]["mult_div"] + con["contadores"]["sum_res"])
print("  intercambios de fila:", con["contadores"]["intercambios"])

archivo = open(os.path.join(ruta_resultados, "traza_eliminacion.txt"), "w")
for nombre, corrida in (("SIN PIVOTEO", sin), ("PIVOTEO PARCIAL ESCALADO", con)):
    archivo.write("=== " + nombre + " ===\n")
    paso = 1
    for matriz in corrida["traza"]:
        archivo.write("\nMatriz despues del paso " + str(paso) + "\n")
        for fila in matriz:
            texto = ""
            for valor in fila:
                texto = texto + "%14.6e" % valor
            archivo.write(texto + "\n")
        paso = paso + 1
    archivo.write("\n")
archivo.close()

escribir_conteo(os.path.join(ruta_resultados, "conteo_flops.md"), sin, con)

x_refinada, digitos = escribir_verificacion(os.path.join(ruta_resultados, "residuo.txt"), A, b, sin, con, A_exacta, b_exacto)

escribir_mantisa(os.path.join(ruta_resultados, "mantisa.txt"), A, b, sin, con, x_refinada)

archivo = open(os.path.join(ruta_resultados, "solucion.csv"), "w")
archivo.write("incognita,valor,digitos_correctos\n")
for i in range(len(x_refinada)):
    archivo.write("x" + str(i + 1) + "," + mp.nstr(x_refinada[i], digitos[i]) + "," + str(digitos[i]) + "\n")
archivo.close()

print()
print("Listo: se escribieron traza_eliminacion.txt, conteo_flops.md, residuo.txt, mantisa.txt y solucion.csv en resultados/")
