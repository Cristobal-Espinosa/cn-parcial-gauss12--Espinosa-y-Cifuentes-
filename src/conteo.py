COTA_MINIMA = 1486


def tabla_pasos(corrida):
    texto = "| Paso | Mult/Div | Suma/Resta | Total |\n"
    texto = texto + "|---|---|---|---|\n"
    for p in corrida["pasos"]:
        total = p["mult_div"] + p["sum_res"]
        texto = texto + "| " + str(p["paso"]) + " | " + str(p["mult_div"]) + " | " + str(p["sum_res"]) + " | " + str(total) + " |\n"
    c = corrida["contadores"]
    total = c["mult_div"] + c["sum_res"]
    texto = texto + "| **Total** | " + str(c["mult_div"]) + " | " + str(c["sum_res"]) + " | " + str(total) + " |\n"
    return texto


def escribir_conteo(ruta, sin, con):
    c_sin = sin["contadores"]
    c_con = con["contadores"]
    total_sin = c_sin["mult_div"] + c_sin["sum_res"]
    total_con = c_con["mult_div"] + c_con["sum_res"]

    texto = "# Conteo de operaciones\n\n"
    texto = texto + "Convencion: una flop por cada operacion aritmetica elemental. "
    texto = texto + "El barrido de cada paso va de la columna k en adelante e incluye la actualizacion de b. "
    texto = texto + "Los numeros salen de los contadores dentro de gauss.py.\n\n"

    texto = texto + "## Sin pivoteo\n\n"
    texto = texto + tabla_pasos(sin) + "\n"

    texto = texto + "## Con pivoteo parcial escalado\n\n"
    texto = texto + tabla_pasos(con) + "\n"

    texto = texto + "## Comparacion con la cota\n\n"
    texto = texto + "| Corrida | Flops contadas | Cota minima | Diferencia |\n"
    texto = texto + "|---|---|---|---|\n"
    dif_sin = (total_sin - COTA_MINIMA) / COTA_MINIMA * 100
    dif_con = (total_con - COTA_MINIMA) / COTA_MINIMA * 100
    texto = texto + "| Sin pivoteo | " + str(total_sin) + " | " + str(COTA_MINIMA) + " | %.2f %% |\n" % dif_sin
    texto = texto + "| Pivoteo escalado | " + str(total_con) + " | " + str(COTA_MINIMA) + " | %.2f %% |\n\n" % dif_con

    texto = texto + "## Sobrecosto del pivoteo (aparte de las flops anteriores)\n\n"
    texto = texto + "| Concepto | Cantidad |\n"
    texto = texto + "|---|---|\n"
    texto = texto + "| Comparaciones en la busqueda del pivote | " + str(c_con["comparaciones"]) + " |\n"
    texto = texto + "| Comparaciones para calcular las escalas de fila | " + str(c_con["comp_escala"]) + " |\n"
    texto = texto + "| Divisiones de escala (valor absoluto de a(i,k) entre s(i)) | " + str(c_con["div_escala"]) + " |\n"
    texto = texto + "| Intercambios de fila | " + str(c_con["intercambios"]) + " |\n"

    archivo = open(ruta, "w")
    archivo.write(texto)
    archivo.close()
