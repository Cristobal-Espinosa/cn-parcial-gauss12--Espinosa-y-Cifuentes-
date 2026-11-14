import copy


def gauss(A, b, pivoteo):
    A = copy.deepcopy(A)
    b = b[:]
    n = len(A)

    cont = {
        "mult_div": 0,
        "sum_res": 0,
        "comparaciones": 0,
        "comp_escala": 0,
        "div_escala": 0,
        "intercambios": 0,
    }
    pasos = []
    traza = []

    max_inicial = 0.0
    for i in range(n):
        for j in range(n):
            if abs(A[i][j]) > max_inicial:
                max_inicial = abs(A[i][j])
    max_actual = max_inicial
    mult_max = 0.0

    if pivoteo:
        s = []
        for i in range(n):
            mayor = abs(A[i][0])
            for j in range(1, n):
                cont["comp_escala"] += 1
                if abs(A[i][j]) > mayor:
                    mayor = abs(A[i][j])
            s.append(mayor)

    for k in range(n - 1):
        md_antes = cont["mult_div"]
        sr_antes = cont["sum_res"]
        comp_antes = cont["comparaciones"]

        if pivoteo:
            fila_p = k
            cont["div_escala"] += 1
            mejor = abs(A[k][k]) / s[k]
            for i in range(k + 1, n):
                cont["div_escala"] += 1
                cont["comparaciones"] += 1
                r = abs(A[i][k]) / s[i]
                if r > mejor:
                    mejor = r
                    fila_p = i
            if fila_p != k:
                A[k], A[fila_p] = A[fila_p], A[k]
                b[k], b[fila_p] = b[fila_p], b[k]
                s[k], s[fila_p] = s[fila_p], s[k]
                cont["intercambios"] += 1

        for i in range(k + 1, n):
            m = A[i][k] / A[k][k]
            cont["mult_div"] += 1
            if abs(m) > mult_max:
                mult_max = abs(m)

            for j in range(k, n):
                A[i][j] = A[i][j] - m * A[k][j]
                cont["mult_div"] += 1
                cont["sum_res"] += 1
            b[i] = b[i] - m * b[k]
            cont["mult_div"] += 1
            cont["sum_res"] += 1

            A[i][k] = 0.0

            for j in range(n):
                if abs(A[i][j]) > max_actual:
                    max_actual = abs(A[i][j])

        pasos.append({
            "paso": k + 1,
            "mult_div": cont["mult_div"] - md_antes,
            "sum_res": cont["sum_res"] - sr_antes,
            "comparaciones": cont["comparaciones"] - comp_antes,
        })
        traza.append(copy.deepcopy(A))

    md_antes = cont["mult_div"]
    sr_antes = cont["sum_res"]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        suma = b[i]
        for j in range(i + 1, n):
            suma = suma - A[i][j] * x[j]
            cont["mult_div"] += 1
            cont["sum_res"] += 1
        x[i] = suma / A[i][i]
        cont["mult_div"] += 1
    pasos.append({
        "paso": "sustitucion",
        "mult_div": cont["mult_div"] - md_antes,
        "sum_res": cont["sum_res"] - sr_antes,
        "comparaciones": 0,
    })

    info = {
        "x": x,
        "U": A,
        "contadores": cont,
        "pasos": pasos,
        "traza": traza,
        "mult_max": mult_max,
        "crecimiento": max_actual / max_inicial,
    }
    return info
