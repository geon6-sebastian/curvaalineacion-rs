# Curva de Alineación

La curva de alineación representa, en forma ideal, la trayectoria del eje de colimación de un teodolito. A diferencia de las secciones normales y centrales, la curva de alineación no es plana sino que tiene torsión, lo que impide soluciones analíticas directas para los problemas geodésicos directo e inverso. En "curvas.py" se implementan las fórmulas de la curva a partir de su función implícita y sus derivadas, formulando sistemas de ecuaciones diferenciales ordinarias para el cálculo de acimut, longitud y área. Las fórmulas están descritas en el documento ["La curva de alineación"](https://doi.org/10.13140/RG.2.2.20608.39684).
El problema inverso se resuelve eficientemente mediante integración numérica con el método de Dormand-Prince. El problema directo requiere un esquema iterativo de Newton-Raphson bivariado con derivadas numéricas, resultando computacionalmente costoso; además no es una implementación robusta y debe ser siempre comprobado con el problema inverso.
Esta forma de resolver la curva de alineación se extiende a los casos de la sección normal primera y la sección central. Estas curvas suelen ser similares entre sí, excepto cuando los puntos terminales están cerca de ser antipodales. En este ejemplo se representan las curvas de alineación (verde), sección normal (magenta) y sección central (negro) con puntos terminales (latitud, longitud): (-30, 0) y (30, 179).

!["Curvas de alineación, sección normal y central"](./fig02.png)


## Tabla de Contenidos

- [Requisitos](#Rrequisitos)
- [Instalación](#Instalación)
- [Uso](#Uso)
- [Ejemplos](#Ejemplos)
- [Licencia](#Licencia)

## Requisitos

- **Python 3.x**

```bash
pip install numba
pip install pandas
pip install geopandas pyshp simplekml shapely
```

En Linux se debe instalar python-venv. Por ejemplo, en Ubuntu y derivados:

```bash
sudo apt install python3-pip
sudo apt install python3.12-venv
```

donde "3.12" se debe modificar de acuerdo a la versión disponible en su distribución. Luego se debe crear y activar un entorno virtual:

```bash
python3 -m venv miscurvas
source ./miscurvas/bin/activate
```

Una vez activado, se pueden instalar los requisitos con "pip install".

## Instalación

**Clonar el repositorio:**

```bash
git clone https://github.com/geon6-sebastian/curvaalineacion.git
cd curvaalineacion
```

---

## Uso

Para ejecutar el script, utiliza el siguiente comando en la terminal:

```bash
python curvas.py [argumentos]
```

### Argumentos

| Argumento                                 | Descripción                                                                                                          | Requerido        | Default                    |
| ----------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ---------------- | -------------------------- |
| -i, --inverso                             | Ejecutar problema inverso                                                                                            | No               | -                          |
| -d, --directo                             | Ejecutar problema directo                                                                                            | No               | -                          |
| -poly, --poly-sup ('coords.csv')          | Calcula la superficie dentro de un polígono dado en un archivo CSV/TXT                                               | No               | -                          |
| -t, --tipo ['align', 'normal', 'central'] |                                                                                                                      | No               | 'central'                  |
| -P1 (latitud longitud)                    | Punto 1: latitud longitud (en grados decimales). Requerido para -i y -d                                              | Si, para -i y -d | -                          |
| -P2 (latitud longitud)                    | Punto 1: latitud longitud (en grados decimales). Requerido para -i                                                   | Si, para -i      | -                          |
| -e (a, inv_f)                             | Elipsoide: semieje_mayor inversa_aplastamiento (por defecto: GRS80)                                                  |                  | GRS80_a, 298.2572221008827 |
| -o, --output ('nombrearchivo')            | Nombre base para guardar salidas (KMZ, SHP, CSV). Este comando SOBREESCRIBE los archivos existentes del mismo nombre | No               | -                          |
| -az (acimut)                              | Acimut inicial (en grados decimales). Requerido para -d                                                              | Si, para -d      | -                          |
| -s (distancia)                            | Distancia (en metros). Requerido para -d                                                                             | Si, para -d      | -                          |
| -mstep, --max-step (paso)                 | Paso máximo de h para Dormand-Prince en grados decimales                                                             | No               | 0.1                        |

---

## Ejemplos

**Ejemplo con puntos cercanos a ser antipodales (Problema inverso, paso 0.1 grados):**

```bash
python curvas.py -i -P1 -30 0 -P2 30 179 -t align -o curva_align
```

Salida:

```bash
========================================
Acimut (deg): 110.7675169475
Distancia (m): 21656598.9243
Area (m2): 2
Latitud Vértice phi0 (deg): 38.5170077016
Longitud Vértice L0 (deg): 135.6005780565
========================================

Generando archivos: curva_align.* ...
Shapefile guardado como 'curva_align_puntos.shp' con 2389 puntos y 5 columnas de datos.
Archivos generados
```

Estos comandos generan las curvas de la figura más arriba.

```bash
python curvas.py -i -P1 -30 0 -P2 30 179 -t normal -o curva_normal
```

```bash
python curvas.py -i -P1 -30 0 -P2 30 179 -t central -o curva_central
```


**Ejemplo básico de problema directo, paso 0.01 grados:**

```bash
python curvas.py -d -P1 -30 -60 -a 30 -s 5000000 -t align -o align_0.01 -mstep 0.01
```

Para distancias muy largas, el algoritmo es altamente inestable cuando el paso máximo es menor a 0.01.

**Cálculo de la superficie de un polígono uniendo vértices con la curva de alineación:**

```bash
python curvas.py -poly poligono.csv -t align -o poligo
```

Donde el contenido de "poligono.csv" es:
```csv
"Y","X"
-73.3370959549,-70.0386475104
-64.4263843966,-8.3982535173
-46.8962976052,-5.0393662652
-29.483791878,-11.5368935014
-10.7399870261,-13.5643624094
6.1376257255,-20.0001368275
16.5713806345,-28.3861075176
2.0613835089,-45.1506945789
-6.5624498119,-56.1496555295
```
---

La salida del comando es:
```bash
==================================================
Calculando polígono con 9 vértices...
--------------------------------------------------
Arista     Distancia (m)   Acimut (deg)
--------------------------------------------------
1 - 2     2527150.9066    99.31718166
2 - 3     1962417.5168    7.61017844
3 - 4     2013087.0434    -18.51916513
4 - 5     2085855.2856    -6.21369770
5 - 6     1998041.0725    -21.20642523
6 - 7     1472208.4403    -37.67719168
7 - 8     2438801.7893    -129.40724881
8 - 9     1550287.6557    -128.03738171
9 - 1     7473628.1239    -175.69485245
--------------------------------------------------
Superficie del polígono (m2): 34975442749769.69
==================================================

Generando archivos: poligo.* ...
Shapefile guardado como 'poligo_puntos.shp' con 2313 puntos y 5 columnas de datos.
Archivos generados
```

!["Polígono poligo.kmz"](./figpoligo.png)

En Meyer, T. H. (2024, p. 140) [Vector-algebra algorithms... ](https://www.mdpi.com/2673-7418/4/2/8) se proporciona el ejemplo
P1 (-40, 165) P2 (45, 0) reportando una distancia de 18671843.56 m. Pero, en nuestro caso al ejecutar
```bash
python curvas.py -i -P1 -40 165 -P2 45 0 -t align
```
se tiene una inconsistencia de unos 3 m:
```bash
========================================
Acimut (deg): -62.3889781191
Distancia (m): 18671840.3838
Area (m2): -34934108256379
Latitud Vértice phi0 (deg): 48.6480754482
Longitud Vértice L0 (deg): 28.2981412449
========================================
```
En el paper [Vector-algebra algorithms... ](https://www.mdpi.com/2673-7418/4/2/8), Meyer escribe "...the arc length of the curve of alignment was computed by recursively bisecting the curve into linear segments and summing the lengths of the segments until the length converged to submillimeter levels.", pero no tengo forma de reproducir el método de cálculo.
## Licencia

Este proyecto está bajo la Licencia MIT.

---
