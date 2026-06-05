# Muestreo Espacial de Puntos

**Complemento QGIS** para muestreo espacial de puntos con Curva de Hilbert.

[![QGIS](https://img.shields.io/badge/QGIS-3.28%20LTR%20%7C%203.44%20LTR%20%7C%204.0-green)](https://qgis.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/License-GPL%20v2-yellow)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-orange)](metadata.txt)

---

## Description

Spatial point sampling with 7 methods: Systematic Hilbert (SH) and Hilbert Groups (GH) use the Hilbert Curve for selection; Simple Random Sampling (SRS), Row-Column Groups (GFC), K-Means Groups (KM), and JSON-based Stratified by Points (Estr_Pts-AL) and by Polygon (Estr_Pol-AL) select randomly within their strata. 

Key features:
- Edge correction (pullback): moves points away from the polygon
boundary so sampling plots remain entirely within the area.
- Minimum inter-point distance enforced with spatial index O(n log n).
- Corrected NNI (IVMC) computed natively in Python — no external
Processing calls. Correction factor based on real polygon area.
- HTML report with quality metrics, proportionality validation per
group, comparative ranking, and total runtime.
- Result layers loaded as invisible by default.
- Supports grids of 50,000+ points (Hilbert Order = 10).
- Compatible with QGIS 3.28 LTR, 3.44 LTR and 4.0 (Qt5 / Qt6).

---
Sampling methods:
- Systematic Hilbert (SH): fixed-interval selection along the Hilbert
curve. Produces dispersed patterns (NNI > 1.2).
- Hilbert Groups (GH): 1D stratification on the Hilbert order.
k=0: automatic group count (k=min(ceil(sqrt(n)), floor(sqrt(N))), min. 2 groups).
- Simple Random Sampling (SRS): random selection without replacement.
Hilbert ordering used for diagnostic reporting only.
- Row-Column Groups (GFC): sequential 1..N numbering of points
following a NW-to-SE row traversal; partitioned into k equal groups
with random selection within each group. No external dependencies.
Produces a row-by-row line layer (Salida: Orden NO-SE).
- K-Means Groups (KM): 2D spatial clustering. Point selection based
entirely on geographic proximity. Requires scikit-learn.
- Stratified by Points JSON (Estr_Pts-AL): variable sample sizes per
stratum defined by a categorical field of the points layer (e.g.
cover type). Random selection within each stratum.
- Stratified by Polygon JSON (Estr_Pol-AL): variable sample sizes per
stratum defined by a field of the polygon layer (e.g. buffer ID,
parcel ID). Random selection within each polygon.

Optional dependencies:
- shapely: bundled with QGIS 3.44 LTR and 4.0. Earlier versions:
pip install shapely (OSGeo4W Shell)
- scikit-learn: required for K-Means Groups only.
pip install scikit-learn (OSGeo4W Shell)

Requirement: projected CRS in metres (e.g. UTM, CRTM05 EPSG:8908).
── ESPAÑOL ──────────────────────────────────────────────────────────────

Complemento QGIS con siete métodos de muestreo espacial de puntos sobre
una capa de puntos (marco muestral) y un polígono de área de estudio.

SH y GH usan la Curva de Hilbert para el ordenamiento y la selección.
AL selecciona aleatoriamente; el ordenamiento Hilbert se usa solo para el
diagnóstico. GFC ordena los puntos por filas NO→SE (1..N) y los divide en
k grupos iguales con selección aleatoria dentro de cada grupo, sin
dependencias externas. KM agrupa por proximidad geográfica 2D con
K-Medias; la Curva de Hilbert no interviene. Estr_Pts-AL y Estr_Pol-AL
permiten tamaños de muestra variables por estrato definidos vía JSON,
usando un campo de la capa de puntos o de la capa de polígonos como
variable de estratificación, con selección aleatoria dentro de cada uno.

Métodos de muestreo:
- Sistemático Hilbert (SH): selección a intervalos regulares a lo
largo de la curva. Produce patrones dispersos (IVMC > 1,2).
- Grupos Hilbert (GH): estratificación 1D sobre el orden de Hilbert.
k=0: número de grupos automático (k=mín(⌈√n⌉, ⌊√N⌋), mín. 2 grupos).
- Aleatorio Simple (AL): selección aleatoria sin reemplazo.
El ordenamiento Hilbert se usa solo para el diagnóstico del reporte.
- Grupos Fila-Columna (GFC): numeración secuencial 1..N de los puntos
siguiendo un recorrido de filas NO→SE; divididos en k grupos iguales
con selección aleatoria dentro de cada grupo. Sin dependencias
externas. Genera capa de líneas por fila (Salida: Orden NO→SE).
- Grupos K-Medias (KM): agrupamiento espacial 2D con K-Medias.
Selección basada en proximidad geográfica. Requiere scikit-learn.
- Estratificado por Puntos JSON (Estr_Pts-AL): tamaños de muestra
variables por estrato definidos por un campo categórico de la capa
de puntos (ej. tipo de cobertura). Selección aleatoria en cada estrato.
- Estratificado por Polígono JSON (Estr_Pol-AL): tamaños de muestra
variables por estrato definidos por un campo de la capa de polígonos
(ej. ID de búfer, ID de parcela). Selección aleatoria en cada polígono.

Características principales:
- Corrección de borde (retracción): aleja puntos de los límites del
área para que las parcelas queden íntegramente dentro del polígono.
- Distancia mínima entre puntos con índice espacial O(n log n).
- IVMC calculado directamente en Python, sin herramientas externas.
Factor de corrección sobre área real del polígono.
- Reporte HTML con métricas de calidad, validación de proporcionalidad
por grupo, ranking comparativo y tiempo total de ejecución.
- Capas de resultados cargadas invisibles por defecto.
- Soporta mallas de hasta 50 000+ puntos (Orden Hilbert = 10).
- Compatible con QGIS 3.28 LTR, 3.44 LTR y 4.0 (Qt5 / Qt6).

Dependencias opcionales:
- shapely: incluido en QGIS 3.44 LTR y 4.0. Versiones anteriores:
pip install shapely (OSGeo4W Shell)
- scikit-learn: requerido únicamente para Grupos K-Medias.
pip install scikit-learn (OSGeo4W Shell)

Requisito: SRC proyectado en metros (ej: CRTM05 EPSG:8908, UTM).

---

## Instalación

### Opción 1 — Desde el Administrador de Complementos de QGIS
*(cuando esté disponible en el repositorio oficial)*

### Opción 2 — Manual desde GitHub

1. Descargar el repositorio como ZIP:
   ```
   https://github.com/jfallas56-CR/Muestreo_Espacial_Puntos/archive/refs/heads/main.zip
   ```

2. En QGIS: **Complementos → Administrar e instalar complementos → Instalar desde ZIP**

3. Seleccionar el archivo descargado y hacer clic en **Instalar complemento**.

4. Activar el complemento en la lista de complementos instalados.

### Opción 3 — Clonar repositorio

```bash
cd %APPDATA%\QGIS\QGIS3\profiles\default\python\plugins   # Windows
# o
cd ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins  # Linux/Mac

git clone https://github.com/jfallas56-CR/Muestreo_Espacial_Puntos.git
```

---

```bash
# OSGeo4W Shell (Windows)
pip install shapely

# Terminal Linux/Mac
pip3 install shapely
```

### scikit-learn — método Grupos K-Medias

```bash
# OSGeo4W Shell (Windows)
pip install scikit-learn

# PowerShell (Windows) — ajustar versión de QGIS
& "C:\Program Files\QGIS 3.x\bin\python3.exe" -m pip install scikit-learn

# Terminal Linux/Mac
pip3 install scikit-learn
```

Verificar en la Consola Python de QGIS (`Ctrl+Alt+P`):
```python
import sklearn; print(sklearn.__version__)
```

---

## Requisitos

- **SRC proyectado en metros** (ej: CRTM05 EPSG:8908, UTM). Un SRC geográfico
  (grados) cancela la ejecución con mensaje de error.

---

## Parámetros Principales

| Parámetro | Descripción |
|-----------|-------------|
| Capa de puntos | Marco muestral completo |
| Área de muestreo | Polígono que define el área de estudio |
| Método | SH / AL / GH / KM |
| Tamaño de muestra (n) | Puntos por iteración |
| Iteraciones | Número de muestras a generar |
| Orden Hilbert (1-11) | Resolución de la curva. Orden 10 para hasta 50 000+ puntos |
| Número de grupos (k) | Para GH y KM. k=0 calcula automáticamente (mín. 2 grupos) |
| Distancia mínima | Separación mínima entre puntos (O(n log n)) |
| Corrección de borde | Retrae puntos cercanos al límite del área |
| Distancia de malla (g) | Para calcular Índice de Cobertura h/g en el reporte |

---

## Salidas

| Capa | Descripción |
|------|-------------|
| `SH_muestra_01` ... | Capas de muestra individuales (GPKG temporal) |
| Salida: Puntos filtrados | Todos los puntos del marco con índice Hilbert |
| Salida: Ruta Hilbert | Línea de la curva de ordenamiento (SH y GH) |
| Salida: Rechazados | Puntos eliminados por distancia mínima |
| Reporte HTML | Métricas IVMC, gráficos, tiempo de ejecución |

---

## Reporte HTML

El reporte incluye:
- Parámetros de configuración utilizados
- Métricas de la Curva de Hilbert (eficiencia R, colisiones, CV, cobertura)
- IVMC corregido por área real del polígono (no bounding box)
- Gráficos comparativos entre iteraciones
- Identificación de las mejores muestras
- Fecha y tiempo total de ejecución

---

## Notación Numérica

Los valores numéricos usan notación española:
- Decimal con **coma**: `1,414`
- Miles con **espacio fino**: `14 747`

---

## Estructura del Repositorio

```
Muestreo_Espacial_Puntos/
├── __init__.py                   # Punto de entrada del complemento
├── muestreo_espacial_plugin.py   # Clase principal del plugin
├── processing_provider.py        # Proveedor de algoritmos Processing
├── Muestreo_Espacial_Puntos.py   # Algoritmo principal
├── metadata.txt                  # Metadatos del complemento (QGIS)
├── icon.png                      # Ícono del complemento
├── LICENSE                       # GPL v2
└── README.md                     # Este archivo
```

---

## Compatibilidad

| QGIS | Qt | Python | Estado |
|------|-----|--------|--------|
| 3.28 LTR | Qt5 | 3.9+ | ✅ Compatible |
| 3.34 LTR | Qt5 | 3.12 | ✅ Compatible |
| 3.44 LTR | Qt5 | 3.12 | ✅ Compatible (Shapely incluido) |
| 4.0 | Qt6 | 3.12+ | ✅ Compatible (Shapely incluido) |

---

## Autor

**Jorge Fallas**
- Email: jfallas56@gmail.com
- GitHub: [jfallas56-CR](https://github.com/jfallas56-CR)

---

## Licencia

GPL v2 — ver [LICENSE](LICENSE)
