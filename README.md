# Muestreo Espacial de Puntos (SH, GH, GFC, AL, KM, Estr_Pts-AL, Estr_Pol-AL)

**Complemento QGIS** para muestreo espacial de puntos con siete métodos independientes.

[![QGIS](https://img.shields.io/badge/QGIS-3.28%20LTR%20%7C%203.44%20LTR%20%7C%204.0-green)](https://qgis.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/License-GPL%20v2-yellow)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-orange)](metadata.txt)

---

## Descripción

La **Curva de Hilbert** es el núcleo de los métodos **SH** y **GH**: al ser una línea continua que recorre el área de estudio preservando la localidad espacial, garantiza que puntos cercanos en el orden de la curva también sean geográficamente próximos. Esto mejora la cobertura y representatividad espacial de SH y GH respecto a la selección completamente aleatoria.

El **Aleatorio Simple (AL)** no prioriza cobertura espacial sino representatividad estadística: es el único método donde las fórmulas del MAS (varianza, intervalos de confianza) son directamente válidas sin corrección de diseño. Los métodos **GFC** y **KM** usan otros criterios de ordenamiento —filas NO→SE y proximidad geográfica 2D, respectivamente—. Los métodos **Estr_Pts-AL** y **Estr_Pol-AL** seleccionan aleatoriamente dentro de estratos definidos vía JSON.

En todos los métodos el ordenamiento Hilbert del marco completo se calcula y queda disponible como referencia de diagnóstico en el reporte HTML.

---

## Métodos de Muestreo

| Abrev. | Método | Descripción |
|--------|--------|-------------|
| **SH** | Sistemático Hilbert | Selección a intervalos regulares sobre la curva de Hilbert. Produce patrones dispersos (IVMC > 1,2). |
| **GH** | Grupos Hilbert | Estratificación 1D sobre el orden de Hilbert. k grupos lineales con selección aleatoria dentro de cada grupo. |
| **GFC** | Grupos Fila-Columna | Numeración secuencial 1..N siguiendo filas NO→SE. Divide en k grupos iguales con selección aleatoria. Sin dependencias externas. |
| **AL** | Aleatorio Simple | Selección aleatoria sin reemplazo. La curva de Hilbert se usa solo para el diagnóstico del reporte. |
| **KM** | Grupos K-Medias | Agrupamiento espacial 2D por proximidad geográfica. La Curva de Hilbert no interviene. Requiere `scikit-learn`. |
| **Estr_Pts-AL** | Estratificado por Puntos (JSON) | Tamaños de muestra variables por estrato definidos en JSON, usando un campo categórico de la capa de **puntos** (ej. tipo de cobertura). Selección aleatoria dentro de cada estrato. |
| **Estr_Pol-AL** | Estratificado por Polígono (JSON) | Tamaños de muestra variables por estrato definidos en JSON, usando un campo de la capa de **polígonos** (ej. id_bufer, id_parcela). Selección aleatoria dentro de cada polígono. |

---

## Instalación

### Opción 1 — Desde el Administrador de Complementos de QGIS
*(cuando esté disponible en el repositorio oficial)*

### Opción 2 — Manual desde GitHub

1. Descargar el repositorio como ZIP:
   https://github.com/jfallas56-CR/Muestreo_Espacial_Puntos/archive/refs/heads/main.zip
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

## Dependencias Opcionales

### Shapely — motor acelerado punto-en-polígono

- **QGIS 3.44 LTR y 4.0+**: incluido automáticamente.
- **QGIS 3.28 – 3.40**: instalar manualmente.

```bash
# OSGeo4W Shell (Windows)
pip install shapely

# Terminal Linux/Mac
pip3 install shapely
```

### scikit-learn — método Grupos K-Medias (KM) únicamente

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
  (grados) cancela la ejecución con mensaje de error antes de iniciar.

---

## Parámetros Principales

| Parámetro | Descripción |
|-----------|-------------|
| Capa de puntos | Marco muestral completo |
| Área de muestreo | Polígono que define el área de estudio |
| Método | SH / GH / GFC / AL / KM / Estr_Pts-AL / Estr_Pol-AL |
| Tamaño de muestra (n) | Puntos por iteración (se ignora en métodos JSON) |
| Iteraciones | Número de muestras a generar |
| Orden Hilbert (1–11) | Resolución de la curva. Orden 10 para hasta 50 000+ puntos. No aplica a GFC. |
| Número de grupos (k) | Para GH, GFC y KM. k=0 calcula automáticamente: k = mín(⌈√n⌉, √N), mín. 2 grupos. |
| Campo de Estrato (Puntos) | Solo Estr_Pts-AL. Campo categórico de la capa de puntos. |
| Campo de Estrato (Polígonos) | Solo Estr_Pol-AL. Campo de la capa de polígonos. |
| Tamaños por estrato (JSON) | Solo Estr_Pts-AL y Estr_Pol-AL. Diccionario `{"Estrato 1": n1, "Estrato 2": n2, ...}`. |
| Distancia mínima | Separación mínima entre puntos de la muestra — índice espacial O(n log n). |
| Corrección de borde | Retrae puntos cercanos al límite del área al interior del búfer. |
| Distancia de malla (g) | Para calcular Índice de Cobertura h/g en el reporte. |

---

## Salidas

| Capa | Descripción |
|------|-------------|
| `SH_muestra_01`, `GH_muestra_01`, `GFC_muestra_01`, `AL_muestra_01`, `KM_muestra_01`, `Estr_Pts-AL_muestra_01`, `Estr_Pol-AL_muestra_01` … | Capas de muestra individuales (GPKG temporal). El prefijo corresponde al método y la palabra clave es configurable. |
| Salida: Puntos filtrados | Todos los puntos del marco dentro del área. Campo `Hilbert_idx_proc` (métodos SH/GH/AL/KM/Estr_Pts-AL/Estr_Pol-AL) o `NW_SE_idx_proc` (método GFC). |
| Salida: Ruta Hilbert | Línea de la curva de ordenamiento (**solo métodos SH y GH**). |
| Salida: Orden NO→SE | Líneas por fila de la malla, una por fila (**solo método GFC**). |
| Salida: Rechazados | Puntos eliminados por la restricción de distancia mínima. |
| Reporte HTML | Métricas IVMC, proporcionalidad por grupo, gráficos comparativos, top-3 muestras recomendadas y tiempo de ejecución. |

---

## Reporte HTML

El reporte incluye tres secciones:

1. **Resumen de ejecución**: parámetros de configuración y diagnóstico de
   pre-ejecución (índice h/g, duplicados, colisiones Hilbert).
2. **Panel de selección de muestras**: Top-3 muestras recomendadas con
   criterios de selección, gráficos comparativos entre iteraciones y tabla
   de resultados por iteración.
3. **Apéndices técnicos** (colapsables): marco de validación de la ruta
   (R, σ/h, L/Nh, CV), guía de interpretación del IVMC, metodología de
   cálculo y validación de proporcionalidad por grupo (GH, GFC, KM,
   Estr_Pts-AL, Estr_Pol-AL).

---

## Campos de Salida — Índice de Ordenamiento

| Método | Campo en capa filtrada | Descripción |
|--------|----------------------|-------------|
| SH, GH, AL, KM, Estr_Pts-AL, Estr_Pol-AL | `Hilbert_idx_proc` | Posición del punto en el orden de la curva de Hilbert (0-based). |
| **GFC** | `NW_SE_idx_proc` | ID secuencial 1..N siguiendo el orden fila NO→SE. 1 = punto más al noroeste. |

---

## Notación Numérica

Los valores numéricos en la interfaz y el reporte usan notación española:

- Decimal con **coma**: `1,414`
- Miles con **espacio fino**: `14 747`

---

## Estructura del Repositorio

```
Muestreo_Espacial_Puntos/
├── __init__.py                       # Entry point (classFactory)
├── muestreo_espacial_plugin.py       # Clase principal del complemento
├── muestreo_espacial_provider.py     # Proveedor de Processing
├── Muestreo_Espacial_Puntos.py       # Algoritmo principal
├── metadata.txt                      # Metadatos del complemento (QGIS)
├── icon.png                          # Ícono del complemento (64×64 PNG)
├── LICENSE                           # GPL v2
└── README.md                         # Este archivo
```

---

## Compatibilidad

| QGIS | Qt | Python | GEOS | Estado |
|------|----|--------|------|--------|
| 3.28 LTR | Qt5 | 3.9+ | 3.9+ | ✅ Compatible |
| 3.44 LTR | Qt5 | 3.12 | 3.14 | ✅ Compatible (Shapely incluido) |
| 4.0 | Qt6 | 3.12+ | 3.14 | ✅ Compatible (Shapely incluido) |

---

## Historial de Versiones

### v1.0.0 — Junio 2026 *(versión inicial publicada)*

Lanzamiento inicial en el QGIS Plugin Repository y en GitHub.

**Métodos de muestreo (siete):** SH, GH, GFC, AL, KM, Estr_Pts-AL y
Estr_Pol-AL. Los dos métodos estratificados por JSON (Estr_Pts-AL y
Estr_Pol-AL) permiten asignar tamaños de muestra variables por estrato
mediante un diccionario JSON, usando un campo categórico de la capa de
puntos o de la capa de polígonos como variable de estratificación.

**Validaciones previas (`checkParameterValues`):**
- SRC geográfico cancela la ejecución con mensaje específico
  (recomendación de CRTM05 EPSG:8908 o UTM zona correspondiente).
- Rango del orden Hilbert (1–11).
- Número de iteraciones ≥ 1.
- Dependencia `scikit-learn` para el método K-Medias.
- Formato JSON para los métodos estratificados.
- Coherencia entre número de grupos y tamaño de muestra (NUM_GROUPS ≤
  SAMPLE_SIZE en métodos agrupados).
- Radio de parcela > 0 cuando la corrección de borde está activa.
- Distancia de malla manual no negativa.

**Compatibilidad y motor de cálculo:**
- Aliases Qt5/Qt6 (`QVariant` → `QMetaType`) para QGIS 3.28 LTR,
  3.44 LTR y 4.0.
- Motor Shapely v2 (`contains_xy` ufunc) / v1 (`vectorized.contains`) /
  GEOS `PreparedGeometry` como fallback automático con log del motor
  activo.
- `makeValid()` post-`unaryUnion()` en las tres ramas de procesamiento
  (estándar, estratos individuales y unión global).
- `releaseCache()` para todas las geometrías preparadas, incluyendo las
  de estratos individuales en el método Estratificado por Polígono.

**Robustez y diagnóstico:**
- Gestión cross-thread de visibilidad de capas con `Qt.QueuedConnection`.
- Logging tipificado en operaciones críticas (escritura del reporte
  HTML, extracción de contornos poligonales, limpieza de temporales).
- Mensajes de error con valor recibido, rango esperado y acción
  sugerida (instalación, reproyección, ajuste de parámetro).

**Integración con la caja de Processing:**
- `group()` / `groupId()` para categorización en el grupo "Muestreo
  Espacial".
- `shortHelpString()` con descripción HTML embedida.
- Ayuda por parámetro (`setHelp`) en todos los campos del diálogo.
- `VERSION` como constante de clase, expuesta en el reporte HTML.

**Reporte HTML:**
- Métricas IVMC, índice de cobertura, CV.
- Comparación entre iteraciones y recomendación automática de las tres
  mejores muestras según el objetivo (dispersión o aleatoriedad).
- Umbrales estadísticos NNI/IVMC como constantes de clase
  (`NNI_RANDOM_LO`, `NNI_RANDOM_HI`, `NNI_DISPERSED_OK`).

---

## Autor

**Jorge Fallas**
- Email: jfallas56@gmail.com
- GitHub: [jfallas56-CR](https://github.com/jfallas56-CR)

---

## Licencia

GPL v2 — ver [LICENSE](LICENSE)
