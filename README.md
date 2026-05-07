# Muestreo Espacial de Puntos

**Complemento QGIS** para muestreo espacial de puntos con Curva de Hilbert.

[![QGIS](https://img.shields.io/badge/QGIS-3.28%20LTR%20%7C%203.44%20LTR%20%7C%204.0-green)](https://qgis.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/License-GPL%20v2-yellow)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-orange)](metadata.txt)

---

## Descripción

Implements four spatial sampling methods (Systematic Hilbert, Random, Hilbert Groups, K-Means) using the Hilbert Curve. / Implementa cuatro métodos de muestreo espacial (Sistemático Hilbert, Aleatorio, Grupos Hilbert, K-Medias) usando la Curva de Hilbert.

---

## Métodos de Muestreo

| Abrev. | Método | Descripción |
|--------|--------|-------------|
| **SH** | Sistemático Hilbert | Selección a intervalos regulares sobre la curva. Produce patrones dispersos (IVMC > 1,2). |
| **AL** | Aleatorio Simple | Selección aleatoria sin reemplazo del marco muestral. |
| **GH** | Grupos Hilbert | Estratificación 1D sobre el orden de Hilbert. k grupos lineales. |
| **KM** | Grupos K-Medias | Agrupamiento espacial 2D. Requiere `scikit-learn`. |

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
