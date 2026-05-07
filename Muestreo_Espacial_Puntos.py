# -*- coding: utf-8 -*-
"""
Complemento QGIS: Muestreo Espacial de Puntos por Curva de Hilbert
==================================================================
Autor:    Jorge Fallas (jfallas56@gmail.com)
Versión:  v1.0.0
Fecha:    2026-05-07
QGIS:     3.28 LTR+ / 3.44 LTR / 4.0 (Qt5 y Qt6) compatible
Python:   3.10+
Deps:     scikit-learn (opcional, método K-Medias)
          shapely (opcional, motor acelerado — incluido en QGIS 3.44+ y 4.0)

Propósito:
    Muestreo espacial de puntos sobre áreas poligonales con cuatro
    métodos: Sistemático Hilbert, Aleatorio Simple, Grupos Hilbert y
    Grupos K-Medias. Genera reporte HTML con métricas de calidad espacial
    (IVMC, índice de cobertura, CV, índice de eficiencia R).

Parámetros de entrada principales:
    INPUT          : Capa de puntos (marco muestral)
    SAMPLING_AREA  : Capa poligonal de área de muestreo
    METHOD         : Método de muestreo (0-3)
    SAMPLE_SIZE    : Tamaño de muestra por iteración
    NUM_ITERATIONS : Número de iteraciones a generar

Salidas:
    OUTPUT_ALL          : Todos los puntos filtrados con índice Hilbert
    OUTPUT_HILBERT_PATH : Ruta Hilbert (Sistemático y Grupos Hilbert)
    OUTPUT_REJECTED     : Puntos rechazados por distancia mínima
    OUTPUT_HTML_REPORT  : Reporte HTML con métricas de calidad

Historial de cambios:
    v1.0.0 (2026-05-07):
        - Notación numérica española (coma decimal, espacio miles) en ayuda
        - Orden Hilbert: documentado soporte hasta 50 000+ puntos
        - Número de grupos k=0: agregado "mínimo 2 grupos" a ayuda y código
        - Tabla atributos: dist_vecino_m (2 dec), IVMC (4 dec), dist_borde (2 dec)
        - Ruta Hilbert: tabla atributos con 5 campos descriptivos (ya no vacía)
        - Puntos filtrados: eliminados campos NULL no utilizados (solo Hilbert_idx
          + campos originales + dist_borde opcional)
        - Motor activo: línea eliminada de cabecera del panel de ayuda
        - VERSION eliminada del log de inicio de processAlgorithm
        - Verificación completa de español en interfaz, ayuda y HTML
    v0.9.0 (2026-05-06):
        Versión previa de desarrollo, corrección de borde,
        motor Shapely/PreparedGeometry, visibilidad de capas, reporte HTML.
"""

import math
import random
import webbrowser
import gc
import copy
import os
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from enum import IntEnum
from datetime import datetime
from pathlib import Path

# --- INICIO: Verificación de Dependencia (K-Means) ---
try:
    from sklearn.cluster import KMeans
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
# --- FIN: Verificación de Dependencia ---

from qgis import processing
from qgis.PyQt.QtCore import QCoreApplication, QVariant, Qt, QObject, pyqtSlot
from qgis.core import (
    QgsProcessing, QgsFeatureSink, QgsProcessingException,
    QgsProcessingAlgorithm, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFeatureSink, QgsProcessingParameterNumber,
    QgsProcessingParameterEnum, QgsProcessingParameterDistance,
    QgsProcessingParameterBoolean, QgsProcessingParameterString,
    QgsFeature, QgsGeometry, QgsPointXY, QgsFields, QgsField,
    QgsWkbTypes, QgsProcessingFeedback, QgsProcessingContext,
    QgsDistanceArea, QgsVectorFileWriter, QgsCoordinateReferenceSystem,
    QgsProcessingParameterFileDestination, QgsMapLayer, QgsVectorLayer,
    QgsCoordinateTransform, QgsProject, QgsFeatureRequest,
    QgsProcessingUtils, QgsRectangle, QgsSpatialIndex,
    QgsProcessingLayerPostProcessorInterface, QgsLayerTreeLayer
)

# --- H-11: Aliases de compatibilidad Qt5 (QVariant) → Qt6 (QMetaType) ---
try:
    from qgis.PyQt.QtCore import QMetaType
    _INT_TYPE = QMetaType.Type.Int
    _STR_TYPE = QMetaType.Type.QString
    _DBL_TYPE = QMetaType.Type.Double
except (ImportError, AttributeError):
    _INT_TYPE = QVariant.Int
    _STR_TYPE = QVariant.String
    _DBL_TYPE = QVariant.Double

try:
    _TypeVectorPolygon = QgsProcessing.TypeVectorPolygon
    _TypeVectorPoint   = QgsProcessing.TypeVectorPoint
except AttributeError:
    _TypeVectorPolygon = QgsProcessing.TypeVector
    _TypeVectorPoint   = QgsProcessing.TypeVector

# --- H-04: Detección de motor Shapely para punto-en-polígono acelerado ---
try:
    import shapely
    import numpy as _np
    _SHAPELY_VERSION = shapely.__version__
    _shapely_major = int(_SHAPELY_VERSION.split('.')[0])
    # Shapely >= 2.0: usar contains_xy (ufunc directo, vectorized deprecado en 2.1)
    # Shapely 1.x: usar shapely.vectorized.contains (API legacy)
    if _shapely_major >= 2:
        import shapely as _shmod
        _shapely_contains_xy = _shmod.contains_xy  # función ufunc nativa
        _shapely_use_v2 = True
    else:
        import shapely.vectorized as _shvec
        _shapely_use_v2 = False
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False
    _SHAPELY_VERSION = None
    _shapely_use_v2 = False

# --- Clases de Configuración y Tipos ---

class SamplingMethod(IntEnum):
    SYSTEMATIC_HILBERT = 0
    SIMPLE_RANDOM = 1
    STRATIFIED_HILBERT = 2
    KMEANS_GROUPS = 3

# H-08: SortingMethod.PEANO eliminado — no implementado. Solo existe HILBERT.
# Si se implementa Peano en el futuro, restaurar el enum y agregar la lógica
# en _sort_points() con elif params.sorting_method == SortingMethod.PEANO.

class DistanceUnit(IntEnum):
    METERS = 0
    KILOMETERS = 1
    FEET = 2
    MILES = 3

@dataclass
class SamplingParameters:
    project_name: str
    perimeter_distance: float
    plot_radius: float
    min_sample_distance: float
    hilbert_order: int
    num_groups: int
    sample_size: int
    num_iterations: int
    method: SamplingMethod
    apply_pullback: bool
    select_inside_buffer: bool
    calc_dist_all_points: bool
    manual_grid_spacing: float

@dataclass
class StratumInfo:
    start: int
    end: int
    size: int

@dataclass
class PointData:
    id: int
    x: float
    y: float
    hilbert_index: int = 0
    original_index: int = 0 
    was_corrected: bool = False
    correction_distance: float = 0.0
    rejected_in_iteration: int = 0
    _geom: Optional[QgsGeometry] = None
    
    def get_geometry(self) -> QgsGeometry:
        if self._geom is None:
            self._geom = QgsGeometry.fromPointXY(QgsPointXY(self.x, self.y))
        return self._geom

    def set_geometry(self, geom: QgsGeometry):
        self._geom = geom
        p = geom.asPoint()
        self.x = p.x()
        self.y = p.y()

    def __deepcopy__(self, memo):
        """B-04: deepcopy seguro - no copia QgsGeometry (objeto C++ de QGIS).
        La geometria se recrea de forma lazy en get_geometry() cuando se necesite."""
        nuevo = PointData(
            id=self.id, x=self.x, y=self.y,
            hilbert_index=self.hilbert_index,
            original_index=self.original_index,
            was_corrected=self.was_corrected,
            correction_distance=self.correction_distance,
            rejected_in_iteration=self.rejected_in_iteration
        )  # _geom queda None; se regenera lazy via get_geometry()
        return nuevo

@dataclass
class SamplePoint:
    feature: QgsFeature 
    original_geometry: QgsGeometry 
    was_corrected: bool = False
    correction_distance: float = 0.0
    rejected_in_iteration: int = 0
    
    def __post_init__(self):
        if self.original_geometry is None:
            self.original_geometry = QgsGeometry(self.feature.geometry())

@dataclass
class KMeansValidationData:
    iteration: int
    cluster_stats: List[Dict]

# --- Clases de Utilidades ---

class HilbertCurveCalculator:
    @staticmethod
    def xy_to_hilbert(x: int, y: int, order: int) -> int:
        n = 1 << order
        d = 0
        s = n >> 1
        while s > 0:
            rx = (x & s) > 0
            ry = (y & s) > 0
            d += s * s * ((3 * rx) ^ ry)
            x, y = HilbertCurveCalculator._rotate_hilbert(s, x, y, rx, ry)
            s >>= 1
        return d

    @staticmethod
    def _rotate_hilbert(s: int, x: int, y: int, rx: bool, ry: bool) -> Tuple[int, int]:
        if not ry:
            if rx:
                x = (s - 1) - x
                y = (s - 1) - y
            x, y = y, x
        return x, y

class GeometryUtils:
    @staticmethod
    def get_polygon_boundary_as_line(
        poly: QgsGeometry,
        feedback: QgsProcessingFeedback
    ) -> Optional[QgsGeometry]:
        """Obtiene el contorno del polígono como línea.
        Intenta 3 métodos en orden de preferencia para máxima compatibilidad
        entre versiones de QGIS (boundary() no existe en todas las versiones).
        """
        if not poly or poly.isNull():
            return None

        # Método 1: boundary() — disponible en QGIS >= 3.x con GEOS >= 3.3
        # Usar hasattr() evita la excepción ruidosa en versiones sin este método.
        if hasattr(poly, 'boundary'):
            try:
                boundary = poly.boundary()
                if boundary and not boundary.isNull() and not boundary.isEmpty():
                    return boundary
            except Exception as e:
                if feedback:
                    feedback.pushWarning(
                        f"[GeometryUtils] boundary() falló (GEOS): {e} — usando método alternativo."
                    )

        # Método 2: convertToType — disponible en todas las versiones de QGIS 3.x
        try:
            result = poly.convertToType(QgsWkbTypes.LineGeometry, True)
            if result and not result.isNull() and not result.isEmpty():
                return result
        except Exception as e:
            if feedback:
                feedback.pushWarning(
                    f"[GeometryUtils] convertToType() falló: {e} — usando método de anillos."
                )

        # Método 3: extraer anillo exterior como QgsPolylineXY → QgsGeometry
        # Fallback final; funciona con cualquier QgsGeometry de tipo polígono.
        try:
            poly_geom = poly.asPolygon()
            if poly_geom and poly_geom[0]:
                return QgsGeometry.fromPolylineXY(poly_geom[0])
            # Multipolígono: usar primer parte
            multi = poly.asMultiPolygon()
            if multi and multi[0] and multi[0][0]:
                return QgsGeometry.fromPolylineXY(multi[0][0])
        except Exception as e:
            if feedback:
                feedback.pushWarning(
                    f"[GeometryUtils] Todos los métodos de contorno fallaron: {e}. "
                    "El campo dist_borde_m no está disponible."
                )
        return None

    @staticmethod
    def is_geometry_valid(geom: QgsGeometry) -> bool:
        if not geom or geom.isNull():
            return False
        try:
            return not geom.isEmpty()
        except Exception:
            return False

class DistanceCalculator:
    def __init__(self, distance_area: QgsDistanceArea):
        self.da = distance_area
    def calculate_distance(self, point1: QgsPointXY, point2: QgsPointXY) -> float:
        if self.da.sourceCrs().isGeographic(): return self.da.measureLine(point1, point2)
        return math.hypot(point2.x() - point1.x(), point2.y() - point1.y())

# --- Algoritmo Principal de QGIS ---

# ---------------------------------------------------------------------------
# Handler de visibilidad de capas — QObject con slot decorado.
# Necesario para Qt.QueuedConnection cross-thread confiable en PyQt5/6.
# ---------------------------------------------------------------------------
class _LayerVisibilityHandler(QObject):
    """Recibe la señal layerWasAdded en el hilo GUI y oculta las capas
    del framework mientras carga las capas de muestra en orden.

    Es un QObject con @pyqtSlot — la única forma confiable de recibir
    señales cross-thread (hilo Processing → hilo GUI) en PyQt5/6.
    """

    def __init__(self, framework_names: set, pending: list, expected: int = 1,
                 report_gen=None, t0: float = 0.0):
        super().__init__()
        self._names      = framework_names
        self._pending    = pending
        self._expected   = max(1, expected)
        self._count      = 0
        self._done       = False
        self._report_gen = report_gen  # ReportGenerator para rewrite al finalizar
        self._t0         = t0          # tiempo de inicio para calcular elapsed real
        from qgis.PyQt.QtCore import QCoreApplication
        self.moveToThread(QCoreApplication.instance().thread())

    @pyqtSlot('QgsMapLayer*')
    def on_layer_added(self, layer):
        if self._done or not layer:
            return
        if layer.name() in self._names:
            self._count += 1
            # Ocultar inmediatamente
            root = QgsProject.instance().layerTreeRoot()
            node = root.findLayer(layer.id())
            if node:
                node.setItemVisibilityChecked(False)

        if self._count >= self._expected and not self._done:
            self._done = True
            try:
                QgsProject.instance().layerWasAdded.disconnect(self.on_layer_added)
            except Exception:
                pass
            self._load_samples()

    def _load_samples(self):
        """Carga capas de muestra en orden 1..N, invisibles.
        Al finalizar, reescribe el HTML con el tiempo real de ejecución
        y abre el navegador — este es el momento más tardío del proceso.
        """
        import time as _time
        root    = QgsProject.instance().layerTreeRoot()
        project = QgsProject.instance()
        for path, name in reversed(self._pending):
            try:
                lyr = QgsVectorLayer(path, name, 'ogr')
                if not lyr.isValid():
                    continue
                project.addMapLayer(lyr, False)
                node = QgsLayerTreeLayer(lyr)
                node.setItemVisibilityChecked(False)
                root.insertChildNode(0, node)
            except Exception:
                pass

        # Calcular tiempo real: _load_samples corre en el hilo GUI,
        # después de que QGIS cargó todas las capas — es el momento correcto.
        if self._report_gen is not None and self._t0 > 0:
            try:
                elapsed = _time.time() - self._t0
                self._report_gen.rewrite_with_time(elapsed)
            except Exception:
                pass
            self._report_gen = None


class MuestreoEspacialPuntos(QgsProcessingAlgorithm):
    # H-07: Constante VERSION centralizada — referenciar en todos los usos
    VERSION = 'v1.0.0'

    INPUT = 'INPUT'
    SAMPLING_AREA = 'SAMPLING_AREA'
    PROJECT_NAME = 'PROJECT_NAME'
    PERIMETER_DISTANCE = 'PERIMETER_DISTANCE'
    SELECT_INSIDE_BUFFER = 'SELECT_INSIDE_BUFFER'
    PLOT_RADIUS = 'PLOT_RADIUS'
    MIN_SAMPLE_DISTANCE = 'MIN_SAMPLE_DISTANCE'
    HILBERT_ORDER = 'HILBERT_ORDER'
    NUM_GROUPS = 'NUM_GROUPS'
    SAMPLE_SIZE = 'SAMPLE_SIZE'
    METHOD = 'METHOD'
    APPLY_PULLBACK = 'APPLY_PULLBACK'
    OUTPUT_ALL = 'OUTPUT_ALL'
    SAMPLE_WORD = 'SAMPLE_WORD'
    NUM_ITERATIONS = 'NUM_ITERATIONS'
    CALC_DIST_ALL_POINTS = 'CALC_DIST_ALL_POINTS'
    DISTANCE_UNITS = 'DISTANCE_UNITS'
    OUTPUT_HTML_REPORT = 'OUTPUT_HTML_REPORT'
    OPEN_REPORT = 'OPEN_REPORT'
    OUTPUT_REJECTED = 'OUTPUT_REJECTED'
    OUTPUT_HILBERT_PATH = 'OUTPUT_HILBERT_PATH'
    MANUAL_GRID_SPACING = 'MANUAL_GRID_SPACING'
    FIELD_GROUP_ID = "grupo_id_proc"

    # H-14: Tolerancia de deduplicación como constante nombrada (10 cm)
    DEDUP_TOLERANCE_M: float = 0.10

    METHOD_NAMES = {
        SamplingMethod.SYSTEMATIC_HILBERT: "Sistemático Hilbert",
        SamplingMethod.SIMPLE_RANDOM: "Aleatorio Simple",
        SamplingMethod.STRATIFIED_HILBERT: "Grupos Hilbert (Grupos 1D)",
        SamplingMethod.KMEANS_GROUPS: "Grupos K-Medias (Conglomerados 2D)"
    }
    
    METHOD_ABBREVIATIONS = {
        SamplingMethod.SYSTEMATIC_HILBERT: "SH",
        SamplingMethod.SIMPLE_RANDOM: "AL",
        SamplingMethod.STRATIFIED_HILBERT: "GH", 
        SamplingMethod.KMEANS_GROUPS: "KM"
    }
    
    UNIT_CONVERSIONS = {
        DistanceUnit.METERS: 1.0, DistanceUnit.KILOMETERS: 1000.0,
        DistanceUnit.FEET: 0.3048, DistanceUnit.MILES: 1609.34
    }
    
    FIELD_SAMPLE_ID = "muestra_id_proc"
    FIELD_HILBERT_IDX = "Hilbert_idx_proc"
    FIELD_ITER_NUM = "iter_num_proc"
    FIELD_METHOD = "metodo_proc"
    FIELD_CORRECTED = "correccion_proc"
    FIELD_DIST_BORDER = "dist_borde_m"
    FIELD_REJECTED_ITER = "rechazado_iter"
    FIELD_NEAREST_NEIGHBOR = "dist_vecino_m"
    FIELD_SAMPLE_IVMC = "IVMC_muestra_proc"
    FIELD_NEIGHBOR_CHECK = "Vecino"

    def __init__(self):
        super().__init__()
        self.has_as_meters = hasattr(QgsProcessingParameterDistance, 'asMeters')
        # Atributos de transferencia al hilo GUI
        self._sample_layers_pending: List[Tuple[str, str]] = []
        self._framework_layer_names: set = set()
        self._visibility_handler = None  # retención contra GC
        self._current_method: SamplingMethod = SamplingMethod.SYSTEMATIC_HILBERT
        self._has_rejected: bool = False
        self._report_gen = None        # ReportGenerator para postProcessAlgorithm
        self._html_path_pending = ''   # path HTML para postProcessAlgorithm
        self._html_open_pending = True # abrir reporte al terminar
        self._t0: float = 0.0          # tiempo de inicio global

    def tr(self, text: str) -> str:
        return QCoreApplication.translate('Processing', text)

    def createInstance(self) -> 'MuestreoEspacialPuntos':
        return MuestreoEspacialPuntos()

    def _connect_layers_added_signal(
        self,
        framework_names: set,
        pending: list,
        expected: int = 1
    ) -> None:
        """Conecta layerWasAdded con un QObject slot para garantizar
        ejecución en el hilo GUI mediante Qt.QueuedConnection.
        'expected' = número real de capas del framework que se generarán.
        """
        handler = _LayerVisibilityHandler(
            framework_names, pending, expected,
            report_gen=self._report_gen,
            t0=self._t0
        )
        self._visibility_handler = handler  # retención contra GC
        self._report_gen = None  # el handler toma la responsabilidad
        project = QgsProject.instance()
        project.layerWasAdded.connect(
            handler.on_layer_added,
            Qt.QueuedConnection
        )

    def postProcessAlgorithm(
        self,
        context: QgsProcessingContext,
        feedback: QgsProcessingFeedback
    ) -> Dict:
        """Invocado por QGIS en el hilo GUI DESPUÉS de cargar todas las capas.
        Calcula el tiempo real completo y reescribe el HTML final con ese dato.
        """
        import time as _time
        elapsed = _time.time() - self._t0 if self._t0 > 0 else 0.0

        if self._report_gen is not None:
            try:
                self._report_gen.rewrite_with_time(elapsed)
            except Exception as e:
                feedback.pushWarning(f"[postProcess] Error actualizando HTML: {e}")
            self._report_gen = None

        self._t0 = 0.0
        return {}

    def name(self) -> str:
        # Algorithm ID para QGIS Processing — único por complemento
        ver = self.VERSION.replace('.', '_').replace(' ', '_')
        return f'muestreo_espacial_puntos_{ver}'

    def displayName(self) -> str:
        return self.tr(
            'Muestreo Espacial de Puntos: Sistemático Hilbert, Grupos Hilbert, Aleatorio, K-Medias'
        )

    def shortHelpString(self) -> str:
        if HAS_SHAPELY:
            api_label = "contains_xy ufunc" if _shapely_use_v2 else "vectorized legacy"
            motor = f"Shapely {_SHAPELY_VERSION} ({api_label})"
        else:
            motor = "PreparedGeometry GEOS — estándar"
        return f"""<p><b>Notación numérica:</b> decimal con coma (3,14) · miles con espacio fino (1 234,56).</p>
<h3>Muestreo Espacial de Puntos {self.VERSION}</h3>
<p class="warn">&#9888; REQUISITO: La capa de puntos y el área de muestreo deben estar en un SRC proyectado en metros (ej: CRTM05 EPSG:8908, UTM). Un SRC geográfico (grados, EPSG:4326) cancela la ejecución con mensaje de error.</p>
<h3>Descripción General</h3>
<p>Algoritmo de muestreo espacial de puntos con cuatro métodos. Usa la curva de Hilbert para crear un ordenamiento espacial del marco muestral completo. Este ordenamiento es la base de los métodos Sistemático Hilbert y Grupos Hilbert, mejorando su distribución y representatividad. Para Aleatorio Simple y K-Medias, el orden no interviene en la selección, pero sí en el Reporte HTML para evaluar la calidad espacial del marco muestral.</p>
<p>El ordenamiento se basa en el rectángulo envolvente (Bounding Box) del área de estudio, alineado con los ejes X e Y. El algoritmo genera un número configurable de muestras iteradas cargadas como capas temporales en QGIS. Se descartan al cerrar el proyecto si no se guardan.</p>
<h3>Capas de Entrada</h3>
<p><b>Capa de puntos:</b> Define el marco muestral completo. El algoritmo verifica duplicados — su presencia puede comprometer el ordenamiento. Se recomienda usar el complemento <i>Crear malla de Puntos</i> para mallas hexagonales. Si activa "Objetos seleccionados solamente", el análisis opera únicamente sobre los puntos elegidos.</p>
<p><b>Área de muestreo (polígono):</b> Define los límites espaciales del estudio. Si activa "Objetos seleccionados solamente", solo se usan los polígonos elegidos.</p>
<h3>Métodos de Muestreo</h3>
<p><b>Sistemático Hilbert (SH):</b> Selección a intervalos regulares a lo largo de la curva. Cada iteración usa un inicio aleatorio, garantizando variabilidad. Produce patrones dispersos (IVMC &gt; 1,2). Genera la capa <i>Salida: Ruta Hilbert</i>.</p>
<p><b>Aleatorio Simple (AL):</b> Selección al azar sin reemplazo. Produce patrones aleatorios (IVMC ≈ 1,0). En mallas hexagonales regulares el resultado permanece disperso (IVMC &gt; 1,1).</p>
<p><b>Grupos Hilbert (GH):</b> Divide los puntos ordenados por la curva en N grupos lineales 1D y selecciona aleatoriamente dentro de cada uno. Genera también la capa <i>Salida: Ruta Hilbert</i>. Requiere definir <i>Número de grupos (k)</i>.</p>
<p><b>Grupos K-Medias (KM):</b> Agrupa los puntos en conglomerados 2D usando scikit-learn y extrae muestra proporcional por grupo. En mallas regulares el resultado es disperso (IVMC &gt; 1,1). Requiere <b>scikit-learn</b> instalado.</p>
<h3>Parámetros Clave</h3>
<p><b>Orden Hilbert (1-11):</b> Resolución de la cuadrícula. Un orden <i>n</i> crea una grilla de (2ⁿ × 2ⁿ) celdas. Valor por defecto (10) es adecuado para mallas de hasta 50 000 puntos o más. Para mallas muy densas (&gt;100 000 puntos) use orden 11. Para muestras pequeñas (&lt;500 puntos) valores bajos (6-8) son suficientes.</p>
<p><b>Distancia al perímetro:</b> Restringe el muestreo a una banda perimetral o al núcleo. Si <i>Muestrear DENTRO del margen</i> está activo, solo se usan puntos en esa franja; de lo contrario, se excluye.</p>
<p><b>Corrección de borde (Retracción):</b> Mueve puntos próximos al límite — cuando la distancia es menor al radio de parcela — hacia el interior, garantizando que las parcelas permanezcan íntegramente dentro del área. Aplicable a todos los métodos. Requiere <i>Radio de parcela</i>.</p>
<p><b>Distancia mínima entre puntos:</b> Garantiza separación mínima usando índice espacial O(n log n). Puntos rechazados van a la capa <i>Salida: Rechazados</i>. Para mallas hexagonales con espaciado <i>d</i>: use <i>d√3</i> para excluir el primer y segundo anillo de vecinos; use <i>d</i> para excluir solo el primero.</p>
<p><b>Distancia de malla conocida (g):</b> Si conoce el espaciado real del muestreo sistemático, ingréselo aquí para calcular el <i>Índice de Cobertura (h/g)</i> en el reporte.</p>
<p><b>Número de grupos (k):</b> Para métodos GH y KM. Si se deja en 0, el algoritmo calcula un valor automático basado en el tamaño de la muestra (fórmula: k = ⌈√n⌉, mínimo 2 grupos). Un solo grupo equivale a muestreo aleatorio.</p>
<h3>Métricas del Reporte HTML</h3>
<p><b>IVMC (Índice de Vecino Más Cercano):</b> La herramienta nativa de QGIS calcula el IVMC sobre el bounding box de la muestra — impreciso en polígonos irregulares y variable entre muestras. Este algoritmo aplica un factor de corrección en dos pasos: calcula el IVMC inicial, luego lo ajusta usando el área real del polígono de estudio, garantizando valores consistentes y comparables entre iteraciones.</p>
<p><b>Diagnóstico del Marco Muestral:</b></p>
<ul>
<li><b>Índice de Compacidad:</b> Evalúa la forma del polígono. Un valor bajo (&lt; 0,6) advierte sobre forma compleja o irregular con potenciales efectos de borde significativos.</li>
<li><b>CV (Coeficiente de Variación):</b> Evalúa el patrón de los puntos de entrada. CV bajo (&lt; 0,5) indica distribución uniforme; CV alto (&gt; 1,0) indica marco agrupado que podría requerir estratificación.</li>
<li><b>IVMC del marco:</b> IVMC &gt; 1,2 indica patrón disperso; IVMC &lt; 0,8 indica agrupamiento.</li>
</ul>
<p>Las recomendaciones del reporte son orientativas y deben adaptarse a cada proyecto.</p>
<h3>Rendimiento — Motor Shapely</h3>
<p>Shapely 2.x usa el API vectorizado (<code>contains_xy</code> ufunc) que opera sobre arrays NumPy completos, acelerando significativamente el filtrado de puntos en áreas complejas. El registro de mensajes indica cuál motor está activo.</p>
<h4>Compatibilidad por versión de QGIS</h4>
<p><b>QGIS 3.44 LTR y QGIS 4.0+:</b> Shapely 2.x está incluido en las instalaciones standalone (.msi) y OSGeo4W — motor Shapely activo sin instalación adicional.</p>
<p><b>QGIS 3.28 LTR – 3.40:</b> Shapely no está incluido por defecto. Instalar manualmente para activar el motor acelerado (opcional — el algoritmo funciona sin él).</p>
<p>Para instalar Shapely en versiones anteriores, en OSGeo4W Shell: <code>pip install shapely</code></p>
<h3>Instalar scikit-learn (para K-Medias)</h3>
<p><b>Opción 1 — OSGeo4W Shell (recomendada):</b> Búsquela en el menú Inicio de Windows como <i>OSGeo4W Shell</i> y ejecute:</p>
<p><code>pip install scikit-learn</code></p>
<p><b>Opción 2 — PowerShell (Windows):</b> PowerShell no tiene <code>pip</code> de QGIS en su PATH. Use la ruta completa al intérprete de QGIS (ajuste la versión según su instalación):</p>
<p><code>&amp; "C:\\Program Files\\QGIS 3.x\\bin\\python3.exe" -m pip install scikit-learn</code></p>
<p>Si instaló QGIS mediante OSGeo4W Network Installer, la ruta es:</p>
<p><code>&amp; "C:\\OSGeo4W\\bin\\python3.exe" -m pip install scikit-learn</code></p>
<p><code>pip3 install scikit-learn</code></p>
<p>Reinicie QGIS tras la instalación. Para verificar en la Consola de Python de QGIS (<code>Ctrl+Alt+P</code>): <code>import sklearn; print(sklearn.__version__)</code></p>
<div class="footer">Autor: Jorge Fallas (jfallas56@gmail.com) — {self.VERSION}</div>"""

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT, self.tr('Capa de puntos de entrada'),
            [_TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.SAMPLING_AREA, self.tr('Área de muestreo'),
            [_TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterString(
            self.PROJECT_NAME, self.tr('Nombre del Proyecto'),
            defaultValue='Proyecto de Muestreo'))
        if not self.has_as_meters:
            self.addParameter(QgsProcessingParameterEnum(
                self.DISTANCE_UNITS, self.tr('Unidades'),
                options=['Metros', 'Km', 'Pies', 'Millas'],
                defaultValue=DistanceUnit.METERS))
        dist_suffix = '' if self.has_as_meters else self.tr(' - Ver unidad arriba')
        self.addParameter(QgsProcessingParameterDistance(
            self.PERIMETER_DISTANCE,
            self.tr('Distancia al perímetro') + dist_suffix,
            parentParameterName=self.INPUT, optional=True, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            self.SELECT_INSIDE_BUFFER,
            self.tr('Muestrear DENTRO del margen'), defaultValue=False))
        self.addParameter(QgsProcessingParameterBoolean(
            self.APPLY_PULLBACK,
            self.tr('Aplicar corrección de borde (Retracción)'),
            defaultValue=False))
        self.addParameter(QgsProcessingParameterDistance(
            self.PLOT_RADIUS,
            self.tr("Radio de parcela (para Retracción)") + dist_suffix,
            parentParameterName=self.INPUT, defaultValue=0.0, optional=True))
        self.addParameter(QgsProcessingParameterDistance(
            self.MIN_SAMPLE_DISTANCE,
            self.tr('Distancia mínima entre puntos') + dist_suffix,
            parentParameterName=self.INPUT, optional=True, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterDistance(
            self.MANUAL_GRID_SPACING,
            self.tr('Distancia de malla conocida (g)') + dist_suffix,
            parentParameterName=self.INPUT, optional=True, defaultValue=0.0))
        method_opts = [
            self.METHOD_NAMES[m]
            for m in [
                SamplingMethod.SYSTEMATIC_HILBERT,
                SamplingMethod.SIMPLE_RANDOM,
                SamplingMethod.STRATIFIED_HILBERT
            ]
        ]
        if HAS_SKLEARN:
            method_opts.append(self.METHOD_NAMES[SamplingMethod.KMEANS_GROUPS])
        self.addParameter(QgsProcessingParameterEnum(
            self.METHOD, self.tr('Método de muestreo'),
            options=method_opts, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.SAMPLE_SIZE, self.tr('Tamaño de muestra'), defaultValue=10))
        self.addParameter(QgsProcessingParameterNumber(
            self.NUM_ITERATIONS, self.tr('Iteraciones'), defaultValue=1))
        self.addParameter(QgsProcessingParameterNumber(
            self.HILBERT_ORDER, self.tr('Orden Hilbert (1-11)'),
            defaultValue=10))
        self.addParameter(QgsProcessingParameterNumber(
            self.NUM_GROUPS, self.tr('Número de grupos (k)'), defaultValue=0))
        self.addParameter(QgsProcessingParameterString(
            self.SAMPLE_WORD, self.tr('Palabra clave salida'),
            defaultValue='muestra', optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_ALL, self.tr('Salida: Puntos filtrados'),
            optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_HILBERT_PATH, self.tr('Salida: Ruta Hilbert'),
            optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            self.CALC_DIST_ALL_POINTS,
            self.tr('Calcular distancia borde'), defaultValue=False))
        self.addParameter(QgsProcessingParameterFileDestination(
            self.OUTPUT_HTML_REPORT, self.tr('Reporte HTML'),
            fileFilter='HTML files (*.html)', optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            self.OPEN_REPORT, self.tr('Abrir reporte'), defaultValue=True))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_REJECTED, self.tr('Salida: Rechazados'),
            optional=True))

    # -----------------------------------------------------------------------
    # H-01 + H-03: checkParameterValues — validaciones previas al proceso
    # -----------------------------------------------------------------------
    def checkParameterValues(
        self, parameters: Dict, context: QgsProcessingContext
    ) -> Tuple[bool, str]:
        # --- H-01: SRC geográfico — cancelar antes de iniciar ---
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source:
            crs = source.sourceCrs()
            if crs.isGeographic():
                return False, (
                    f"[!] ERROR CRÍTICO — SRC geográfico detectado: "
                    f"{crs.authid()} (unidades en grados).\n"
                    "Este algoritmo requiere un SRC proyectado en metros "
                    "(ej: CRTM05 EPSG:8908, UTM zona correspondiente).\n"
                    "Reproyecte la capa de puntos antes de continuar."
                )
        area = self.parameterAsSource(parameters, self.SAMPLING_AREA, context)
        if area:
            crs_a = area.sourceCrs()
            if crs_a.isGeographic():
                return False, (
                    f"[!] ERROR CRÍTICO — SRC geográfico en área de muestreo: "
                    f"{crs_a.authid()} (unidades en grados).\n"
                    "Reproyecte el área de muestreo a un SRC proyectado en metros."
                )

        # --- H-03: Validaciones de parámetros numéricos ---
        order = self.parameterAsInt(parameters, self.HILBERT_ORDER, context)
        if not (1 <= order <= 11):
            return False, (
                f"[!] Orden Hilbert inválido: {order}. "
                "Debe estar entre 1 y 11."
            )

        sample_size = self.parameterAsInt(parameters, self.SAMPLE_SIZE, context)
        if sample_size < 1:
            return False, (
                "[!] Tamaño de muestra debe ser mayor o igual a 1."
            )

        num_iter = self.parameterAsInt(parameters, self.NUM_ITERATIONS, context)
        if num_iter < 1:
            return False, (
                "[!] El número de iteraciones debe ser mayor o igual a 1."
            )

        # --- H-03: K-Medias sin sklearn ---
        m_idx = self.parameterAsInt(parameters, self.METHOD, context)
        # Si sklearn no está y el índice corresponde a K-Medias (última opción)
        if not HAS_SKLEARN and HAS_SKLEARN is False:
            method_list = [
                SamplingMethod.SYSTEMATIC_HILBERT,
                SamplingMethod.SIMPLE_RANDOM,
                SamplingMethod.STRATIFIED_HILBERT,
                SamplingMethod.KMEANS_GROUPS
            ]
            if m_idx == 3:
                return False, (
                    "[!] El método 'Grupos K-Medias' requiere la biblioteca "
                    "'scikit-learn', que no está instalada.\n"
                    "Instale con: pip install scikit-learn (OSGeo4W Shell) o\n"
                    "& \"C:\\Program Files\\QGIS 3.x\\bin\\python3.exe\" -m pip install scikit-learn (PowerShell)\n"
                    "O seleccione otro método de muestreo."
                )

        return super().checkParameterValues(parameters, context)

    def processAlgorithm(self, parameters: Dict, context: QgsProcessingContext, feedback: QgsProcessingFeedback) -> Dict:
        import time as _time
        self._t0 = _time.time()
        _t0 = self._t0
        feedback.pushInfo(f"--- Iniciando {self.displayName()} ---")
        # H-04: Log de motor activo al inicio
        if HAS_SHAPELY:
            api_msg = "contains_xy ufunc" if _shapely_use_v2 else "vectorized.contains (legacy)"
            feedback.pushInfo(f"[Motor] Shapely {_SHAPELY_VERSION} disponible — {api_msg}.")
        else:
            feedback.pushInfo("[Motor] Shapely no disponible — usando PreparedGeometry GEOS (fallback automático).")
        temp_files = []
        try:
            layers, params, paths = self._setup(parameters, context)
            
            # H-01: Salvaguarda adicional para entornos programáticos donde
            # checkParameterValues() no se ejecuta automáticamente.
            if layers.get('source'):
                crs = layers['source'].sourceCrs()
                if crs.isGeographic():
                    raise QgsProcessingException(
                        f"[!] SRC geográfico detectado en tiempo de ejecución: "
                        f"{crs.authid()}. Use un SRC proyectado en metros."
                    )

            data = self._prepare_data(layers, params, context, feedback)
            if data is None:
                feedback.reportError(
                    "AVISO: No se generaron muestras — la capa de entrada está "
                    "vacía o no hay puntos en el área seleccionada."
                )
                return {}
                
            path_eval = self._evaluate_path(
                data['sorted_features'], data['area_for_filter'],
                params, data['distance_calculator'],
                data['duplicate_count'], feedback
            )
            raw_poly = parameters[self.SAMPLING_AREA]
            
            iter_res = self._run_iterations(data, params, parameters, context, feedback, raw_poly)
            temp_files.extend(iter_res.get('temp_files', []))
            
            res = self._create_outputs(
                iter_res, data, params, paths, parameters,
                context, feedback, path_eval, _t0
            )
            
            elapsed = _time.time() - _t0
            mins, secs = divmod(int(elapsed), 60)
            tiempo_str = f"{mins} min {secs:02d} s" if mins else f"{secs} s"
            feedback.pushInfo(f"Proceso completado en {tiempo_str}.")
            return res
        except Exception as e:
            import traceback
            feedback.reportError(traceback.format_exc())
            raise QgsProcessingException(f"Error: {e}")
        finally:
            # Liberar instancias retenidas del post-procesador de visibilidad
            _InvisibleLayerPostProcessor.clear_instances()
            gc.collect()
            for f in temp_files:
                try:
                    Path(f).unlink(missing_ok=True)
                except Exception:
                    pass

    def _setup(self, parameters, context):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        area = self.parameterAsSource(parameters, self.SAMPLING_AREA, context)
        if self.has_as_meters:
            p_dist, radius, min_d, grid_g = (self.parameterAsDistance(parameters, p, context) for p in [self.PERIMETER_DISTANCE, self.PLOT_RADIUS, self.MIN_SAMPLE_DISTANCE, self.MANUAL_GRID_SPACING])
        else:
            factor = self.UNIT_CONVERSIONS.get(DistanceUnit(self.parameterAsInt(parameters, self.DISTANCE_UNITS, context)), 1.0)
            p_dist, radius, min_d, grid_g = (self.parameterAsDouble(parameters, p, context) * factor for p in [self.PERIMETER_DISTANCE, self.PLOT_RADIUS, self.MIN_SAMPLE_DISTANCE, self.MANUAL_GRID_SPACING])
        m_idx = self.parameterAsInt(parameters, self.METHOD, context)
        m_enum = list(SamplingMethod)[m_idx] if m_idx < len(SamplingMethod) else SamplingMethod.SYSTEMATIC_HILBERT
        if m_enum == SamplingMethod.KMEANS_GROUPS and not HAS_SKLEARN: m_enum = SamplingMethod.SYSTEMATIC_HILBERT
        params = SamplingParameters(
            project_name=self.parameterAsString(parameters, self.PROJECT_NAME, context),
            perimeter_distance=p_dist, plot_radius=radius, min_sample_distance=min_d,
            hilbert_order=self.parameterAsInt(parameters, self.HILBERT_ORDER, context),
            num_groups=self.parameterAsInt(parameters, self.NUM_GROUPS, context),
            sample_size=self.parameterAsInt(parameters, self.SAMPLE_SIZE, context),
            num_iterations=self.parameterAsInt(parameters, self.NUM_ITERATIONS, context),
            method=m_enum,
            apply_pullback=self.parameterAsBoolean(parameters, self.APPLY_PULLBACK, context),
            select_inside_buffer=self.parameterAsBoolean(parameters, self.SELECT_INSIDE_BUFFER, context),
            calc_dist_all_points=self.parameterAsBoolean(parameters, self.CALC_DIST_ALL_POINTS, context),
            manual_grid_spacing=grid_g
        )
        html_raw = self.parameterAsString(parameters, self.OUTPUT_HTML_REPORT, context)
        if not html_raw or html_raw == 'TEMPORARY_OUTPUT': html_raw = QgsProcessingUtils.generateTempFilename('reporte.html')
        paths = {
            'output_folder': Path(QgsProcessingUtils.tempFolder()),
            'sample_word': self.parameterAsString(parameters, self.SAMPLE_WORD, context) or 'muestra',
            'html_report_path': html_raw,
            'open_report': self.parameterAsBoolean(parameters, self.OPEN_REPORT, context)
        }
        # Guardar para postProcessAlgorithm()
        self._html_path_pending = html_raw
        self._html_open_pending = self.parameterAsBoolean(parameters, self.OPEN_REPORT, context)
        return {'source': source, 'area': area}, params, paths

    def _prepare_data(self, layers, params, context, feedback):
        source, area = layers['source'], layers['area']
        if not source: return None
        
        count = source.featureCount()
        feedback.pushInfo(f"Entidades en capa de entrada: {count}")
        if count == 0: 
            feedback.reportError("AVISO: Capa de puntos vacía.")
            return None

        valid = [f.geometry() for f in area.getFeatures() if GeometryUtils.is_geometry_valid(f.geometry())]
        if not valid: 
            feedback.reportError("AVISO: Área de muestreo inválida.")
            return None
            
        orig_area = QgsGeometry.unaryUnion(valid)

        # H-02: Validar y reparar geometría post-unaryUnion (§24.1)
        if orig_area is None or orig_area.isEmpty():
            feedback.reportError("AVISO: unaryUnion() produjo geometría vacía. Verifique el área de muestreo.")
            return None
        if not orig_area.isGeosValid():
            feedback.pushWarning(
                "[H-02] Geometría inválida post-unaryUnion — aplicando makeValid(). "
                "Verifique la topología del área de muestreo."
            )
            orig_area = orig_area.makeValid()
            if orig_area is None or orig_area.isEmpty():
                feedback.reportError("ERROR: No se pudo reparar la geometría del área de muestreo.")
                return None
            feedback.pushInfo("[H-02] Geometría reparada con makeValid().")

        s_crs, a_crs = source.sourceCrs(), area.sourceCrs()
        
        filter_area = QgsGeometry(orig_area)
        if s_crs != a_crs:
            feedback.pushInfo(f"Reproyectando área de {a_crs.authid()} a {s_crs.authid()}...")
            transformer = QgsCoordinateTransform(a_crs, s_crs, QgsProject.instance())
            filter_area.transform(transformer)
        
        filter_geometry = QgsGeometry(filter_area)
        if params.perimeter_distance > 0:
            inner = filter_area.buffer(-params.perimeter_distance, 12)
            filter_geometry = inner if params.select_inside_buffer else filter_area.difference(inner)
        
        # H-04: Preparar geometría una sola vez para acelerar consultas masivas.
        # prepareGeometry() construye índice GEOS interno: ~500 pts/s → ~92 000 pts/s.
        _geom_prepared = False
        if HAS_SHAPELY:
            # Motor Shapely: convertir a shapely para ufunc vectorizado
            try:
                import shapely.wkb as _shwkb
                _sh_filter_geom = _shwkb.loads(
                    bytes(filter_geometry.asWkb()), hex=False
                )
                feedback.pushInfo("[Motor] Geometría de filtro convertida a Shapely.")
            except Exception as e:
                feedback.pushWarning(f"[Motor] Fallo conversión Shapely: {e} — usando PreparedGeometry.")
                _sh_filter_geom = None
        else:
            _sh_filter_geom = None

        # PreparedGeometry siempre como base (fallback y para conteos individuales)
        if hasattr(filter_geometry, 'prepareGeometry'):
            filter_geometry.prepareGeometry()
            _geom_prepared = True

        feedback.pushInfo("Escaneando puntos...")
        filtered_points_data = []
        seen_coords = set()
        duplicate_count = 0
        # H-14: Tolerancia de deduplicación como constante — 10 cm
        grid_multiplier = 1.0 / self.DEDUP_TOLERANCE_M
        
        request = QgsFeatureRequest()
        request.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        # H-04: Estrategia de filtrado según motor disponible
        if HAS_SHAPELY and _sh_filter_geom is not None:
            # Motor Shapely vectorizado: acumular coordenadas, filtrar en batch
            all_pts_raw: List[Tuple[int, float, float]] = []
            for feature in source.getFeatures(request):
                if feedback.isCanceled():
                    break
                if not feature.hasGeometry():
                    continue
                geom = feature.geometry()
                if geom.isEmpty():
                    continue
                p = geom.asPoint()
                all_pts_raw.append((feature.id(), p.x(), p.y()))

            feedback.pushInfo(f"Filtrando {len(all_pts_raw):,} puntos con Shapely vectorizado...")
            if all_pts_raw:
                try:
                    xs = _np.array([r[1] for r in all_pts_raw])
                    ys = _np.array([r[2] for r in all_pts_raw])
                    # Shapely >= 2.0: contains_xy(geom, xs, ys) — vectorizado sin deprecación
                    # Shapely 1.x: vectorized.contains(geom, xs, ys) — API legacy
                    if _shapely_use_v2:
                        mask = _shapely_contains_xy(_sh_filter_geom, xs, ys)
                    else:
                        mask = _shvec.contains(_sh_filter_geom, xs, ys)
                    for inside, (fid, px, py) in zip(mask, all_pts_raw):
                        if inside:
                            key = (int(px * grid_multiplier), int(py * grid_multiplier))
                            if key in seen_coords:
                                duplicate_count += 1
                            else:
                                seen_coords.add(key)
                                filtered_points_data.append(PointData(id=fid, x=px, y=py))
                except Exception as e:
                    feedback.pushWarning(
                        f"[Motor] Error en Shapely vectorizado: {e} — "
                        "reintentando con PreparedGeometry."
                    )
                    # Fallback a PreparedGeometry con los datos ya recopilados
                    for fid, px, py in all_pts_raw:
                        if filter_geometry.contains(QgsPointXY(px, py)):
                            key = (int(px * grid_multiplier), int(py * grid_multiplier))
                            if key in seen_coords:
                                duplicate_count += 1
                            else:
                                seen_coords.add(key)
                                filtered_points_data.append(PointData(id=fid, x=px, y=py))
            scanned = len(all_pts_raw)
        else:
            # Motor PreparedGeometry GEOS: contains(QgsPointXY) — sin objeto temporal
            scanned = 0
            feature_iterator = source.getFeatures(request)
            for feature in feature_iterator:
                scanned += 1
                if feedback.isCanceled():
                    return None
                if scanned % 50000 == 0:
                    feedback.pushInfo(f"Procesados {scanned:,}...")
                if not feature.hasGeometry():
                    continue
                geom = feature.geometry()
                if geom.isEmpty():
                    continue
                p = geom.asPoint()
                # H-04: contains(QgsPointXY) — sin crear QgsGeometry temporal por punto
                if filter_geometry.contains(p):
                    key = (int(p.x() * grid_multiplier), int(p.y() * grid_multiplier))
                    if key in seen_coords:
                        duplicate_count += 1
                    else:
                        seen_coords.add(key)
                        filtered_points_data.append(PointData(id=feature.id(), x=p.x(), y=p.y()))

        # H-04: Liberar caché PreparedGeometry (siempre, incluso si hubo excepción)
        if _geom_prepared and hasattr(filter_geometry, 'releaseCache'):
            try:
                filter_geometry.releaseCache()
            except Exception:
                pass
        
        if not filtered_points_data: 
             feedback.reportError(f"AVISO: 0 puntos dentro del área (de {scanned} escaneados).")
             return None
        
        feedback.pushInfo(f"Ordenando {len(filtered_points_data)} puntos válidos...")
        sorted_pts = self._sort_points(filtered_points_data, filter_area, params.hilbert_order)
        
        dist_area = QgsDistanceArea(); dist_area.setSourceCrs(s_crs, context.transformContext())
        return {'sorted_features': sorted_pts, 'filtered_features': filtered_points_data, 'original_area': orig_area, 'area_for_filter': filter_area, 'source_crs': s_crs, 'distance_calculator': DistanceCalculator(dist_area), 'duplicate_count': duplicate_count, 'initial_point_count': scanned, 'filtered_point_count': len(filtered_points_data), 'source_layer_for_hydration': source, 'area_layer': area}

    def _sort_points(self, points: List, area: QgsGeometry, order: int) -> List:
        """Ordena puntos según el índice de la curva de Hilbert.
        H-08: SortingMethod.PEANO eliminado — solo se implementa Hilbert."""
        bbox = area.boundingBox()
        min_x, min_y = bbox.xMinimum(), bbox.yMinimum()
        range_x, range_y = (bbox.width() or 1.0), (bbox.height() or 1.0)
        h_size = 1 << order
        for i, p in enumerate(points):
            p.original_index = i
            nx = int(math.floor((h_size - 1) * (p.x - min_x) / range_x))
            ny = int(math.floor((h_size - 1) * (p.y - min_y) / range_y))
            p.hilbert_index = HilbertCurveCalculator.xy_to_hilbert(max(0, min(nx, h_size - 1)), max(0, min(ny, h_size - 1)), order)
        points.sort(key=lambda x: x.hilbert_index)
        return points

    def _run_iterations(self, data, params, parameters, ctx, fb, poly_path):
        sorted_pts = data['sorted_features']; pool = data['filtered_features']
        res = {'generated_files_info': [], 'nna_results': [], 'all_rejected_points': [], 'temp_files_to_delete': [], 'hilbert_index_map': {p.id: i for i, p in enumerate(sorted_pts)}, 'strata': {}, 'kmeans_cluster_map': {}, 'kmeans_validation_data': [], 'k_efectivo': 0}
        
        strata = {}
        if params.method == SamplingMethod.STRATIFIED_HILBERT:
            k_efectivo = params.num_groups or self._suggest_groups(len(sorted_pts), params.sample_size)
            if params.num_groups == 0:
                fb.pushInfo(f"[Grupos Hilbert] k=0 → k automático calculado: {k_efectivo} grupos "
                            f"(muestra={params.sample_size}, total={len(sorted_pts)}).")
            strata = self._create_strata(len(sorted_pts), k_efectivo, fb)
            res['k_efectivo'] = k_efectivo
            res['strata'] = strata
        elif params.method == SamplingMethod.KMEANS_GROUPS:
            k_km = params.num_groups or self._suggest_groups(len(pool), params.sample_size)
            if params.num_groups == 0:
                fb.pushInfo(f"[K-Medias] k=0 → k automático calculado: {k_km} grupos "
                            f"(muestra={params.sample_size}, total={len(pool)}).")
            res['k_efectivo'] = k_km
            if HAS_SKLEARN: res['kmeans_cluster_map'] = self._get_kmeans_cluster_map(pool, k_km, fb)

        for i in range(params.num_iterations):
            if fb.isCanceled(): break
            it = i + 1
            
            # NOMBRES UNIFICADOS (FIX IndexError/AttributeError)
            if params.method == SamplingMethod.SIMPLE_RANDOM: 
                samp = self._select_random_sample(pool, params.sample_size, fb)
            elif params.method == SamplingMethod.STRATIFIED_HILBERT: 
                samp = self._select_stratified_sample(sorted_pts, params.sample_size, strata, fb, it)  # usa k_efectivo via strata
            elif params.method == SamplingMethod.KMEANS_GROUPS:
                samp = self._select_kmeans_sample(pool, params.sample_size, k_km, fb, it)
                if HAS_SKLEARN: res['kmeans_validation_data'].append(KMeansValidationData(it, self._get_kmeans_cluster_stats(pool, params.sample_size, k_km, fb, it)))
            else: 
                samp = self._select_systematic_sample(sorted_pts, params.sample_size, fb, it)
            
            if not samp: continue
            # P-04: deepcopy solo cuando se va a modificar la geometria (pullback)
            # Con __deepcopy__ de B-04 ya es seguro en ambos casos, pero evitar si no es necesario
            sc = [copy.deepcopy(p) for p in samp] if params.apply_pullback else list(samp)
            if params.apply_pullback: self._apply_pullback_correction(sc, data['original_area'], params.plot_radius, fb)
            if params.min_sample_distance > 0:
                sc, rej = self._enforce_minimum_distance(sc, params.min_sample_distance, data['distance_calculator'], params.sample_size, fb, it)
                res['all_rejected_points'].extend(rej)
            if not sc: continue
            
            nn_temp = self._calc_nn_temp(sc)
            contig = sum(1 for d in nn_temp.values() if d and (d < params.manual_grid_spacing or math.isclose(d, params.manual_grid_spacing, rel_tol=1e-4))) if params.manual_grid_spacing > 0 else 0
            nna, tmp = self._run_nna_for_sample(sc, data, parameters, ctx, fb, nn_temp, poly_path)
            if tmp: res['temp_files_to_delete'].append(tmp)  # None con NNA directa — no opera
            if nna:
                nna['contiguous_count'] = contig
                res['nna_results'].append(nna)
            # Guardar nna junto a la muestra (None si falló) para mantener
            # correspondencia 1:1 garantizada en _create_outputs (evita desorden).
            res['generated_files_info'].append((sc, nn_temp, nna))
        return res

    def _run_nna_for_sample(self, sample, data, parameters, ctx, fb, nn_map, poly_path):
        """Calcula IVMC directamente en Python — sin processing.run().

        Reemplaza 'qgis:nearestneighbouranalysis' eliminando:
          - Escritura de archivo .gpkg temporal por iteración
          - Invocación del framework de Processing (no thread-safe)
          - Lectura del HTML de resultado de QGIS

        Fórmula equivalente a la herramienta nativa de QGIS:
          d_obs  = media de distancias al vecino más cercano (ya en nn_map)
          d_esp  = 0.5 / sqrt(n / A)   donde A = área del polígono
          IVMC   = d_obs / d_esp
          Z      = (d_obs - d_esp) / (0.26136 / sqrt(n² / A))

        Factor de corrección de área: igual al aplicado antes con la herramienta
        nativa — ajusta por la relación entre área real y bounding box.
        """
        tmp = None  # ya no se genera archivo temporal
        n = len(sample)
        if n < 2:
            return None, tmp

        # Usar nn_map pre-calculado (O(n log n) con QgsSpatialIndex)
        dists = [d for d in nn_map.values() if d is not None and d > 0]
        if not dists:
            return None, tmp

        d_obs = sum(dists) / len(dists)
        area_geom = data.get('area_for_filter')
        if not area_geom or area_geom.isEmpty():
            return None, tmp

        A = area_geom.area()
        if A <= 0:
            return None, tmp

        # Distancia esperada para distribución aleatoria
        d_esp = 0.5 / math.sqrt(n / A)

        # IVMC base
        ivmc_base = d_obs / d_esp if d_esp > 0 else 0

        # Factor de corrección: área real vs bounding box (igual que versión nativa)
        ab = area_geom.boundingBox().area()
        if ab > 0:
            ivmc_corr = ivmc_base * math.sqrt(ab / A)
        else:
            ivmc_corr = ivmc_base

        # Z-score (para referencia, no se muestra en reporte actualmente)
        try:
            z = (d_obs - d_esp) / (0.26136 / math.sqrt(n * n / A))
        except (ZeroDivisionError, ValueError):
            z = 0.0

        # Conteo de polígonos muestreados (lógica sin cambios)
        polys_count = 0
        try:
            area_layer = data.get('area_layer')
            if area_layer:
                crs_transform = QgsCoordinateTransform(
                    data['source_crs'], area_layer.sourceCrs(), QgsProject.instance()
                )
                poly_feats = [f for f in area_layer.getFeatures() if f.hasGeometry()]
                poly_dict  = {f.id(): f for f in poly_feats}
                poly_index = QgsSpatialIndex()
                for f in poly_feats:
                    poly_index.addFeature(f)
                unique_polys_hit = set()
                for p in sample:
                    pt_geom_tr = QgsGeometry(p.get_geometry())
                    pt_geom_tr.transform(crs_transform)
                    for fid in poly_index.intersects(pt_geom_tr.boundingBox()):
                        poly_feat = poly_dict.get(fid)
                        if poly_feat and poly_feat.geometry().intersects(pt_geom_tr):
                            unique_polys_hit.add(fid)
                polys_count = len(unique_polys_hit)
        except Exception as e:
            if fb:
                fb.pushWarning(f"[NNA] Error contando polígonos: {e}")

        return {
            'indice_nn':        ivmc_corr,
            'n_puntos':         n,
            'polygons_hit_count': polys_count,
            'OBSERVED_MD':      d_obs,
            'EXPECTED_MD':      d_esp,
            'Z_SCORE':          z,
        }, tmp

    def _evaluate_path(self, points, area, params, calc, dups, fb):
        if len(points) < 2: return {}
        dists = [calc.calculate_distance(QgsPointXY(points[i-1].x, points[i-1].y), QgsPointXY(points[i].x, points[i].y)) for i in range(1, len(points))]
        tot = sum(dists); mean = tot / len(dists)
        h = math.sqrt(area.area() / len(points)) if len(points) > 0 else 0
        var = sum([(d-mean)**2 for d in dists])/(len(dists)-1) if len(dists)>1 else 0
        std = math.sqrt(var)
        
        h_idxs = set()
        bbox = area.boundingBox(); mx, my = bbox.xMinimum(), bbox.yMinimum(); rx, ry = (bbox.width() or 1), (bbox.height() or 1)
        h_size = 1 << params.hilbert_order
        for p in points:
            nx = int((h_size-1)*(p.x-mx)/rx); ny = int((h_size-1)*(p.y-my)/ry)
            h_idxs.add(HilbertCurveCalculator.xy_to_hilbert(nx, ny, params.hilbert_order))
        coll = len(points) - len(h_idxs)
        g = params.manual_grid_spacing if params.manual_grid_spacing > 0 else 0
        perim = area.length(); compac = (4*math.pi*area.area())/(perim**2) if perim > 0 else 0
        return {"efficiency_index_R": mean/h if h>0 else 0, "cv": std/mean if mean>0 else 0, "order": params.hilbert_order, "collisions": coll, "collision_ratio": (coll/len(points))*100, "area_km2": area.area()/1e6, "mesh_distance_h": h, "grid_spacing_g": g, "num_points": len(points), "total_length_km": tot/1000, "std_dev_norm": std/h if h>0 else 0, "length_norm": tot/(len(points)*h) if h>0 else 0, "indice_compacidad": compac, "coverage_index": h/g if g>0 else 0, "duplicate_count": dups}

    # --- FIX INDEX ERROR: Validación estricta de límites ---
    def _select_systematic_sample(self, pts, k, fb, it):
        n = len(pts)
        if n == 0: return []
        if k >= n: return list(pts)
        
        step = n / k
        start = random.uniform(0, step)
        indices = []
        
        for i in range(k):
            raw_val = start + i * step
            val = int(round(raw_val))
            
            # Clamp para asegurar que nunca exceda el último índice válido
            if val >= n: val = n - 1
            if val < 0: val = 0
            
            indices.append(val)
        
        # Eliminar duplicados y ordenar para mantener consistencia
        unique_indices = sorted(list(set(indices)))
        return [pts[i] for i in unique_indices]

    def _select_random_sample(self, pts, k, fb): return random.sample(pts, min(len(pts), k))

    def _select_stratified_sample(self, pts, k, strata, fb, it):
        s=[]; base, rem=divmod(min(k, len(pts)), len(strata)); keys=list(strata.keys()); random.shuffle(keys)
        for i, key in enumerate(keys):
            info=strata[key]; sub=pts[info.start:info.end]; sz=base+(1 if i<rem else 0)
            if sub: s.extend(random.sample(sub, min(sz, len(sub))))
        return s

    def _select_kmeans_sample(self, pts, k, g, fb, it):
        if not HAS_SKLEARN: return []
        # B-02: n_clusters real puede ser < g si len(pts) < g
        n_real = min(g, len(pts))
        coords=[[p.x,p.y] for p in pts]; km=KMeans(n_clusters=n_real, random_state=it).fit(coords)
        clus={i:[] for i in range(n_real)};
        for p,l in zip(pts, km.labels_): clus[l].append(p)
        # Eliminar clusters vacios antes de la asignacion proporcional
        clus = {k: v for k, v in clus.items() if v}
        s=[]; total=len(pts); alloc=[]
        for c_id, c_pts in clus.items():
            prop=len(c_pts)/total; ideal=prop*k; base=int(ideal); alloc.append((c_id, c_pts, base, ideal-base))
        alloc.sort(key=lambda x:x[3], reverse=True); rem=k-sum(x[2] for x in alloc)
        for c_id, c_pts, base, _ in alloc:
            sz = base + (1 if rem > 0 else 0); rem -= 1
            if sz > 0: s.extend(random.sample(c_pts, min(len(c_pts), sz)))
        return s

    def _create_strata(self, tot, num, fb):
        s={}; base, rem=divmod(tot, num); st=0
        for i in range(num):
            sz=base+1 if i<rem else base
            if sz>0: s[i]=StratumInfo(st, st+sz, sz); st+=sz
        return s
    
    # H-09: _get_kmeans_cluster_map implementado — retorna dict fid→cluster_id
    def _get_kmeans_cluster_map(self, pts: List, g: int, fb: QgsProcessingFeedback) -> Dict:
        """Asigna a cada punto su cluster K-Medias para poblar FIELD_GROUP_ID."""
        if not HAS_SKLEARN or not pts:
            return {}
        try:
            n_real = min(g, len(pts))
            coords = [[p.x, p.y] for p in pts]
            km = KMeans(n_clusters=n_real, random_state=42, n_init='auto').fit(coords)
            return {p.id: int(lbl) for p, lbl in zip(pts, km.labels_)}
        except Exception as e:
            fb.pushWarning(f"[KMeans cluster map] Error: {e}")
            return {}

    # H-09: _get_kmeans_cluster_stats implementado — estadísticas por cluster
    def _get_kmeans_cluster_stats(
        self, pts: List, k: int, g: int, fb: QgsProcessingFeedback, it: int
    ) -> List[Dict]:
        """Calcula estadísticas de tamaño y dispersión por cluster K-Medias."""
        if not HAS_SKLEARN or not pts:
            return []
        try:
            n_real = min(g, len(pts))
            coords = [[p.x, p.y] for p in pts]
            km = KMeans(n_clusters=n_real, random_state=it, n_init='auto').fit(coords)
            stats: List[Dict] = []
            for cid in range(n_real):
                members = [pts[i] for i, lbl in enumerate(km.labels_) if lbl == cid]
                if not members:
                    continue
                cx, cy = km.cluster_centers_[cid]
                dists = [math.hypot(p.x - cx, p.y - cy) for p in members]
                stats.append({
                    'cluster_id': cid,
                    'size': len(members),
                    'center_x': cx,
                    'center_y': cy,
                    'mean_dist_to_center': sum(dists) / len(dists),
                    'max_dist_to_center': max(dists),
                })
            return stats
        except Exception as e:
            fb.pushWarning(f"[KMeans cluster stats it={it}] Error: {e}")
            return []


    def _apply_pullback_correction(self, pts: List, area: QgsGeometry, rad: float, fb: QgsProcessingFeedback) -> None:
        """Retrae puntos fuera del búfer interior al punto más cercano dentro de él."""
        if rad <= 0:
            return
        try:
            safe = area.buffer(-rad, 12)
        except Exception as e:
            fb.pushWarning(f'[Pullback] Error calculando búfer interior (-{rad} m): {e}')
            return
        if safe is None or safe.isEmpty():
            fb.pushWarning(f'[Pullback] Búfer interior vacío con radio {rad} m — área demasiado pequeña.')
            return
        for p in pts:
            if not safe.contains(p.get_geometry()):
                try:
                    p.set_geometry(safe.nearestPoint(p.get_geometry()))
                    p.was_corrected = True
                except Exception as e:
                    # H-06: registrar fallo de corrección — no silencioso
                    fb.pushWarning(
                        f'[Pullback] No se pudo retraer punto id={p.id} '
                        f'({p.x:.4f}, {p.y:.4f}): {e}'
                    )

    def _enforce_minimum_distance(self, pts, d, calc, k, fb, it):
        """P-05: O(n log n) con QgsSpatialIndex en vez de O(n^2)."""
        random.shuffle(pts)
        k_list, r_list = [], []
        idx = QgsSpatialIndex()
        accepted_geoms = {}
        for p in pts:
            pt_geom = p.get_geometry()
            # Buscar vecino mas cercano en los ya aceptados
            neighbors = idx.nearestNeighbor(pt_geom, 1)
            too_close = False
            if neighbors:
                nn_geom = accepted_geoms.get(neighbors[0])
                if nn_geom is not None:
                    dist = nn_geom.distance(pt_geom)
                    too_close = dist < d
            if not too_close:
                f = QgsFeature(); f.setId(p.id); f.setGeometry(pt_geom)
                idx.addFeature(f)
                accepted_geoms[p.id] = pt_geom
                k_list.append(p)
            else:
                p.rejected_in_iteration = it
                r_list.append(p)
        return k_list, r_list

    def _calc_nn_temp(self, pts):
        """P-01: O(n log n) con QgsSpatialIndex en vez de O(n^2)."""
        if not pts: return {}
        if len(pts) == 1: return {pts[0].id: None}
        # Construir indice espacial
        idx = QgsSpatialIndex()
        feats_tmp = {}
        for p in pts:
            f = QgsFeature(); f.setId(p.id); f.setGeometry(p.get_geometry())
            idx.addFeature(f)
            feats_tmp[p.id] = p.get_geometry()
        m = {}
        for p in pts:
            # nearestNeighbor devuelve lista; pedir 2 porque el primero es el propio punto
            neighbors = idx.nearestNeighbor(p.get_geometry(), 2)
            nn_id = next((nid for nid in neighbors if nid != p.id), None)
            if nn_id is not None and nn_id in feats_tmp:
                m[p.id] = feats_tmp[nn_id].distance(p.get_geometry())
            else:
                m[p.id] = None
        return m
    
    def _suggest_groups(self, tot: int, samp: int) -> int:
        """Calcula el número de grupos automático cuando k=0.

        Fórmula: k = ceil(sqrt(sample_size)), acotado por sqrt(total_puntos).
        Garantiza al menos ceil(sqrt(n)) puntos por grupo y nunca más grupos
        que raíz del total (evita grupos vacíos).
        Mínimo: 2 grupos (un solo grupo equivale a muestreo aleatorio).
        Ejemplo: samp=100, tot=14747 → k=min(10, 121)=10
                 samp=10,  tot=14747 → k=min(4,  121)=4
        """
        import math
        base = max(2, math.ceil(math.sqrt(samp)))
        return min(base, max(2, int(math.sqrt(tot))))

    def _create_outputs(self, res, data, params, paths, parameters, ctx, fb, path_eval, t0: float = 0.0):
        source = data['source_layer_for_hydration']; fields = self._fields(source.fields())
        idx_map = res['hilbert_index_map']
        boundary = GeometryUtils.get_polygon_boundary_as_line(data['area_for_filter'], fb)
        
        all_ids = set()
        for s, *_ in res['generated_files_info']: all_ids.update(p.id for p in s)
        if res['all_rejected_points']: all_ids.update(p.id for p in res['all_rejected_points'])
        if not all_ids: return {}
        req = QgsFeatureRequest().setFilterFids(list(all_ids))
        fmap = {f.id(): QgsFeature(f) for f in source.getFeatures(req)}

        # Generar archivos de muestra en orden 1..N.
        # Se acumulan en self._sample_layers_pending para que el post-procesador
        # de OUTPUT_ALL los cargue en orden desde el hilo GUI.
        # NO se registran en el contexto — se cargan directamente desde GUI.
        self._sample_layers_pending = []  # reset para nueva ejecución
        self._current_method = params.method
        self._has_rejected = bool(res.get('all_rejected_points'))

        for i, entry in enumerate(res['generated_files_info']):
            samp, nn = entry[0], entry[1]
            nna = entry[2] if len(entry) > 2 else None
            it = i + 1

            ivmc = nna.get('indice_nn') if nna else None
            cont = nna.get('contiguous_count') if nna else None

            nm = f"{self.METHOD_ABBREVIATIONS[params.method]}_{paths['sample_word']}_{it:02d}"
            if params.apply_pullback: nm += "_RETR"
            if cont is not None: nm += f"_Cont_{cont}"
            if ivmc: nm += f"_IVMC_{f'{ivmc:.3f}'.replace('.',',')}"

            if nna is not None:
                nna['muestra'] = nm

            hyd_samp = []
            for p in samp:
                if p.id in fmap:
                    f = QgsFeature(fmap[p.id]); f.setGeometry(p.get_geometry())
                    hyd_samp.append(SamplePoint(f, QgsGeometry(f.geometry()), p.was_corrected, p.correction_distance, None))

            tmp = QgsProcessingUtils.generateTempFilename(f"{nm}.gpkg")
            ok = self._write_sample_file(
                hyd_samp, Path(tmp), fields, idx_map, it, params.method,
                boundary, data['source_crs'], ctx, fb, params, ivmc, nn_map=nn
            )
            if ok:
                # Acumular en orden 1..N — el post-procesador los carga en orden
                self._sample_layers_pending.append((tmp, nm))

        if paths['html_report_path']:
            gen = ReportGenerator(res['nna_results'], params.project_name, params, self.METHOD_NAMES, path_eval, data['initial_point_count'], data['filtered_point_count'], res['kmeans_validation_data'], t0=t0)
            gen.write(paths['html_report_path'], fb, paths['open_report'])
            self._report_gen = gen  # retener para postProcessAlgorithm()
        
        result = {}
        
        if res['all_rejected_points']:
            hyd_rej = []
            for p in res['all_rejected_points']:
                if p.id in fmap:
                    f = QgsFeature(fmap[p.id]); f.setGeometry(p.get_geometry())
                    hyd_rej.append(SamplePoint(f, QgsGeometry(f.geometry()), rejected_in_iteration=p.rejected_in_iteration))
            self._write_rejected(hyd_rej, source, parameters, ctx, result)

        self._create_all_points(parameters, ctx, data['sorted_features'], source, boundary, params.calc_dist_all_points, data['source_crs'], result, res['strata'], res['kmeans_cluster_map'])
        
        # Ruta Hilbert: Sistemático Hilbert Y Grupos Hilbert
        # (ambos usan sorted_pts ordenado por curva de Hilbert)
        if params.method in (SamplingMethod.SYSTEMATIC_HILBERT, SamplingMethod.STRATIFIED_HILBERT):
            self._create_hilbert_path(parameters, params, ctx, data['sorted_features'], data['source_crs'], path_eval, result)

        # Asignar post-procesadores a TODOS los sinks de una vez,
        # después de que el framework haya registrado todos sus LayerDetails.
        self._assign_post_processors(ctx)

        return result

    # -----------------------------------------------------------------------
    # H-05: Carga invisible — usa QgsProcessingLayerPostProcessorInterface
    # para operar en el hilo principal de GUI, nunca desde el hilo de fondo.
    # -----------------------------------------------------------------------
    def _load_invisible_layer_sink(
        self,
        path: str,
        name: str,
        context: QgsProcessingContext,
        feedback: Optional[QgsProcessingFeedback]
    ) -> None:
        """Asigna _InvisibleLayerPostProcessor a una capa sink del framework.

        El 'path' recibido puede ser el ID del sink (ej. "memory:...") o el
        path real del archivo. El framework registra los LayerDetails usando
        el path real. Se intenta la búsqueda directa primero; si falla,
        se itera sobre todos los LayerDetails buscando por nombre de capa.
        """
        try:
            processor = _InvisibleLayerPostProcessor.create()
            # Intento 1: búsqueda directa por path/ID
            if context.willLoadLayerOnCompletion(path):
                context.layerToLoadOnCompletionDetails(path).setPostProcessor(processor)
                return
            # Intento 2: iterar todos los LayerDetails y buscar por nombre
            # Necesario cuando 'path' es el ID del sink y el framework
            # registró los detalles con el path real del archivo temporal.
            try:
                for registered_path, details in context.layersToLoadOnCompletion().items():
                    if details.name == name:
                        details.setPostProcessor(processor)
                        return
            except Exception:
                pass
            # Intento 3: si el context tiene el path como clave directa
            # (algunas versiones de QGIS usan el ID del sink como clave)
            try:
                all_keys = list(context.layersToLoadOnCompletion().keys())
                if all_keys:
                    # Buscar por coincidencia parcial del nombre en la clave
                    for k in all_keys:
                        if name.lower() in k.lower() or k == path:
                            context.layersToLoadOnCompletion()[k].setPostProcessor(processor)
                            return
            except Exception:
                pass
            if feedback:
                feedback.pushWarning(
                    f"[H-05] No se encontró LayerDetails para '{name}' — "
                    "la capa podría cargarse visible."
                )
        except Exception as e:
            if feedback:
                feedback.pushWarning(f"[H-05] Error asignando post-procesador para '{name}': {e}.")

    def _fields(self, src: QgsFields) -> QgsFields:
        f = QgsFields(src)
        for n in [self.FIELD_SAMPLE_ID, self.FIELD_HILBERT_IDX, self.FIELD_ITER_NUM]:
            f.append(QgsField(n, _INT_TYPE))
        for n in [self.FIELD_METHOD, self.FIELD_CORRECTED, self.FIELD_NEIGHBOR_CHECK]:
            f.append(QgsField(n, _STR_TYPE))
        for n in [self.FIELD_NEAREST_NEIGHBOR, self.FIELD_SAMPLE_IVMC, self.FIELD_DIST_BORDER]:
            f.append(QgsField(n, _DBL_TYPE))
        return f
        
    def _create_output_fields(self, src, dist):
        """Campos para la capa 'Todos los puntos filtrados'.
        Solo incluye FIELD_HILBERT_IDX (único campo que se rellena en esta capa).
        Los campos de muestra (ID, iteración, método, IVMC, etc.) son exclusivos
        de las capas de muestra individuales y se omiten aquí para evitar NULLs.
        """
        f = QgsFields(src)
        # Solo agregar el índice Hilbert — campo relevante para todos los puntos
        f.append(QgsField(self.FIELD_HILBERT_IDX, _INT_TYPE))
        # dist_borde solo si se solicitó y aplica
        if dist:
            f.append(QgsField(self.FIELD_DIST_BORDER, _DBL_TYPE))
        return f

    def _write_sample_file(self, pts, path, flds, imap, it, method, boundary, crs, ctx, fb, params, ivmc, nn_map=None):
        opts = QgsVectorFileWriter.SaveVectorOptions(); opts.driverName="GPKG"; opts.layerName=path.stem
        w = QgsVectorFileWriter.create(str(path), flds, QgsWkbTypes.Point, crs, ctx.transformContext(), opts)
        if w.hasError() != QgsVectorFileWriter.NoError: return False
        
        pts.sort(key=lambda sp: imap.get(sp.feature.id(), float('inf')))
        # P-02: no recalcular O(n^2); usar nn_map pre-calculado en _calc_nn_temp
        geoms = [p.feature.geometry() for p in pts]
        col_map = {field.name(): i for i, field in enumerate(flds)}
        
        for i, p in enumerate(pts):
            # P-02: distancia al vecino mas cercano del mapa pre-calculado
            if nn_map is not None:
                final_nn = nn_map.get(p.feature.id()) or 0
            else:
                # Fallback O(n^2) por compatibilidad si no se pasa nn_map
                curr = geoms[i]; min_d = float('inf')
                for j, other in enumerate(geoms):
                    if i == j: continue
                    d = curr.distance(other)
                    if d < min_d: min_d = d
                final_nn = min_d if min_d != float('inf') else 0
            
            n_check = "NULL"
            if params.manual_grid_spacing > 0:
                 is_neighbor = final_nn < params.manual_grid_spacing or math.isclose(final_nn, params.manual_grid_spacing, rel_tol=1e-4)
                 n_check = "Verdadero" if is_neighbor else "Falso"
            elif params.manual_grid_spacing == 0: n_check = "N/A (g=0)"

            dist_border = 0
            if boundary and p.feature.geometry():
                 try: dist_border = float(p.feature.geometry().distance(boundary))
                 except Exception: pass
            
            attrs = [None] * len(flds)
            for field in p.feature.fields():
                 idx = col_map.get(field.name())
                 if idx is not None: attrs[idx] = p.feature.attribute(field.name())
            
            if self.FIELD_SAMPLE_ID in col_map: attrs[col_map[self.FIELD_SAMPLE_ID]] = i + 1
            if self.FIELD_HILBERT_IDX in col_map: attrs[col_map[self.FIELD_HILBERT_IDX]] = imap.get(p.feature.id(), -1)
            if self.FIELD_ITER_NUM in col_map: attrs[col_map[self.FIELD_ITER_NUM]] = it
            if self.FIELD_METHOD in col_map: attrs[col_map[self.FIELD_METHOD]] = self.METHOD_ABBREVIATIONS[method]
            if self.FIELD_CORRECTED in col_map: attrs[col_map[self.FIELD_CORRECTED]] = "Verdadero" if p.was_corrected else "Falso"
            if self.FIELD_NEAREST_NEIGHBOR in col_map: attrs[col_map[self.FIELD_NEAREST_NEIGHBOR]] = round(final_nn, 2) if final_nn else 0.0
            if self.FIELD_NEIGHBOR_CHECK in col_map: attrs[col_map[self.FIELD_NEIGHBOR_CHECK]] = n_check
            if self.FIELD_SAMPLE_IVMC in col_map and ivmc: attrs[col_map[self.FIELD_SAMPLE_IVMC]] = round(float(ivmc), 4)
            if self.FIELD_DIST_BORDER in col_map: attrs[col_map[self.FIELD_DIST_BORDER]] = round(float(dist_border), 2)

            f = QgsFeature(); f.setGeometry(p.feature.geometry()); f.setAttributes(attrs)
            w.addFeature(f)
        del w
        return True

    def _write_rejected(self, pts, src, params, ctx, res):
        # H-11: usar alias de tipo compatible Qt5/Qt6
        flds = QgsFields(src.fields())
        flds.append(QgsField(self.FIELD_REJECTED_ITER, _INT_TYPE))
        sink, did = self.parameterAsSink(
            params, self.OUTPUT_REJECTED, ctx, flds, src.wkbType(), src.sourceCrs()
        )
        if sink:
            for p in pts:
                f = QgsFeature(flds)
                f.setGeometry(p.feature.geometry())
                for a in src.fields():
                    f.setAttribute(a.name(), p.feature.attribute(a.name()))
                f.setAttribute(self.FIELD_REJECTED_ITER, p.rejected_in_iteration)
                sink.addFeature(f, QgsFeatureSink.FastInsert)
            res[self.OUTPUT_REJECTED] = did
            # H-10: Se elimina la llamada incorrecta a _load_invisible_layer_sink(did, ...)
            # 'did' es un ID de sink, no una ruta de archivo. La capa de rechazados
            # se carga por el mecanismo estándar de QgsProcessingContext.

    def _assign_post_processors(self, ctx: QgsProcessingContext) -> None:
        """Conecta layerWasAdded para ocultar capas en el hilo GUI.

        IMPORTANTE: framework_names es dinámico — solo incluye las capas
        que realmente se van a generar según el método y los parámetros.
        Si se incluye 'Ruta Hilbert' pero el método es Aleatorio, el
        handler espera 3 capas pero solo llegan 2 → _load_samples nunca
        se ejecuta → las muestras no se cargan.
        """
        # Capas que siempre se generan
        framework_names = {
            self.tr('Salida: Puntos filtrados'),
            'Salida: Puntos filtrados',
        }
        # OUTPUT_REJECTED solo si hay puntos rechazados (min_distance > 0)
        # Pero como no sabemos si habrá rechazados, siempre incluirlo como
        # OPCIONAL — el handler no debe bloquearse esperando esta capa.
        # OUTPUT_HILBERT_PATH solo si el método es Sistemático Hilbert
        hilbert_names = set()
        if self._current_method in (SamplingMethod.SYSTEMATIC_HILBERT, SamplingMethod.STRATIFIED_HILBERT):
            hilbert_names = {
                self.tr('Salida: Ruta Hilbert'),
                'Salida: Ruta Hilbert',
            }
            framework_names |= hilbert_names

        rejected_names = {
            self.tr('Salida: Rechazados'),
            'Salida: Rechazados',
        }
        framework_names |= rejected_names

        self._framework_layer_names = framework_names
        pending = list(self._sample_layers_pending)

        # Contar cuántas capas realmente se van a generar
        expected = 1  # OUTPUT_ALL siempre
        if self._current_method in (SamplingMethod.SYSTEMATIC_HILBERT, SamplingMethod.STRATIFIED_HILBERT):
            expected += 1  # OUTPUT_HILBERT_PATH
        if self._has_rejected:
            expected += 1  # OUTPUT_REJECTED

        try:
            self._connect_layers_added_signal(framework_names, pending, expected)
        except Exception:
            pass


    def _create_all_points(self, params, ctx, pts, src, bnd, dist, crs, res, strata, kmap):
        flds = self._create_output_fields(src.fields(), dist)
        # H-11: usar alias de tipo compatible Qt5/Qt6
        if strata or kmap:
            flds.append(QgsField(self.FIELD_GROUP_ID, _INT_TYPE))
        sink, did = self.parameterAsSink(params, self.OUTPUT_ALL, ctx, flds, QgsWkbTypes.Point, crs)
        if sink:
            imap = {p.id: i for i, p in enumerate(pts)}
            req = QgsFeatureRequest().setFilterFids([p.id for p in pts])
            col_map = {field.name(): idx for idx, field in enumerate(flds)}
            has_border = self.FIELD_DIST_BORDER in col_map
            for f in src.getFeatures(req):
                nf = QgsFeature(flds)
                nf.setGeometry(f.geometry())
                for a in src.fields():
                    if a.name() in col_map:
                        nf.setAttribute(a.name(), f.attribute(a.name()))
                nf.setAttribute(self.FIELD_HILBERT_IDX, imap.get(f.id(), -1))
                if has_border and bnd and f.hasGeometry():
                    try:
                        nf.setAttribute(self.FIELD_DIST_BORDER,
                                        round(float(f.geometry().distance(bnd)), 2))
                    except Exception:
                        pass
                sink.addFeature(nf, QgsFeatureSink.FastInsert)
            res[self.OUTPUT_ALL] = did
            # El post-procesador se asignará en _assign_post_processors()
            # llamado desde _create_outputs() después de crear todos los sinks.

    def _create_hilbert_path(self, parameters, params, ctx, pts, crs, path_eval, res):
        """parameters = dict del framework (para parameterAsSink)
           params     = SamplingParameters dataclass (para hilbert_order, method)
        """
        # Campos descriptivos de la Ruta Hilbert
        flds = QgsFields()
        flds.append(QgsField("id",           _INT_TYPE))
        flds.append(QgsField("descripcion",  _STR_TYPE))
        flds.append(QgsField("orden",        _INT_TYPE))
        flds.append(QgsField("n_puntos",     _INT_TYPE))
        flds.append(QgsField("metodo",       _STR_TYPE))
        sink, did = self.parameterAsSink(
            parameters, self.OUTPUT_HILBERT_PATH, ctx, flds, QgsWkbTypes.LineString, crs
        )
        if sink:
            line = QgsGeometry.fromPolylineXY([QgsPointXY(p.x, p.y) for p in pts])
            f = QgsFeature(flds)
            f.setGeometry(line)
            f.setAttribute(0, 1)
            f.setAttribute(1, "Curva de Hilbert — ordenamiento espacial del marco muestral")
            f.setAttribute(2, params.hilbert_order)
            f.setAttribute(3, len(pts))
            f.setAttribute(4, self.METHOD_NAMES.get(params.method, ""))
            sink.addFeature(f, QgsFeatureSink.FastInsert)
            res[self.OUTPUT_HILBERT_PATH] = did
            # El post-procesador se asignará en _assign_post_processors()
            # llamado desde _create_outputs() después de crear todos los sinks.

# ---------------------------------------------------------------------------
# H-05: Post-procesador para cargar capas como invisibles en el hilo de GUI.
# Se instancia FUERA de MuestreoEspacialPuntos para que QGIS pueda encontrarlo
# mediante el sistema de plugins de Qt.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Post-procesadores de capas — se ejecutan en el hilo principal de Qt.
# ---------------------------------------------------------------------------

class _InvisibleLayerPostProcessor(QgsProcessingLayerPostProcessorInterface):
    """Oculta una capa sink del framework (OUTPUT_ALL, OUTPUT_HILBERT_PATH, etc.)
    al finalizar el proceso. Se ejecuta en el hilo GUI — seguro para Qt.
    Cada instancia es independiente para evitar GC prematuro.
    """
    _instances: list = []

    def __init__(self):
        super().__init__()

    @classmethod
    def create(cls) -> '_InvisibleLayerPostProcessor':
        inst = cls()
        cls._instances.append(inst)
        return inst

    @classmethod
    def clear_instances(cls) -> None:
        cls._instances.clear()

    def postProcessLayer(self, layer, context, feedback):
        if not layer or not layer.isValid():
            return
        try:
            root = QgsProject.instance().layerTreeRoot()
            node = root.findLayer(layer.id())
            if node:
                node.setItemVisibilityChecked(False)
        except Exception as e:
            if feedback:
                feedback.pushWarning(f"[PostProcessor] No se pudo ocultar '{layer.name()}': {e}")
        finally:
            try:
                _InvisibleLayerPostProcessor._instances.remove(self)
            except ValueError:
                pass


class _SampleLayerLoader(QgsProcessingLayerPostProcessorInterface):
    """Post-procesador especial que se adjunta a OUTPUT_ALL (o a una capa
    centinela vacía si OUTPUT_ALL no está activo).

    Al ser invocado por QGIS en el hilo GUI, carga TODAS las capas de muestra
    acumuladas en self._pending en orden 1..N, invisibles, insertándolas en
    el árbol con la posición correcta.

    Recibe la lista de (path, name) en orden ascendente (muestra_01 primero).
    Para que muestra_01 quede ARRIBA en el árbol, se insertan en orden INVERSO
    usando insertChildNode(0, ...) — cada inserción empuja las anteriores abajo.
    """
    # Retención contra GC — una única instancia activa a la vez
    _instance = None

    def __init__(self, pending: List[Tuple[str, str]]):
        """
        Args:
            pending: lista de (path_gpkg, nombre_capa) en orden 1..N.
        """
        super().__init__()
        self._pending = pending  # [(path, name), ...] orden 1..N

    @classmethod
    def create(cls, pending: List[Tuple[str, str]]) -> '_SampleLayerLoader':
        cls._instance = cls(pending)
        return cls._instance

    def postProcessLayer(self, layer, context, feedback):
        """Carga todas las capas pendientes en orden 1..N, invisibles.

        Algoritmo de orden en árbol:
          - Se itera en orden INVERSO (N→1) usando insertChildNode(0, nodo).
          - Cada inserción en posición 0 empuja las anteriores hacia abajo.
          - Resultado final: muestra_01 en tope, muestra_N al fondo.
        """
        from qgis.core import QgsLayerTreeLayer

        root = QgsProject.instance().layerTreeRoot()

        # Ocultar también la capa centinela (OUTPUT_ALL / Todos Puntos)
        if layer and layer.isValid():
            try:
                sentinel_node = root.findLayer(layer.id())
                if sentinel_node:
                    sentinel_node.setItemVisibilityChecked(False)
            except Exception:
                pass

        # Cargar capas de muestra en orden inverso (N..1)
        loaded_count = 0
        for path, name in reversed(self._pending):
            try:
                lyr = QgsVectorLayer(path, name, "ogr")
                if not lyr.isValid():
                    if feedback:
                        feedback.pushWarning(f"[SampleLoader] Capa inválida — omitida: {name}")
                    continue
                # Añadir al proyecto sin insertar en árbol (addToLegend=False)
                QgsProject.instance().addMapLayer(lyr, False)
                # Crear nodo de árbol y colocarlo en posición 0
                tree_node = QgsLayerTreeLayer(lyr)
                tree_node.setItemVisibilityChecked(False)
                root.insertChildNode(0, tree_node)
                loaded_count += 1
            except Exception as e:
                if feedback:
                    feedback.pushWarning(f"[SampleLoader] Error cargando '{name}': {e}")

        if feedback:
            feedback.pushInfo(
                f"[SampleLoader] {loaded_count}/{len(self._pending)} capas de muestra "
                "cargadas — ocultas, orden 1..N."
            )
        _SampleLayerLoader._instance = None  # liberar retención


# ---------------------------------------------------------------------------
# Generador del reporte HTML de resultados
# ---------------------------------------------------------------------------
class ReportGenerator:
    def __init__(self, results, project, params, methods, path_eval, init, filt, k_val=None, t0: float = 0.0):
        self.res = results; self.proj = project; self.params = params; self.meth = methods
        self.eval = path_eval; self.init = init; self.filt = filt; self.k_val = k_val
        self._t0 = t0  # tiempo de inicio — elapsed se calcula al escribir el HTML

    def _fmt_elapsed(self, secs: float) -> str:
        """Formatea segundos en texto legible."""
        s = int(round(secs))
        if s < 60:
            return f"{s} s"
        m, s = divmod(s, 60)
        if m < 60:
            return f"{m} min {s:02d} s"
        h, m = divmod(m, 60)
        return f"{h} h {m:02d} min {s:02d} s"

    def write(self, path, fb, open_it):
        """Genera el HTML con tiempo provisional. El tiempo final y la apertura
        del navegador ocurren en postProcessAlgorithm() vía rewrite_with_time()."""
        if not path: return
        if path == 'TEMPORARY_OUTPUT': path = QgsProcessingUtils.generateTempFilename('reporte.html')
        self._html_path = path
        self._html_fb   = fb
        self._html_open = open_it
        try:
            # Escribir con placeholder — se reemplazará en rewrite_with_time()
            self._elapsed_at_write = 0.0
            with open(path, 'w', encoding='utf-8') as f:
                f.write(self._html())
            # NO abrir aquí — se abre en rewrite_with_time() con tiempo correcto
        except Exception as e:
            fb.pushWarning(f"Error generando HTML: {e}")

    def rewrite_with_time(self, elapsed: float) -> None:
        """Reescribe el HTML con el tiempo real de ejecución y abre el navegador.
        Se llama desde postProcessAlgorithm() — después de cargar todas las capas.
        """
        self._elapsed_at_write = elapsed
        try:
            with open(self._html_path, 'w', encoding='utf-8') as f:
                f.write(self._html())
            if self._html_open:
                webbrowser.open(Path(self._html_path).as_uri())
        except Exception as e:
            if self._html_fb:
                self._html_fb.pushWarning(f"Error actualizando HTML con tiempo: {e}")

    def _html(self):
        rows = ""
        for r in self.res:
            nni = r.get('indice_nn'); s = self.params.sample_size; o = r.get('n_puntos', 0)
            rej = ((s-o)/s)*100 if s>0 else 0
            pat, cal = ("N/A", "N/A") if nni is None else (
                ("Agrupado Fuerte", "No Aceptable") if nni<0.8 else 
                ("Agrupado Moderado", "No Aceptable") if nni<0.95 else 
                ("Aleatorio", "Óptimo") if nni<=1.05 else 
                ("Disperso Moderado", "Bueno") if nni<=1.2 else 
                ("Disperso Fuerte", "Aceptable") if nni<=1.5 else ("Extremadamente Disperso", "Aceptable")
            )
            nni_txt = f"{nni:.4f}".replace('.', ',') if nni is not None else "NULL"
            rej_txt = f"{rej:.2f}".replace('.', ',')
            # Se corrigió el acceso a 'polygons_hit_count'
            polys = r.get('polygons_hit_count', 'N/A')
            contig = r.get('contiguous_count', 'N/A')
            rows += f"<tr><td>{r.get('muestra')}</td><td>{s}</td><td>{o}</td><td>{rej_txt}%</td><td>{polys}</td><td>{nni_txt}</td><td>{pat}</td><td>{cal}</td><td>{contig}</td></tr>"

        # Datos para JS
        max_cont = max([r.get('contiguous_count', 0) or 0 for r in self.res]) if self.res else 0
        max_poly = max([r.get('polygons_hit_count', 0) or 0 for r in self.res]) if self.res else 0
        sugg_poly = math.ceil(max_poly * 1.1) if max_poly > 0 else 10
        sugg_cont = math.ceil(max_cont * 1.1) if max_cont > 0 else 5
        
        js_labels = [r.get('muestra','') for r in self.res]
        js_nni = [r.get('indice_nn') or 'null' for r in self.res]
        js_eff = [r.get('n_puntos',0) for r in self.res]
        js_poly = [r.get('polygons_hit_count') or 'null' for r in self.res]
        js_cont = [r.get('contiguous_count') or 'null' for r in self.res]
        js_rej = [((self.params.sample_size - r.get('n_puntos', 0)) / self.params.sample_size) * 100 for r in self.res]

        def fc(v): return f"{v:.2f}".replace('.', ',') if v is not None else "0,00"
        def fi(v): return f"{v:,}".replace(',', ' ') if v is not None else "0"

        return f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Reporte — Muestreo Espacial de Puntos</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/chartjs-plugin-annotation/2.2.1/chartjs-plugin-annotation.min.js"></script>
    <style>
        body {{ font-family: sans-serif; margin: 2em; background-color: #f4f4f9; color: #333; line-height: 1.6;}}
        .container {{ max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 0 15px rgba(0,0,0,0.1); }}
        h1, h2, h3 {{ color: #0056b3; border-bottom: 2px solid #ddd; padding-bottom: 10px; margin-top: 1.5em;}}
        h3 {{ border-bottom: 1px solid #eee; color: #444;}}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 2em; box-shadow: 0 2px 3px rgba(0,0,0,0.1);}}
        th, td {{ text-align: left; padding: 12px; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #0056b3; color: white; font-weight: bold;}}
        tr:nth-child(even) {{background-color: #f9f9f9;}}
        tr:hover {{background-color: #f1f7ff;}}
        .summary, .best-samples {{ margin-top: 2em; padding: 25px; border-left: 5px solid; border-radius: 5px; }}
        .summary {{ background-color: #e8f5e9; border-color: #4CAF50; }}
        .best-samples {{ background-color: #e3f2fd; border-color: #2196F3; }}
        .config-table td {{ padding: 8px; border: 1px solid #ddd; }}
        .chart-container {{ position: relative; margin: auto; height: 50vh; width: 80vw; margin-top: 2em; margin-bottom: 4em; }}
        .footer {{ margin-top: 3em; padding-top: 1em; border-top: 1px solid #ddd; font-size: 0.9em; color: #777; text-align: center; }}
        .interpretation-table {{ width: 100%; margin: 2em auto; }}
        .interpretation-table th {{ background-color: #4CAF50; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Análisis de Resultados — Muestreo Espacial de Puntos</h1>
        <p><strong>Proyecto:</strong> {self.proj}</p>
        <p><strong>Fecha de generación:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} &nbsp;|&nbsp; <strong>Ejecución completada:</strong> {self._fmt_elapsed(getattr(self, "_elapsed_at_write", 0.0))}</p>
        
        <h3>Parámetros de Configuración Utilizados</h3>
        <table class="config-table">
            <tbody><tr><td><strong>Método de muestreo:</strong></td><td>{self.meth.get(self.params.method)}</td>
                <td><strong>Puntos marco muestral (inicial):</strong></td><td>{fi(self.init)}</td></tr>
            <tr><td><strong>Tamaño de muestra solicitado (n):</strong></td><td>{self.params.sample_size}</td>
                <td><strong>Puntos marco muestral (filtrado):</strong></td><td>{fi(self.filt)}</td></tr>
            <tr><td><strong>Iteraciones:</strong></td><td>{self.params.num_iterations}</td>
                <td><strong>Distancia mínima:</strong></td><td>{self.params.min_sample_distance}</td></tr>
            <tr><td><strong>Orden Hilbert:</strong></td><td>{self.params.hilbert_order}</td>
                <td><strong>Distancia al perímetro:</strong></td><td>{self.params.perimeter_distance}</td></tr>
            <tr><td><strong>Grupos (k):</strong></td><td>{self.params.num_groups if self.params.num_groups > 0 else f'Auto ({self.res[0].get("k_efectivo", "N/A") if self.res else "N/A"})'}</td>
                <td><strong>Corrección de Retracción:</strong></td><td>{'Sí' if self.params.apply_pullback else 'No aplica'}</td></tr>
            <tr><td><strong></strong></td><td></td>
                <td><strong>Radio parcela (Retracción):</strong></td><td>{self.params.plot_radius if self.params.apply_pullback else 'No aplica'}</td></tr>
        </tbody></table>
        
        <h2>Marco de Validación Cuantitativo de la Curva de Hilbert</h2>
        <div class="summary" style="background-color:#e3f2fd; border-color:#2196F3;">
            <h4>Escala Espacial de Referencia y Definiciones</h4>
            <ul>
                <li><b>Distancia media de malla (h):</b> Es la separación teórica esperada entre puntos si estuvieran distribuidos de manera equiespaciada sobre el área efectiva. Se calcula como <b>h = Raíz (A/N)</b>.</li>
                <li><b>Distancia media de la ruta (d):</b> Es el promedio de las distancias entre puntos consecutivos en la ruta ordenada por Hilbert.</li>
            </ul>
            <p><b>Valores de Referencia:</b> Área Efectiva (A): {fc(self.eval.get('area_km2'))} km² | Puntos (N): {fi(self.eval.get('num_points'))} | Distancia media de malla (h): {fc(self.eval.get('mesh_distance_h'))} m</p>
        </div>
        <table class="config-table">
            <tbody><tr><th colspan="4" style="text-align:center; background-color:#4CAF50;">Calidad y Cobertura del Ordenamiento</th></tr>
            <tr>
                <td><strong>Orden de Hilbert (n):</strong></td><td>{self.eval.get('order')}</td>
                <td><strong>Colisiones Detectadas:</strong></td><td>{self.eval.get('collisions')} ({fc(self.eval.get('collision_ratio'))}%)</td>
            </tr>
            <tr>
                <td><strong>Espaciado de Malla (g):</strong></td><td>{fc(self.eval.get('grid_spacing_g'))} m {'(Digitado)' if self.params.manual_grid_spacing > 0 else '(Estimado)'}</td>
                <td><strong>Índice de Cobertura (h/g):</strong></td><td>{fc(self.eval.get('coverage_index')) if self.params.manual_grid_spacing > 0 else "No Aplica (no se digitó 'g')"}</td>
            </tr>
            <tr>
                <td><strong>Puntos Duplicados Detectados:</strong></td><td colspan="3">{self.eval.get('duplicate_count')}</td>
            </tr>
            <tr>
                <td><strong>Índice de Compacidad (Forma):</strong></td><td>{self.eval.get('indice_compacidad',0):.3f}</td>
                <td colspan="2"><strong>Patrón de Forma:</strong> {self._get_compacidad_desc(self.eval.get('indice_compacidad',0))}</td>
            </tr>
        </tbody></table>
        
        <table class="config-table">
            <tbody><tr><th colspan="4" style="text-align:center; background-color:#4CAF50;">Métricas Normalizadas y Diagnóstico</th></tr>
            <tr><th>Métrica</th><th>Rango (Excelente)</th><th>Valor Calculado</th><th>Diagnóstico</th></tr>
            <tr style="background-color:{self._get_color(self.eval.get('efficiency_index_R',0), 0.8, 1.5)};"><td>Índice de eficiencia (R = d/h)</td><td>0.8 - 1.5</td><td>{fc(self.eval.get('efficiency_index_R'))}</td><td>{self._get_diag(self.eval.get('efficiency_index_R',0), 0.8, 1.5)}</td></tr>
            <tr style="background-color:{self._get_color(self.eval.get('std_dev_norm',0), 0.5, 1.0)};"><td>Desviación estándar normalizada (sigma/h)</td><td>0.5 - 1.0</td><td>{fc(self.eval.get('std_dev_norm'))}</td><td>{self._get_diag(self.eval.get('std_dev_norm',0), 0.5, 1.0)}</td></tr>
            <tr style="background-color:{self._get_color(self.eval.get('length_norm',0), 0.9, 1.6)};"><td>Longitud total normalizada (L/Nh)</td><td>0.9 - 1.6</td><td>{fc(self.eval.get('length_norm'))}</td><td>{self._get_diag(self.eval.get('length_norm',0), 0.9, 1.6)}</td></tr>
        </tbody></table>
        
        <table class="config-table">
            <tbody><tr><th colspan="4" style="text-align:center; background-color:#4CAF50;">Índice y Descripción del Patrón Espacial</th></tr>
            <tr><th>Índice</th><th>Valor Calculado</th><th>Descripción del Patrón</th><th>Valoración</th></tr>
            <tr>
                <td>Coeficiente de Variación (CV = sigma/d)</td>
                <td>{fc(self.eval.get('cv'))}</td>
                <td>{self._get_cv_desc(self.eval.get('cv',0))}</td>
                <td>{self._get_diag(self.eval.get('cv',0), 0.5, 0.8)}</td>
            </tr>
        </tbody></table>
        
        <h3>Guía de Interpretación del Diagnóstico de Ruta</h3>
        <div class="summary" style="border-color: #ff9800; background-color: #fff3e0;">
            <p><b>Nota Aclaratoria General:</b> Las escalas de valoración presentadas en esta guía (ej. "Excelente", "Adecuado") son referencias diseñadas para estandarizar el análisis. Es fundamental entender que no son umbrales absolutos. Por ejemplo, un patrón "Disperso" puede ser el objetivo deseado en un muestreo sistemático Hilbert, mientras que un patrón "Agrupado" puede ser el reflejo correcto de un paisaje fragmentado. La interpretación de cualquier métrica debe considerar siempre los <b>objetivos de su muestreo</b> y las <b>características particulares de su área de estudio</b>.</p>
        </div>
        
        <div class="summary" style="border-color: #ffc107; background-color: #fff8e1;">
            <h4>Índice de Compacidad Poligonal</h4>
            <p>Este índice describe la forma del polígono de estudio (1.0 = círculo). Valores bajos indican formas complejas que pueden aumentar los efectos de borde y, en consecuencia, el valor del Coeficiente de Variación (CV).</p>
            <table class="interpretation-table">
                <thead><tr><th>Rango del Índice</th><th>Descripción de la Forma</th></tr></thead>
                <tbody>
                    <tr><td>&gt; 0.85</td><td>Compacto</td></tr>
                    <tr><td>0.6 - 0.85</td><td>Moderadamente Compacto</td></tr>
                    <tr><td>0.4 - 0.6</td><td>Alargado o Irregular</td></tr>
                    <tr><td>&lt; 0.4</td><td>Muy Alargado o Fragmentado</td></tr>
                </tbody>
            </table>
            <h4>Diagnóstico de Métricas Normalizadas (Excelente, Adecuado, Revisión, Deficiente)</h4>
            <p>Un diagnóstico "Deficiente" o "Revisión" es una alerta sobre una desviación del comportamiento ideal. Se recomienda verificar lo siguiente:</p>
            <ol>
                <li><strong>Orden de la Curva de Hilbert:</strong> Verifique si el orden utilizado es adecuado.</li>
                <li><strong>Puntos Duplicados:</strong> Verifique si existen puntos con coordenadas idénticas.</li>
            </ol>
            <h4>Escala de Validación de Rutas Hilbert</h4>
            <table class="interpretation-table">
                <thead><tr style="background-color: #333; color: white;"><th>Métrica</th><th>Excelente (Verde)</th><th>Adecuado (Amarillo)</th><th>Revisión (Naranja)</th><th>Deficiente (Rojo)</th></tr></thead>
                <tbody>
                    <tr><td><b>Índice de Cobertura (h/g)</b></td><td>&gt; 1.5</td><td>0.7 - 1.5</td><td colspan="2">&lt; 0.7 (Baja Resolución)</td></tr>
                    <tr><td><b>Índice de eficiencia (R)</b></td><td>0.8 - 1.5</td><td>1.5 - 2.0</td><td>2.0 - 2.5</td><td>&gt; 2.5</td></tr>
                    <tr><td><b>Desviación estándar normalizada (sigma/h)</b></td><td>0.5 - 1.0</td><td>1.0 - 1.5</td><td>1.5 - 2.0</td><td>&gt; 2.0</td></tr>
                    <tr><td><b>Longitud total normalizada (L/Nh)</b></td><td>0.9 - 1.6</td><td>1.6 - 2.2</td><td>2.2 - 2.8</td><td>&gt; 2.8</td></tr>
                    <tr><td><b>Coef. de variación (CV)</b></td><td>&lt; 0.50</td><td>0.50 - 0.80</td><td>0.80 - 1.20</td><td>>= 1.20</td></tr>
                </tbody>
            </table>
        </div>
        
        <h3>Guía de Interpretación del Análisis de Vecino Más Cercano (AVMC)</h3>
        <p>El <strong>Índice de Vecino Más Cercano (IVMC)</strong> es una medida que compara la distancia media observada entre cada punto y su vecino más cercano, con la distancia media esperada para una distribución hipotética aleatoria.</p>
        <table class="interpretation-table">
            <caption>Clasificación del Patrón Espacial según el IVMC</caption>
            <thead><tr><th>Rango IVMC (Aprox.)</th><th>Patrón</th><th>Calidad</th></tr></thead>
            <tbody>
                <tr><td>IVMC &lt; 0.8</td><td>Agrupado Fuerte</td><td>No Aceptable</td></tr>
                <tr><td>0.8 <= IVMC &lt; 0.95</td><td>Agrupado Moderado</td><td>No Aceptable</td></tr>
                <tr><td>0.95 <= IVMC <= 1.05</td><td>Aleatorio</td><td>Óptimo</td></tr>
                <tr><td>1.05 &lt; IVMC <= 1.2</td><td>Disperso Moderado</td><td>Bueno</td></tr>
                <tr><td>1.2 &lt; IVMC <= 1.5</td><td>Disperso Fuerte</td><td>Aceptable</td></tr>
                <tr><td>IVMC &gt; 1.5</td><td>Extremadamente Disperso</td><td>Aceptable</td></tr>
            </tbody>
        </table>

        <h2>Resultados Detallados por Muestra Generada</h2>
        <table>
            <tbody><tr><th>Muestra</th><th>Puntos Solicitados</th><th>Puntos Obtenidos</th><th>% Rechazo</th><th>Polígonos Muestreados</th><th>IVMC (Corregido)</th><th>Patrón</th><th>Calidad</th><th>Número de Puntos Contiguos</th></tr>
            {rows}
        </tbody></table>

        <h2>Visualización Comparativa de Resultados</h2>
        <div class="chart-container"><canvas id="nniChart"></canvas></div>
        <div class="chart-container"><canvas id="rejectionChart"></canvas></div>
        <div class="chart-container"><canvas id="effectiveChart"></canvas></div>
        <div class="chart-container"><canvas id="polygonsChart"></canvas></div>
        <div class="chart-container"><canvas id="contiguousChart"></canvas></div>

        <div class="best-samples">
            <h2>Mejores Muestras Identificadas</h2>
            {self._get_best()}
        </div>

        <div class="footer">
            <p>Reporte generado: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}</p>
            <p>Complemento QGIS: <b>Muestreo Espacial de Puntos por Curva de Hilbert</b> &nbsp;|&nbsp; Versión {MuestreoEspacialPuntos.VERSION}</p>
            <p>Desarrollado por Jorge Fallas (jfallas56@gmail.com)</p>
        </div>
    </div>
    <script>
        const labels = {js_labels};
        new Chart(document.getElementById('nniChart'), {{ type: 'bar', data: {{ labels: labels, datasets: [{{ label: 'IVMC', data: {js_nni}, backgroundColor: 'rgba(54, 162, 235, 0.6)', spanGaps: true }}] }}, options: {{ scales: {{ y: {{ beginAtZero: false }} }} }} }});
        new Chart(document.getElementById('rejectionChart'), {{ type: 'bar', data: {{ labels: labels, datasets: [{{ label: '% Rechazo', data: {js_rej}, backgroundColor: 'rgba(255, 99, 132, 0.6)' }}] }}, options: {{ scales: {{ y: {{ min: 0, max: 50 }} }} }} }});
        new Chart(document.getElementById('effectiveChart'), {{ type: 'bar', data: {{ labels: labels, datasets: [{{ label: 'Puntos', data: {js_eff}, backgroundColor: 'rgba(75, 192, 192, 0.6)' }}] }}, options: {{ scales: {{ y: {{ min: 0 }} }} }} }});
        
        // Fixed: Ensure max suggested value is valid for Polygons
        const maxPoly = {max_poly};
        const suggMaxPoly = maxPoly > 0 ? Math.ceil(maxPoly * 1.1) : 10;
        
        new Chart(document.getElementById('polygonsChart'), {{ 
            type: 'bar', 
            data: {{ 
                labels: labels, 
                datasets: [{{ 
                    label: 'Polígonos', 
                    data: {js_poly}, 
                    backgroundColor: 'rgba(255, 159, 64, 0.6)',
                    spanGaps: true
                }}] 
            }},
            options: {{
                scales: {{
                    y: {{
                        min: 0,
                        suggestedMax: suggMaxPoly
                    }}
                }}
            }}
        }});

        // Fixed: Ensure max suggested value is valid for Contiguous
        const maxCont = {max_cont};
        const suggMaxCont = maxCont > 0 ? Math.ceil(maxCont * 1.1) : 5;

        new Chart(document.getElementById('contiguousChart'), {{ 
            type: 'line', 
            data: {{ 
                labels: labels, 
                datasets: [{{ 
                    label: 'Contiguos', 
                    data: {js_cont}, 
                    borderColor: 'rgb(153, 102, 255)', 
                    fill: true,
                    spanGaps: true
                }}] 
            }},
            options: {{
                scales: {{
                    y: {{
                        beginAtZero: true,
                        suggestedMax: suggMaxCont
                    }}
                }}
            }}
        }});
    </script>
</body></html>"""

    def _get_compacidad_desc(self, v):
        if v>0.85: return "Compacto"
        if v>0.6: return "Moderadamente Compacto"
        if v>0.4: return "Alargado o Irregular"
        return "Muy Alargado o Fragmentado"

    def _get_diag(self, v, min_v, max_v):
        """B-01: 4 categorias alineadas con la tabla de interpretacion del HTML."""
        if min_v <= v <= max_v: return "Excelente"
        rng = max_v - min_v
        ext_lo = min_v - rng * 0.5
        ext_hi = max_v + rng * 0.5
        if ext_lo <= v <= ext_hi: return "Adecuado"
        ext2_lo = min_v - rng
        ext2_hi = max_v + rng
        if ext2_lo <= v <= ext2_hi: return "Revision"
        return "Deficiente"
        
    def _get_color(self, v, min_v, max_v):
        """B-01: 4 colores coherentes con _get_diag."""
        if min_v <= v <= max_v: return "#e8f5e9"   # verde - Excelente
        rng = max_v - min_v
        if (min_v - rng*0.5) <= v <= (max_v + rng*0.5): return "#fff9c4"  # amarillo - Adecuado
        if (min_v - rng)     <= v <= (max_v + rng):     return "#fff3e0"  # naranja - Revision
        return "#ffebee"  # rojo - Deficiente

    def _get_cv_desc(self, v):
        if v<0.5: return "Uniforme / Regular"
        if v<0.8: return "Relativamente Homogéneo"
        if v<1.2: return "Moderadamente Agrupado"
        return "Agrupado / Fragmentado"

    def _get_best(self):
        if not self.res: return "<p>No hay datos.</p>"
        valid = [r for r in self.res if r.get('indice_nn') is not None]
        if not valid: return "<p>Sin IVMC válido.</p>"
        
        is_dispersion = self.params.min_sample_distance > 0 or self.params.method == SamplingMethod.SYSTEMATIC_HILBERT
        if is_dispersion:
            # Fix: Handle potential None values in sort key
            top3 = sorted(valid, key=lambda x: (x.get('indice_nn') or 0, x.get('polygons_hit_count') or 0), reverse=True)[:3]
            txt = "Objetivo: Dispersión (IVMC más alto)."
        else:
            optimal = [r for r in valid if 0.95 <= (r.get('indice_nn') or 0) <= 1.05]
            if not optimal:
                top3 = sorted(valid, key=lambda x: abs((x.get('indice_nn') or 99) - 1.0))[:3]
                txt = "Objetivo: Aleatoriedad (IVMC ~ 1.0)."
            else:
                top3 = sorted(optimal, key=lambda x: (x.get('polygons_hit_count') or 0), reverse=True)[:3]
                txt = "Objetivo: Aleatoriedad (Óptimos)."
                
        items = ""
        for s in top3:
            # Fix: Handle potential None in f-string
            ivmc = s.get('indice_nn')
            ivmc_str = f"{ivmc:.4f}".replace('.', ',') if ivmc is not None else "NULL"
            items += f"<li><strong>{s.get('muestra')}</strong>: IVMC {ivmc_str}, Polígonos {s.get('polygons_hit_count')}</li>"
        return f"<p>{txt}</p><ol>{items}</ol>"