# -*- coding: utf-8 -*-
"""
Complemento QGIS: Muestreo Espacial de Puntos
==================================================================
Autor:    Jorge Fallas (jfallas56@gmail.com)
Versión:  v1.0.0
Fecha:    2026-06-05
QGIS:     3.28 LTR+ / 3.44 LTR / 4.0 (Qt5 y Qt6) compatible
Python:   3.10+
Deps:     scikit-learn (opcional, método K-Medias)
          shapely (opcional, motor acelerado — incluido en QGIS 3.44+ y 4.0)
Licencia: GNU GPL v2

Propósito:
    Muestreo espacial de puntos sobre áreas poligonales con siete
    métodos: SH (Sistemático Hilbert), GH (Grupos Hilbert),
    GFC (Grupos Fila-Columna), AL (Aleatorio Simple), KM (K-Medias),
    Estr_Pts-AL (Estratificado por Puntos Aleatorio) y
    Estr_Pol-AL (Estratificado por Polígono Aleatorio).

Características principales (versión inicial v1.0.0):
    - checkParameterValues() con validación de SRC geográfico, rangos
      del orden Hilbert (1-11), iteraciones, dependencias K-Medias,
      formato JSON estratificado, coherencia NUM_GROUPS vs SAMPLE_SIZE,
      PLOT_RADIUS > 0 cuando APPLY_PULLBACK=True, y MANUAL_GRID_SPACING.
    - Aliases Qt5/Qt6 (QVariant → QMetaType) para compatibilidad con
      QGIS 3.28 LTR, 3.44 LTR y 4.0.
    - Motor Shapely v2 (contains_xy ufunc) / v1 (vectorized.contains) /
      GEOS PreparedGeometry como fallback automático.
    - makeValid() post-unaryUnion() en todas las ramas (estándar,
      estratos individuales y unión global).
    - releaseCache() para todas las geometrías preparadas, incluyendo
      las de estratos en STRATIFIED_BY_POLYGON.
    - _LayerVisibilityHandler con Qt.QueuedConnection para gestión
      cross-thread de visibilidad de capas.
    - VERSION como constante de clase.
    - group()/groupId() para categorización en Processing.
    - Umbrales NNI/IVMC como constantes de clase
      (NNI_RANDOM_LO/HI, NNI_DISPERSED_OK, SIGMA_H_UNIFORM).
    - Mensajes de error con valor recibido, rango esperado y acción
      sugerida (ej. SRC recomendado: CRTM05 EPSG:8908).
    - Logging tipificado en operaciones críticas (escritura del reporte
      HTML, extracción de contornos, limpieza de archivos temporales).
"""

import math
import random
import webbrowser
import gc
import copy
import json
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
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

from qgis.PyQt.QtCore import QCoreApplication, Qt, QObject, pyqtSlot
from qgis.core import (
    QgsProcessing, QgsFeatureSink, QgsProcessingException,
    QgsProcessingAlgorithm, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFeatureSink, QgsProcessingParameterNumber,
    QgsProcessingParameterEnum, QgsProcessingParameterDistance,
    QgsProcessingParameterBoolean, QgsProcessingParameterString,
    QgsProcessingParameterField,
    QgsFeature, QgsGeometry, QgsPointXY, QgsFields, QgsField,
    QgsWkbTypes, QgsProcessingFeedback, QgsProcessingContext,
    QgsDistanceArea, QgsVectorFileWriter,
    QgsProcessingParameterFileDestination, QgsVectorLayer,
    QgsCoordinateTransform, QgsProject, QgsFeatureRequest,
    QgsProcessingUtils, QgsSpatialIndex,
    QgsLayerTreeLayer
)

# --- Aliases de compatibilidad Qt5 (QVariant) → Qt6 (QMetaType) ---
try:
    from qgis.PyQt.QtCore import QMetaType
    _INT_TYPE = QMetaType.Type.Int
    _STR_TYPE = QMetaType.Type.QString
    _DBL_TYPE = QMetaType.Type.Double
except (ImportError, AttributeError):
    from qgis.PyQt.QtCore import QVariant
    _INT_TYPE = QVariant.Int
    _STR_TYPE = QVariant.String
    _DBL_TYPE = QVariant.Double

try:
    _TypeVectorPolygon = QgsProcessing.TypeVectorPolygon
    _TypeVectorPoint = QgsProcessing.TypeVectorPoint
except AttributeError:
    _TypeVectorPolygon = QgsProcessing.TypeVector
    _TypeVectorPoint = QgsProcessing.TypeVector

# --- Detección de motor Shapely para punto-en-polígono acelerado ---
try:
    import shapely
    import numpy as _np
    _SHAPELY_VERSION = shapely.__version__
    _shapely_major = int(_SHAPELY_VERSION.split('.')[0])
    if _shapely_major >= 2:
        import shapely as _shmod
        _shapely_contains_xy = _shmod.contains_xy
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
    GROUPS_ROW_COL = 4
    STRATIFIED_BY_FIELD = 5
    STRATIFIED_BY_POLYGON = 6


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
    stratum_field: str = ""
    polygon_stratum_field: str = ""
    strata_json: Dict[str, int] = None


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
    stratum_value: str = ""
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
        nuevo = PointData(
            id=self.id, x=self.x, y=self.y,
            hilbert_index=self.hilbert_index,
            original_index=self.original_index,
            was_corrected=self.was_corrected,
            correction_distance=self.correction_distance,
            rejected_in_iteration=self.rejected_in_iteration,
            stratum_value=self.stratum_value
        )
        return nuevo


@dataclass
class SamplePoint:
    feature: QgsFeature
    original_geometry: QgsGeometry
    was_corrected: bool = False
    correction_distance: float = 0.0
    rejected_in_iteration: int = 0
    stratum_value: str = ""

    def __post_init__(self):
        if self.original_geometry is None:
            self.original_geometry = QgsGeometry(self.feature.geometry())


@dataclass
class KMeansValidationData:
    iteration: int
    cluster_stats: List[Dict]
    proportionality: Optional[Dict] = None


@dataclass
class HilbertGroupValidationData:
    iteration: int
    proportionality: Dict

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
    def get_polygon_boundary_as_line(poly: QgsGeometry, feedback: QgsProcessingFeedback) -> Optional[QgsGeometry]:
        if not poly or poly.isNull():
            return None
        if hasattr(poly, 'boundary'):
            try:
                boundary = poly.boundary()
                if boundary and not boundary.isNull() and not boundary.isEmpty():
                    return boundary
            except Exception as e:
                if feedback:
                    feedback.pushDebugInfo(f"[GeometryUtils.boundary] fallback: {e}")
        try:
            result = poly.convertToType(QgsWkbTypes.LineGeometry, True)
            if result and not result.isNull() and not result.isEmpty():
                return result
        except Exception as e:
            if feedback:
                feedback.pushDebugInfo(f"[GeometryUtils.convertToType] fallback: {e}")
        try:
            poly_geom = poly.asPolygon()
            if poly_geom and poly_geom[0]:
                return QgsGeometry.fromPolylineXY(poly_geom[0])
            multi = poly.asMultiPolygon()
            if multi and multi[0] and multi[0][0]:
                return QgsGeometry.fromPolylineXY(multi[0][0])
        except Exception as e:
            if feedback:
                feedback.pushWarning(f"[GeometryUtils] Fallo extrayendo contorno: {e}.")
        return None

    @staticmethod
    def is_geometry_valid(geom: QgsGeometry) -> bool:
        if not geom or geom.isNull():
            return False
        try:
            return not geom.isEmpty() and geom.isGeosValid()
        except Exception:
            return False


class DistanceCalculator:
    def __init__(self, distance_area: QgsDistanceArea): self.da = distance_area

    def calculate_distance(self, point1: QgsPointXY, point2: QgsPointXY) -> float:
        if self.da.sourceCrs().isGeographic():
            return self.da.measureLine(point1, point2)
        return math.hypot(point2.x() - point1.x(), point2.y() - point1.y())

# --- Handler de visibilidad cross-thread ---


class _LayerVisibilityHandler(QObject):
    def __init__(self, framework_names: set, pending: list, expected: int = 1, report_gen=None, t0: float = 0.0):
        super().__init__()
        self._names = framework_names
        self._pending = pending
        self._expected = max(1, expected)
        self._count = 0
        self._done = False
        self._report_gen = report_gen
        self._t0 = t0
        from qgis.PyQt.QtCore import QCoreApplication
        self.moveToThread(QCoreApplication.instance().thread())

    @pyqtSlot('QgsMapLayer*')
    def on_layer_added(self, layer):
        if self._done or not layer:
            return
        if layer.name() in self._names:
            self._count += 1
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
        import time as _time
        root = QgsProject.instance().layerTreeRoot()
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

        if self._report_gen is not None and self._t0 > 0:
            try:
                elapsed = _time.time() - self._t0
                self._report_gen.rewrite_with_time(elapsed)
            except Exception:
                pass
            self._report_gen = None


class MuestreoEspacialPuntos(QgsProcessingAlgorithm):
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
    STRATUM_FIELD = 'STRATUM_FIELD'
    POLYGON_STRATUM_FIELD = 'POLYGON_STRATUM_FIELD'
    STRATA_JSON = 'STRATA_JSON'
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
    DEDUP_TOLERANCE_M: float = 0.10

    METHOD_NAMES = {
        SamplingMethod.SYSTEMATIC_HILBERT: "Sistemático Hilbert",
        SamplingMethod.SIMPLE_RANDOM: "Aleatorio Simple",
        SamplingMethod.STRATIFIED_HILBERT: "Grupos Hilbert (Grupos 1D)",
        SamplingMethod.KMEANS_GROUPS: "Grupos K-Medias (Conglomerados 2D)",
        SamplingMethod.GROUPS_ROW_COL: "Grupos Fila-Col (NO→SE)",
        SamplingMethod.STRATIFIED_BY_FIELD: "Estratificado por Puntos (JSON) [Estr_Pts-AL]",
        SamplingMethod.STRATIFIED_BY_POLYGON: "Estratificado por Polígono (JSON) [Estr_Pol-AL]"
    }

    METHOD_ABBREVIATIONS = {
        SamplingMethod.SYSTEMATIC_HILBERT: "SH",
        SamplingMethod.SIMPLE_RANDOM: "AL",
        SamplingMethod.STRATIFIED_HILBERT: "GH",
        SamplingMethod.KMEANS_GROUPS: "KM",
        SamplingMethod.GROUPS_ROW_COL: "GFC",
        SamplingMethod.STRATIFIED_BY_FIELD: "Estr_Pts-AL",
        SamplingMethod.STRATIFIED_BY_POLYGON: "Estr_Pol-AL"
    }

    UNIT_CONVERSIONS = {
        DistanceUnit.METERS: 1.0, DistanceUnit.KILOMETERS: 1000.0,
        DistanceUnit.FEET: 0.3048, DistanceUnit.MILES: 1609.34
    }

    FIELD_SAMPLE_ID = "muestra_id_proc"
    FIELD_HILBERT_IDX = "Hilbert_idx_proc"
    FIELD_NWSE_IDX = "NW_SE_idx_proc"
    FIELD_ITER_NUM = "iter_num_proc"
    FIELD_METHOD = "metodo_proc"
    FIELD_CORRECTED = "correccion_proc"
    FIELD_DIST_BORDER = "dist_borde_m"
    FIELD_REJECTED_ITER = "rechazado_iter"
    FIELD_NEAREST_NEIGHBOR = "dist_vecino_m"
    FIELD_SAMPLE_IVMC = "IVMC_muestra_proc"
    FIELD_NEIGHBOR_CHECK = "Vecino"
    FIELD_STRATUM = "estrato_proc"

    # Umbrales estadísticos para clasificación NNI/IVMC en el reporte HTML.
    # Cambiarlos aquí los actualiza en todos los puntos del código.
    NNI_RANDOM_LO: float = 0.95   # Borde inferior del rango "aleatorio puro"
    NNI_RANDOM_HI: float = 1.05   # Borde superior del rango "aleatorio puro"
    NNI_DISPERSED_OK: float = 1.20  # Umbral "disperso aceptable" (método SH)
    SIGMA_H_UNIFORM: float = 0.50   # Umbral σ/h "uniforme"

    def __init__(self):
        super().__init__()
        self.has_as_meters = hasattr(QgsProcessingParameterDistance, 'asMeters')
        self._sample_layers_pending: List[Tuple[str, str]] = []
        self._framework_layer_names: set = set()
        self._visibility_handler = None
        self._current_method: SamplingMethod = SamplingMethod.SYSTEMATIC_HILBERT
        self._has_rejected: bool = False
        self._report_gen = None
        self._hilbert_path_created = False
        self._html_path_pending = ''
        self._html_open_pending = True
        self._t0: float = 0.0

    def tr(self, text: str) -> str:
        return QCoreApplication.translate('Processing', text)

    def createInstance(self) -> 'MuestreoEspacialPuntos':
        return MuestreoEspacialPuntos()

    def _connect_layers_added_signal(self, framework_names: set, pending: list, expected: int = 1) -> None:
        handler = _LayerVisibilityHandler(framework_names, pending, expected, report_gen=self._report_gen, t0=self._t0)
        self._visibility_handler = handler
        self._report_gen = None
        QgsProject.instance().layerWasAdded.connect(handler.on_layer_added, Qt.QueuedConnection)

    def postProcessAlgorithm(self, context: QgsProcessingContext, feedback: QgsProcessingFeedback) -> Dict:
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

    def icon(self):
        from qgis.PyQt.QtGui import QIcon
        import os as _os
        _icon_path = _os.path.join(_os.path.dirname(__file__), "icon.png")
        if _os.path.isfile(_icon_path):
            return QIcon(_icon_path)
        return super().icon()

    def name(self) -> str:
        ver = self.VERSION.replace('.', '_').replace(' ', '_')
        return f'muestreo_espacial_puntos_{ver}'

    def displayName(self) -> str:
        return self.tr('Muestreo Espacial de Puntos')

    def group(self) -> str:
        return self.tr('Muestreo Espacial')

    def groupId(self) -> str:
        return 'muestreo_espacial'

    def shortHelpString(self) -> str:
        if HAS_SHAPELY:
            api_label = "contains_xy ufunc" if _shapely_use_v2 else "vectorized legacy"
            motor = f"Shapely {_SHAPELY_VERSION} ({api_label})"
        else:
            motor = "PreparedGeometry GEOS — estándar"
        return f"""<p><b>Descripción:</b> Genera muestras espaciales sobre una capa de puntos (marco muestral) delimitada por un polígono de área de estudio. Produce hasta <i>n</i> iteraciones independientes, cada una como capa vectorial nombrada con sus métricas de calidad (IVMC, índice de cobertura, CV). El reporte HTML compara las iteraciones y recomienda las tres mejores según el objetivo estadístico (dispersión o aleatoriedad). Compatible con QGIS 3.28 LTR, 3.44 LTR y 4.0 (Qt5/Qt6).</p>  
<p><b>Notación numérica:</b> decimal con coma (3,14) · miles con espacio fino (1 234,56).</p>
<p>&#9432; Para crear la malla de puntos use el complemento <b>Crear_Malla_Puntos</b>: <a href="https://github.com/jfallas56-CR/Crear-malla-de-Puntos">https://github.com/jfallas56-CR/Crear-malla-de-Puntos</a></p>  
<h3>Muestreo Espacial de Puntos</h3>
<p>&#9888; <b>REQUISITO:</b> SRC proyectado en metros (ej: CRTM05 EPSG:8908, UTM). Un SRC geográfico cancela la ejecución.</p>  
<h3>Descripción General</h3>
<p>La <b>Curva de Hilbert</b> es el núcleo de los métodos <b>SH</b> y <b>GH</b>: al ser una línea continua que recorre el área de estudio preservando la localidad espacial, garantiza que puntos cercanos en el orden de la curva también sean geográficamente próximos. Esto mejora la cobertura y representatividad espacial de SH y GH respecto a la selección completamente aleatoria.</p>  
<p>El <b>Aleatorio Simple (AL)</b> no busca cobertura espacial sino representatividad estadística: es el único método donde las fórmulas del MAS (varianza, intervalos de confianza) son directamente válidas sin corrección de diseño. Los métodos <b>GFC</b> (orden de filas NO→SE) y <b>KM</b> (proximidad geográfica 2D) usan otros criterios de ordenamiento espacial. Los métodos <b>Estr_Pts-AL</b> y <b>Estr_Pol-AL</b> seleccionan aleatoriamente dentro de estratos definidos vía JSON.</p>  
<p>En todos los métodos el ordenamiento Hilbert del marco completo se calcula y queda disponible como referencia de diagnóstico en el reporte HTML.</p>  
<h3>Métodos de Muestreo</h3>
<p><b>Sistemático Hilbert (SH):</b> Selección a paso fijo con inicio aleatorio a lo largo de la Curva de Hilbert. Al preservar la localidad espacial, garantiza cobertura dispersa sobre toda el área (IVMC &gt; 1,2). Produce mayor representatividad geográfica que AL.</p>  
<p><b>Aleatorio Simple (AL):</b> Selección al azar sin reemplazo, sin ordenamiento espacial previo. Único método donde las fórmulas del MAS (varianza, intervalos de confianza) son directamente válidas. Menor cobertura espacial que SH y GH, pero validez estadística directa.</p>  
<p><b>Grupos Hilbert (GH):</b> Divide el marco en k grupos 1D según el orden de la Curva de Hilbert —preservando localidad espacial— y selecciona aleatoriamente dentro de cada grupo. Combina cobertura espacial con selección aleatoria por estrato.</p>  
<p><b>Grupos Fila-Col (GFC):</b> Estratificación en k grupos por orden de filas NO→SE. Cada punto recibe un ID secuencial 1..N.</p>  
<p><b>Grupos K-Medias (KM):</b> Agrupa por proximidad geográfica 2D y selecciona dentro de cada grupo. Requiere <b>scikit-learn</b>.</p>  
<p><b>Estratificado por Puntos (JSON) [Estr_Pts-AL]:</b> Selecciona tamaños de muestra variables usando las categorías de un campo de la capa de <i>puntos</i> (Ej: Tipo de Cobertura). La selección al interior de cada estrato es <b>Aleatoria</b>.</p>  
<p><b>Estratificado por Polígono (JSON) [Estr_Pol-AL]:</b> Usa un atributo de la capa de <i>polígonos</i> para aislar áreas separadas espacialmente (Ej: Búferes individuales, ID Parcela). La selección al interior de cada polígono es <b>Aleatoria</b>.</p>  
<h3>Parámetros Clave</h3>
<p><b>Campo de Estrato (Puntos o Polígonos):</b> El atributo que divide los estratos según el método JSON elegido.</p>
<p><b>Tamaños en JSON:</b> Diccionario donde la clave es el estrato y el valor es el número de puntos a muestrear. Ej. <code>{{"Bufer 1": 100, "Bufer 2": 150}}</code>.</p>  
<p><b>Orden Hilbert (1-11):</b> Resolución de la cuadrícula (2ⁿ × 2ⁿ celdas). Defecto 10: adecuado para hasta 50 000+ puntos.</p>  
<p><b>Distancia mínima entre puntos:</b> Separación mínima garantizada con índice espacial O(n log n). Puntos excluidos van a la capa <i>Rechazados</i>.</p>  
<p><b>Corrección de borde (Retracción):</b> Mueve puntos a distancia ≥ radio de parcela del límite. Requiere <i>Radio de parcela</i>.</p>  
<h3>Métricas del Reporte HTML</h3>
<p><b>IVMC:</b> <code>IVMC = (d_obs/d_esp) × √(A_bbox/A_real)</code>, donde <code>d_esp = 0,5/√(n/A)</code>.</p>
<ul>
<li><b>Malla regular (compacta):</b> aleatoriedad pura no alcanzable. Valores esperados: rectangular ≈ 2,17; hexagonal ≈ 2,33.</li>  
<li><b>Malla regular fragmentada (multipolígono):</b> los valores IVMC son estructuralmente altos (típicamente 4–7) por el factor de corrección √(A_bbox/A_real). Use el ratio IVMC_método/IVMC_AL para comparar.</li>  
</ul>
<h3>Motor y Rendimiento</h3>
<p>Motor activo: <b>{motor}</b>.</p>
<div class="footer">Autor: Jorge Fallas (jfallas56@gmail.com) — {self.VERSION}</div>"""

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Capa de puntos de entrada'),
                [_TypeVectorPoint]))
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.SAMPLING_AREA,
                self.tr('Área de muestreo (Polígonos)'),
                [_TypeVectorPolygon]))
        self.addParameter(
            QgsProcessingParameterString(
                self.PROJECT_NAME,
                self.tr('Nombre del Proyecto'),
                defaultValue='Proyecto de Muestreo'))
        if not self.has_as_meters:
            self.addParameter(
                QgsProcessingParameterEnum(
                    self.DISTANCE_UNITS,
                    self.tr('Unidades'),
                    options=[
                        'Metros',
                        'Km',
                        'Pies',
                        'Millas'],
                    defaultValue=DistanceUnit.METERS))

        dist_suffix = '' if self.has_as_meters else self.tr(' - Ver unidad arriba')
        self.addParameter(
            QgsProcessingParameterDistance(
                self.PERIMETER_DISTANCE,
                self.tr('Distancia al perímetro')
                + dist_suffix,  # noqa: W503
                parentParameterName=self.INPUT,
                optional=True,
                defaultValue=0.0))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.SELECT_INSIDE_BUFFER,
                self.tr('Muestrear DENTRO del margen'),
                defaultValue=False))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.APPLY_PULLBACK,
                self.tr('Aplicar corrección de borde (Retracción)'),
                defaultValue=False))
        self.addParameter(
            QgsProcessingParameterDistance(
                self.PLOT_RADIUS,
                self.tr("Radio de parcela (para Retracción)")
                + dist_suffix,  # noqa: W503
                parentParameterName=self.INPUT,
                defaultValue=0.0,
                optional=True))
        self.addParameter(
            QgsProcessingParameterDistance(
                self.MIN_SAMPLE_DISTANCE,
                self.tr('Distancia mínima entre puntos')
                + dist_suffix,  # noqa: W503
                parentParameterName=self.INPUT,
                optional=True,
                defaultValue=0.0))
        self.addParameter(
            QgsProcessingParameterDistance(
                self.MANUAL_GRID_SPACING,
                self.tr('Distancia de malla conocida (g)')
                + dist_suffix,  # noqa: W503
                parentParameterName=self.INPUT,
                optional=True,
                defaultValue=0.0))

        method_opts = [
            self.METHOD_NAMES[m] for m in [
                SamplingMethod.SYSTEMATIC_HILBERT,
                SamplingMethod.SIMPLE_RANDOM,
                SamplingMethod.STRATIFIED_HILBERT,
                SamplingMethod.GROUPS_ROW_COL,
                SamplingMethod.STRATIFIED_BY_FIELD,
                SamplingMethod.STRATIFIED_BY_POLYGON
            ]
        ]
        if HAS_SKLEARN:
            method_opts.insert(4, self.METHOD_NAMES[SamplingMethod.KMEANS_GROUPS])

        self.addParameter(
            QgsProcessingParameterEnum(
                self.METHOD,
                self.tr('Método de muestreo'),
                options=method_opts,
                defaultValue=0))

        self.addParameter(
            QgsProcessingParameterField(
                self.STRATUM_FIELD,
                self.tr('Campo de Estrato en PUNTOS (Solo Estr_Pts-AL)'),
                type=QgsProcessingParameterField.Any,
                parentLayerParameterName=self.INPUT,
                optional=True))
        self.addParameter(
            QgsProcessingParameterField(
                self.POLYGON_STRATUM_FIELD,
                self.tr('Campo de Estrato en POLÍGONOS (Solo Estr_Pol-AL)'),
                type=QgsProcessingParameterField.Any,
                parentLayerParameterName=self.SAMPLING_AREA,
                optional=True))
        self.addParameter(
            QgsProcessingParameterString(
                self.STRATA_JSON,
                self.tr('Tamaños por estrato en JSON (Solo Estr_Pts-AL y Estr_Pol-AL)'),
                defaultValue='{\n  "Estrato 1": 100,\n  "Estrato 2": 150\n}',
                multiLine=True,
                optional=True))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.SAMPLE_SIZE,
                self.tr('Tamaño de muestra general (Se ignora en JSON)'),
                defaultValue=10))
        self.addParameter(QgsProcessingParameterNumber(self.NUM_ITERATIONS, self.tr('Iteraciones'), defaultValue=1))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.HILBERT_ORDER,
                self.tr('Orden Hilbert (1-11)'),
                defaultValue=10))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.NUM_GROUPS,
                self.tr('Número de grupos (k)'),
                defaultValue=0))
        self.addParameter(
            QgsProcessingParameterString(
                self.SAMPLE_WORD,
                self.tr('Palabra clave salida'),
                defaultValue='muestra',
                optional=True))
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_ALL,
                self.tr('Salida: Puntos filtrados'),
                optional=True))
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_HILBERT_PATH,
                self.tr('Salida: Ruta Hilbert'),
                optional=True))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CALC_DIST_ALL_POINTS,
                self.tr('Calcular distancia borde'),
                defaultValue=False))
        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT_HTML_REPORT,
                self.tr('Reporte HTML'),
                fileFilter='HTML files (*.html)',
                optional=True))
        self.addParameter(QgsProcessingParameterBoolean(self.OPEN_REPORT, self.tr('Abrir reporte'), defaultValue=True))
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_REJECTED,
                self.tr('Salida: Rechazados'),
                optional=True))

        _helps = {
            self.INPUT: ('<b>Marco muestral:</b> capa de puntos sobre la que se aplica el muestreo. Debe estar en un <b>SRC proyectado en metros</b>. Un SRC geográfico (grados) cancela la ejecución antes de iniciar.'),  
            self.SAMPLING_AREA: ('<b>Área de muestreo:</b> polígono o multipolígono que delimita la zona de interés. Geometrías inválidas se reparan con makeValid() antes de procesar.'),  
            self.PROJECT_NAME: ('<b>Nombre del proyecto:</b> etiqueta visible en el encabezado del reporte HTML.'),
            self.PERIMETER_DISTANCE: ('<b>Distancia al perímetro:</b> excluye una franja del área. Combinado con <i>Muestrear DENTRO del margen</i> define el área hábil.'),  
            self.SELECT_INSIDE_BUFFER: ('<b>Muestrear DENTRO del margen:</b> invierte la zona de selección. Solo tiene efecto si <i>Distancia al perímetro</i> > 0.'),  
            self.APPLY_PULLBACK: ('<b>Corrección de borde (Retracción):</b> mueve puntos a menos de <i>Radio de parcela</i> del límite al punto más cercano dentro del búfer interior.'),  
            self.PLOT_RADIUS: ('<b>Radio de parcela:</b> distancia mínima al límite para <i>Corrección de borde</i>.'),
            self.MIN_SAMPLE_DISTANCE: ('<b>Distancia mínima entre puntos:</b> separación mínima garantizada. Puntos rechazados van a <i>Salida: Rechazados</i>.'),  
            self.MANUAL_GRID_SPACING: ('<b>Distancia de malla conocida (g):</b> espaciado real de la malla de puntos para calcular el Índice de Cobertura (h/g).'),  
            self.METHOD: ('<b>Método de muestreo:</b><br>• SH, GH, AL, KM, GFC.<br>• Estr_Pts-AL: Estratificado Aleatorio guiado por campo de puntos.<br>• Estr_Pol-AL: Estratificado Aleatorio guiado por polígonos separados.'),  
            self.SAMPLE_SIZE: ('<b>Tamaño de muestra general:</b> puntos a seleccionar por iteración. Ignorado si se usan métodos JSON.'),  
            self.NUM_ITERATIONS: ('<b>Iteraciones:</b> número de muestras independientes a generar. Cada una produce una capa vectorial separada.'),  
            self.HILBERT_ORDER: ('<b>Orden Hilbert (1-11):</b> resolución de la cuadrícula para el ordenamiento espacial.'),  
            self.NUM_GROUPS: ('<b>Número de grupos (k):</b> para métodos GH, KM, GFC. k=0 es automático.'),
            self.SAMPLE_WORD: ('<b>Palabra clave salida:</b> prefijo en el nombre de cada capa de muestra.'),
            self.CALC_DIST_ALL_POINTS: ('<b>Calcular distancia borde:</b> agrega <i>dist_borde_m</i> a la capa <i>Salida: Puntos filtrados</i>.'),  
            self.OUTPUT_ALL: ('<b>Salida: Puntos filtrados:</b> todos los puntos del marco dentro del área. Incluye índices y atributos espaciales.'),  
            self.OUTPUT_HILBERT_PATH: ('<b>Salida: Ruta Hilbert:</b> línea que conecta los puntos del marco en el orden procesado (Hilbert o NO-SE).'),  
            self.OUTPUT_HTML_REPORT: ('<b>Reporte HTML:</b> documento avanzado con métricas de calidad y diagnóstico por iteración.'),  
            self.OPEN_REPORT: ('<b>Abrir reporte:</b> abre el HTML en el navegador al finalizar.'),
            self.OUTPUT_REJECTED: ('<b>Salida: Rechazados:</b> puntos excluidos por <i>Distancia mínima entre puntos</i>.'),  
            self.STRATUM_FIELD: ('<b>Campo de Estrato (Puntos):</b> Columna de la capa de puntos que define las categorías o clases para el muestreo estratificado. Solo tiene efecto si se elige el método Estratificado por Puntos [Estr_Pts-AL].'),  
            self.POLYGON_STRATUM_FIELD: ('<b>Campo de Estrato (Polígonos):</b> Columna de la capa de Área de Muestreo (ej. id_bufer o nombre_zona). Aísla cada polígono y asigna la muestra. Solo tiene efecto si se elige el método Estratificado por Polígono [Estr_Pol-AL].'),  
            self.STRATA_JSON: ('<b>Tamaños en JSON:</b> Formato de texto para asignar la muestra. Ej:<br><code>{"Bufer 1": 120, "Bufer 2": 50}</code><br>Usado por los métodos Estr_Pts-AL y Estr_Pol-AL. El Tamaño de Muestra general numérico se ignora; la cuota será la suma de los valores definidos aquí.'),  
        }
        for _pname, _htxt in _helps.items():
            _pd = self.parameterDefinition(_pname)
            if _pd:
                _pd.setHelp(_htxt)

    def checkParameterValues(self, parameters: Dict, context: QgsProcessingContext) -> Tuple[bool, str]:
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source and source.sourceCrs().isGeographic():
            return False, (
                f"[!] ERROR CRÍTICO — SRC geográfico detectado en la capa de puntos: "
                f"{source.sourceCrs().authid()} (unidades en grados). "
                f"Este algoritmo requiere un SRC proyectado en metros "
                f"(ej: CRTM05 EPSG:8908, UTM zona correspondiente). "
                f"Reproyecte la capa antes de continuar."
            )
        area = self.parameterAsSource(parameters, self.SAMPLING_AREA, context)
        if area and area.sourceCrs().isGeographic():
            return False, (
                f"[!] ERROR CRÍTICO — SRC geográfico en el área de muestreo: "
                f"{area.sourceCrs().authid()} (unidades en grados). "
                f"Reproyecte a un SRC en metros (ej: CRTM05 EPSG:8908, UTM)."
            )

        order = self.parameterAsInt(parameters, self.HILBERT_ORDER, context)
        if not (1 <= order <= 11):
            return False, f"[!] Orden Hilbert inválido: {order}. Debe estar entre 1 y 11."
        num_iter = self.parameterAsInt(parameters, self.NUM_ITERATIONS, context)
        if num_iter < 1:
            return False, f"[!] Número de iteraciones inválido: {num_iter}. Debe ser >= 1."

        _avail = [
            SamplingMethod.SYSTEMATIC_HILBERT,
            SamplingMethod.SIMPLE_RANDOM,
            SamplingMethod.STRATIFIED_HILBERT,
            SamplingMethod.GROUPS_ROW_COL,
            SamplingMethod.STRATIFIED_BY_FIELD,
            SamplingMethod.STRATIFIED_BY_POLYGON]
        if HAS_SKLEARN:
            _avail.insert(4, SamplingMethod.KMEANS_GROUPS)

        m_idx = self.parameterAsInt(parameters, self.METHOD, context)
        m_enum = _avail[m_idx] if m_idx < len(_avail) else SamplingMethod.SYSTEMATIC_HILBERT

        if m_enum == SamplingMethod.KMEANS_GROUPS and not HAS_SKLEARN:
            return False, (
                "[!] El método K-Medias requiere scikit-learn. "
                "Instálelo con: pip install scikit-learn (OSGeo4W Shell). "
                "Alternativamente, seleccione otro método de muestreo."
            )

        if m_enum in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
            strat_fld = self.parameterAsString(
                parameters,
                self.STRATUM_FIELD,
                context) if m_enum == SamplingMethod.STRATIFIED_BY_FIELD else self.parameterAsString(
                parameters,
                self.POLYGON_STRATUM_FIELD,
                context)
            strat_json = self.parameterAsString(parameters, self.STRATA_JSON, context)
            if not strat_fld:
                return False, f"[!] Método JSON requiere seleccionar el Campo de Estrato correspondiente ({
                    'Puntos' if m_enum == SamplingMethod.STRATIFIED_BY_FIELD else 'Polígonos'})."
            try:
                j_data = json.loads(strat_json)
                if not isinstance(j_data, dict):
                    return False, "[!] El JSON debe ser un diccionario (ej. {\"Bufer 1\": 100, \"Bufer 2\": 50})."
                for k, v in j_data.items():
                    if not isinstance(v, int) or v < 0:
                        return False, f"[!] El tamaño para el estrato '{k}' debe ser un entero >= 0 (recibido: {v!r})."
            except Exception as e:
                return False, f"[!] JSON inválido en 'Tamaños por estrato': {e}"
        else:
            sample_size = self.parameterAsInt(parameters, self.SAMPLE_SIZE, context)
            if sample_size < 1:
                return False, f"[!] Tamaño de muestra inválido: {sample_size}. Debe ser >= 1."

            # Coherencia número de grupos vs. tamaño de muestra (métodos agrupados).
            if m_enum in (SamplingMethod.STRATIFIED_HILBERT,
                          SamplingMethod.KMEANS_GROUPS,
                          SamplingMethod.GROUPS_ROW_COL):
                num_groups = self.parameterAsInt(parameters, self.NUM_GROUPS, context)
                if num_groups < 0:
                    return False, (
                        f"[!] Número de grupos inválido: {num_groups}. "
                        f"Debe ser >= 0 (use 0 para cálculo automático: k = mín(⌈√n⌉, ⌊√N⌋))."
                    )
                if num_groups > sample_size and num_groups > 0:
                    return False, (
                        f"[!] Número de grupos ({num_groups}) mayor que tamaño de muestra ({sample_size}). "
                        f"Reduzca el número de grupos, aumente el tamaño de muestra, "
                        f"o use k=0 para cálculo automático."
                    )

        # Corrección de borde requiere radio de parcela > 0.
        apply_pullback = self.parameterAsBoolean(parameters, self.APPLY_PULLBACK, context)
        if apply_pullback:
            try:
                plot_radius = self.parameterAsDouble(parameters, self.PLOT_RADIUS, context)
                if plot_radius <= 0:
                    return False, (
                        f"[!] Corrección de borde activada (apply_pullback=True) "
                        f"pero 'Radio de parcela' = {plot_radius}. "
                        f"Asigne un radio > 0 m o desactive la corrección de borde."
                    )
            except Exception:
                pass  # Parámetro no presente o tipo incompatible; ignorar.

        # Distancia de malla manual no puede ser negativa.
        try:
            grid_g = self.parameterAsDouble(parameters, self.MANUAL_GRID_SPACING, context)
            if grid_g < 0:
                return False, (
                    f"[!] Distancia de malla (g) inválida: {grid_g}. "
                    f"Debe ser >= 0 (use 0 para detección automática a partir del marco muestral)."
                )
        except Exception:
            pass

        return super().checkParameterValues(parameters, context)

    def processAlgorithm(
            self,
            parameters: Dict,
            context: QgsProcessingContext,
            feedback: QgsProcessingFeedback) -> Dict:
        import time as _time
        self._t0 = _time.time()
        _t0 = self._t0
        feedback.setProgress(0)
        feedback.pushInfo(f"--- Iniciando {self.displayName()} ---")
        if HAS_SHAPELY:
            feedback.pushInfo(f"[Motor] Shapely {_SHAPELY_VERSION} disponible.")
        temp_files = []
        try:
            layers, params, paths = self._setup(parameters, context)
            feedback.setProgress(5)

            data = self._prepare_data(layers, params, context, feedback)
            if data is None:
                feedback.reportError("AVISO: No se generaron muestras — marco muestral vacío o cruce inválido.")
                return {}

            feedback.setProgress(45)
            path_eval = self._evaluate_path(
                data['sorted_features'],
                data['area_for_filter'],
                params,
                data['distance_calculator'],
                data['duplicate_count'],
                feedback)
            feedback.setProgress(50)

            iter_res = self._run_iterations(data, params, parameters, context, feedback, parameters[self.SAMPLING_AREA])
            temp_files.extend(iter_res.get('temp_files', []))
            feedback.setProgress(90)

            res = self._create_outputs(iter_res, data, params, paths, parameters, context, feedback, path_eval, _t0)

            elapsed = _time.time() - _t0
            mins, secs = divmod(int(elapsed), 60)
            feedback.setProgress(100)
            feedback.pushInfo(f"Proceso completado en {mins} min {secs:02d} s.")
            return res
        except Exception as e:
            import traceback
            feedback.reportError(traceback.format_exc())
            raise QgsProcessingException(f"Error: {e}")
        finally:
            gc.collect()
            for f in temp_files:
                try:
                    Path(f).unlink(missing_ok=True)
                except Exception as e:
                    feedback.pushDebugInfo(f"[cleanup] No se pudo eliminar {f}: {e}")

    def _setup(self, parameters, context):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        area = self.parameterAsSource(parameters, self.SAMPLING_AREA, context)
        if self.has_as_meters:
            p_dist, radius, min_d, grid_g = (self.parameterAsDistance(parameters, p, context) for p in [
                                             self.PERIMETER_DISTANCE, self.PLOT_RADIUS, self.MIN_SAMPLE_DISTANCE, self.MANUAL_GRID_SPACING])  
        else:
            factor = self.UNIT_CONVERSIONS.get(
                DistanceUnit(
                    self.parameterAsInt(
                        parameters,
                        self.DISTANCE_UNITS,
                        context)),
                1.0)
            p_dist, radius, min_d, grid_g = (self.parameterAsDouble(parameters, p, context) * factor for p in [
                                             self.PERIMETER_DISTANCE, self.PLOT_RADIUS, self.MIN_SAMPLE_DISTANCE, self.MANUAL_GRID_SPACING])  

        _avail = [
            SamplingMethod.SYSTEMATIC_HILBERT,
            SamplingMethod.SIMPLE_RANDOM,
            SamplingMethod.STRATIFIED_HILBERT,
            SamplingMethod.GROUPS_ROW_COL,
            SamplingMethod.STRATIFIED_BY_FIELD,
            SamplingMethod.STRATIFIED_BY_POLYGON]
        if HAS_SKLEARN:
            _avail.insert(4, SamplingMethod.KMEANS_GROUPS)

        m_idx = self.parameterAsInt(parameters, self.METHOD, context)
        m_enum = _avail[m_idx] if m_idx < len(_avail) else SamplingMethod.SYSTEMATIC_HILBERT

        _hp = self.parameterDefinition(self.OUTPUT_HILBERT_PATH)
        if _hp:
            if m_enum == SamplingMethod.GROUPS_ROW_COL:
                _hp.setDescription(self.tr("Salida: Orden NO→SE"))
            else:
                _hp.setDescription(self.tr("Salida: Ruta Hilbert"))

        strata_json_dict = {}
        stratum_field = ""
        polygon_stratum_field = ""

        if m_enum in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
            if m_enum == SamplingMethod.STRATIFIED_BY_FIELD:
                stratum_field = self.parameterAsString(parameters, self.STRATUM_FIELD, context)
            else:
                polygon_stratum_field = self.parameterAsString(parameters, self.POLYGON_STRATUM_FIELD, context)

            raw_json = self.parameterAsString(parameters, self.STRATA_JSON, context)
            try:
                j_data = json.loads(raw_json)
                strata_json_dict = {str(k): int(v) for k, v in j_data.items()}
                s_size = sum(strata_json_dict.values())
            except Exception:
                s_size = self.parameterAsInt(parameters, self.SAMPLE_SIZE, context)
        else:
            s_size = self.parameterAsInt(parameters, self.SAMPLE_SIZE, context)

        params = SamplingParameters(
            project_name=self.parameterAsString(parameters, self.PROJECT_NAME, context),
            perimeter_distance=p_dist, plot_radius=radius, min_sample_distance=min_d,
            hilbert_order=self.parameterAsInt(parameters, self.HILBERT_ORDER, context),
            num_groups=self.parameterAsInt(parameters, self.NUM_GROUPS, context),
            sample_size=s_size,
            num_iterations=self.parameterAsInt(parameters, self.NUM_ITERATIONS, context),
            method=m_enum,
            apply_pullback=self.parameterAsBoolean(parameters, self.APPLY_PULLBACK, context),
            select_inside_buffer=self.parameterAsBoolean(parameters, self.SELECT_INSIDE_BUFFER, context),
            calc_dist_all_points=self.parameterAsBoolean(parameters, self.CALC_DIST_ALL_POINTS, context),
            manual_grid_spacing=grid_g,
            stratum_field=stratum_field,
            polygon_stratum_field=polygon_stratum_field,
            strata_json=strata_json_dict
        )
        html_raw = self.parameterAsString(parameters, self.OUTPUT_HTML_REPORT, context)
        if not html_raw or html_raw == 'TEMPORARY_OUTPUT':
            html_raw = QgsProcessingUtils.generateTempFilename('reporte.html')
        paths = {
            'output_folder': Path(QgsProcessingUtils.tempFolder()),
            'sample_word': self.parameterAsString(parameters, self.SAMPLE_WORD, context) or 'muestra',
            'html_report_path': html_raw,
            'open_report': self.parameterAsBoolean(parameters, self.OPEN_REPORT, context)
        }
        self._html_path_pending = html_raw
        self._html_open_pending = self.parameterAsBoolean(parameters, self.OPEN_REPORT, context)
        return {'source': source, 'area': area}, params, paths

    def _prepare_data(self, layers, params, context, feedback):
        source, area = layers['source'], layers['area']
        if not source:
            return None
        count = source.featureCount()
        if count == 0:
            return None

        s_crs, a_crs = source.sourceCrs(), area.sourceCrs()
        transformer = QgsCoordinateTransform(a_crs, s_crs, QgsProject.instance()) if s_crs != a_crs else None

        poly_strata_geoms = {}
        filter_geometry = QgsGeometry()

        # --- Lógica de Geometría Diferenciada para Estr_Pol-AL vs Normal ---
        if params.method == SamplingMethod.STRATIFIED_BY_POLYGON:
            fld_idx = area.fields().lookupField(params.polygon_stratum_field)
            if fld_idx == -1:
                feedback.reportError(
                    f"AVISO: Campo '{
                        params.polygon_stratum_field}' no encontrado en el área de muestreo.")
                return None

            temp_geoms = defaultdict(list)
            for f in area.getFeatures():
                if not f.hasGeometry() or not GeometryUtils.is_geometry_valid(f.geometry()):
                    continue
                val = str(f.attribute(fld_idx))
                geom = QgsGeometry(f.geometry())
                if transformer:
                    geom.transform(transformer)
                temp_geoms[val].append(geom)

            feedback.pushInfo(
                f"ℹ Estratos (polígonos) detectados en el campo '{
                    params.polygon_stratum_field}': {
                    list(
                        temp_geoms.keys())}")

            for val, geoms in temp_geoms.items():
                u_geom = QgsGeometry.unaryUnion(geoms)
                if u_geom and not u_geom.isEmpty() and not u_geom.isGeosValid():
                    feedback.pushWarning(
                        f"[geom] Estrato '{val}': geometría inválida post-unaryUnion — aplicando makeValid().")
                    u_geom = u_geom.makeValid()
                if params.perimeter_distance > 0:
                    inner = u_geom.buffer(-params.perimeter_distance, 12)
                    u_geom = inner if params.select_inside_buffer else u_geom.difference(inner)
                else:
                    _eps = u_geom.buffer(-0.005, 4)
                    if _eps and not _eps.isNull() and not _eps.isEmpty():
                        u_geom = _eps

                if u_geom and not u_geom.isEmpty():
                    if hasattr(u_geom, 'prepareGeometry'):
                        u_geom.prepareGeometry()
                    poly_strata_geoms[val] = u_geom

            if not poly_strata_geoms:
                feedback.reportError("AVISO: No quedaron polígonos válidos tras aplicar el margen.")
                return None

            filter_geometry = QgsGeometry.unaryUnion(list(poly_strata_geoms.values()))
            if filter_geometry and not filter_geometry.isEmpty() and not filter_geometry.isGeosValid():
                feedback.pushWarning("[geom] Unión global de estratos inválida — aplicando makeValid().")
                filter_geometry = filter_geometry.makeValid()
            filter_area = QgsGeometry(filter_geometry)
            orig_area = QgsGeometry(filter_geometry)
            total_polygons = len(poly_strata_geoms)

        else:
            valid = [f.geometry() for f in area.getFeatures() if GeometryUtils.is_geometry_valid(f.geometry())]
            if not valid:
                return None
            orig_area = QgsGeometry.unaryUnion(valid)
            if orig_area is None or orig_area.isEmpty():
                return None
            if not orig_area.isGeosValid():
                orig_area = orig_area.makeValid()

            filter_area = QgsGeometry(orig_area)
            if transformer:
                filter_area.transform(transformer)

            filter_geometry = QgsGeometry(filter_area)
            if params.perimeter_distance > 0:
                inner = filter_area.buffer(-params.perimeter_distance, 12)
                filter_geometry = inner if params.select_inside_buffer else filter_area.difference(inner)
            else:
                _eps = filter_geometry.buffer(-0.005, 4)
                if _eps and not _eps.isNull() and not _eps.isEmpty():
                    filter_geometry = _eps

            total_polygons = len(valid)

        _geom_prepared = False
        if HAS_SHAPELY:
            try:
                import shapely.wkb as _shwkb
                _sh_filter_geom = _shwkb.loads(bytes(filter_geometry.asWkb()), hex=False)
            except Exception:
                _sh_filter_geom = None
        else:
            _sh_filter_geom = None

        if hasattr(filter_geometry, 'prepareGeometry'):
            filter_geometry.prepareGeometry()
            _geom_prepared = True

        filtered_points_data = []
        seen_coords = set()
        duplicate_count = 0
        grid_multiplier = 1.0 / self.DEDUP_TOLERANCE_M
        request = QgsFeatureRequest().setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        try:
            if HAS_SHAPELY and _sh_filter_geom is not None:
                all_pts_raw = []
                for feature in source.getFeatures(request):
                    if feedback.isCanceled():
                        break
                    if not feature.hasGeometry() or feature.geometry().isEmpty():
                        continue
                    p = feature.geometry().asPoint()
                    all_pts_raw.append((feature.id(), p.x(), p.y()))

                if all_pts_raw:
                    try:
                        xs = _np.array([r[1] for r in all_pts_raw])
                        ys = _np.array([r[2] for r in all_pts_raw])
                        mask = _shapely_contains_xy(
                            _sh_filter_geom, xs, ys) if _shapely_use_v2 else _shvec.contains(
                            _sh_filter_geom, xs, ys)
                        for inside, (fid, px, py) in zip(mask, all_pts_raw):
                            if inside:
                                p_qgs = QgsPointXY(px, py)
                                assigned_stratum = ""
                                if params.method == SamplingMethod.STRATIFIED_BY_POLYGON:
                                    for val, geom in poly_strata_geoms.items():
                                        if geom.contains(p_qgs):
                                            assigned_stratum = val
                                            break
                                    if not assigned_stratum:
                                        continue

                                key = (int(px * grid_multiplier), int(py * grid_multiplier))
                                if key in seen_coords:
                                    duplicate_count += 1
                                else:
                                    seen_coords.add(key)
                                    filtered_points_data.append(
                                        PointData(id=fid, x=px, y=py, stratum_value=assigned_stratum))
                    except Exception:
                        for fid, px, py in all_pts_raw:
                            p_qgs = QgsPointXY(px, py)
                            if filter_geometry.contains(p_qgs):
                                assigned_stratum = ""
                                if params.method == SamplingMethod.STRATIFIED_BY_POLYGON:
                                    for val, geom in poly_strata_geoms.items():
                                        if geom.contains(p_qgs):
                                            assigned_stratum = val
                                            break
                                    if not assigned_stratum:
                                        continue

                                key = (int(px * grid_multiplier), int(py * grid_multiplier))
                                if key in seen_coords:
                                    duplicate_count += 1
                                else:
                                    seen_coords.add(key)
                                    filtered_points_data.append(
                                        PointData(id=fid, x=px, y=py, stratum_value=assigned_stratum))
                scanned = len(all_pts_raw)
            else:
                scanned = 0
                for feature in source.getFeatures(request):
                    scanned += 1
                    if feedback.isCanceled():
                        return None
                    if not feature.hasGeometry() or feature.geometry().isEmpty():
                        continue
                    p = feature.geometry().asPoint()
                    if filter_geometry.contains(p):
                        assigned_stratum = ""
                        if params.method == SamplingMethod.STRATIFIED_BY_POLYGON:
                            for val, geom in poly_strata_geoms.items():
                                if geom.contains(p):
                                    assigned_stratum = val
                                    break
                            if not assigned_stratum:
                                continue

                        key = (int(p.x() * grid_multiplier), int(p.y() * grid_multiplier))
                        if key in seen_coords:
                            duplicate_count += 1
                        else:
                            seen_coords.add(key)
                            filtered_points_data.append(
                                PointData(
                                    id=feature.id(),
                                    x=p.x(),
                                    y=p.y(),
                                    stratum_value=assigned_stratum))
        finally:
            if _geom_prepared and hasattr(filter_geometry, 'releaseCache'):
                try:
                    filter_geometry.releaseCache()
                except Exception as e:
                    if feedback:
                        feedback.pushDebugInfo(f"[releaseCache] filter_geometry: {e}")
            # Liberar también las geometrías GEOS preparadas de los estratos
            # individuales (previene acumulación en STRATIFIED_BY_POLYGON).
            if poly_strata_geoms:
                for _val, _g in poly_strata_geoms.items():
                    if hasattr(_g, 'releaseCache'):
                        try:
                            _g.releaseCache()
                        except Exception as e:
                            if feedback:
                                feedback.pushDebugInfo(f"[releaseCache] estrato '{_val}': {e}")

        if not filtered_points_data:
            return None

        if params.method == SamplingMethod.STRATIFIED_BY_FIELD:
            feedback.pushInfo(f"Extrayendo valores de estrato del campo '{params.stratum_field}'...")
            req = QgsFeatureRequest().setFilterFids([p.id for p in filtered_points_data])
            fld_idx = source.fields().lookupField(params.stratum_field)
            if fld_idx == -1:
                feedback.reportError(f"AVISO: Campo '{params.stratum_field}' no encontrado en la capa de puntos.")
                return None
            req.setSubsetOfAttributes([fld_idx])
            attr_map = {f.id(): f.attribute(fld_idx) for f in source.getFeatures(req)}
            for p in filtered_points_data:
                val = attr_map.get(p.id)
                p.stratum_value = str(val) if val is not None else "NULL"

            unique_found = set(p.stratum_value for p in filtered_points_data)
            feedback.pushInfo(
                f"ℹ Estratos (puntos) detectados en el campo '{
                    params.stratum_field}': {
                    list(unique_found)}")

        sorted_pts = self._sort_points(filtered_points_data, filter_area, params.hilbert_order)
        dist_area = QgsDistanceArea()
        dist_area.setSourceCrs(s_crs, context.transformContext())
        return {
            'sorted_features': sorted_pts,
            'filtered_features': filtered_points_data,
            'original_area': orig_area,
            'area_for_filter': filter_area,
            'source_crs': s_crs,
            'distance_calculator': DistanceCalculator(dist_area),
            'duplicate_count': duplicate_count,
            'initial_point_count': scanned,
            'filtered_point_count': len(filtered_points_data),
            'source_layer_for_hydration': source,
            'area_layer': area,
            'total_polygon_count': total_polygons}

    def _sort_points(self, points: List, area: QgsGeometry, order: int) -> List:
        bbox = area.boundingBox()
        min_x, min_y = bbox.xMinimum(), bbox.yMinimum()
        range_x, range_y = (bbox.width() or 1.0), (bbox.height() or 1.0)
        h_size = 1 << order
        for i, p in enumerate(points):
            p.original_index = i
            nx = int(math.floor((h_size - 1) * (p.x - min_x) / range_x))
            ny = int(math.floor((h_size - 1) * (p.y - min_y) / range_y))
            p.hilbert_index = HilbertCurveCalculator.xy_to_hilbert(
                max(0, min(nx, h_size - 1)), max(0, min(ny, h_size - 1)), order)
        points.sort(key=lambda x: x.hilbert_index)
        return points

    def _sort_points_row_col(self, points: List, area: QgsGeometry) -> List:
        if not points:
            return []
        bbox = area.boundingBox()
        y_max = bbox.yMaximum()
        A = area.area()
        N = max(len(points), 1)
        cell_h = max(math.sqrt(A / N), 1e-6)
        cell_h = self._calc_row_cell_h(points, fallback_h=cell_h)
        pts_copy = list(points)
        pts_copy.sort(key=lambda p: (int((y_max - p.y) / cell_h), -p.y, p.x))
        for seq, p in enumerate(pts_copy):
            p.original_index = seq + 1
        return pts_copy

    @staticmethod
    def _calc_row_cell_h(points, fallback_h: float = 1e-6) -> float:
        ys = sorted(set(round(p.y, 4) for p in points))
        if len(ys) < 2:
            return fallback_h
        all_gaps = [b - a for a, b in zip(ys, ys[1:]) if b > a]
        gaps = [g for g in all_gaps if g >= 1.0]
        if not gaps:
            gaps = sorted(all_gaps)
        if not gaps:
            return fallback_h
        gaps.sort()
        n = len(gaps)
        median_gap = (gaps[n // 2] if n % 2 == 1 else (gaps[n // 2 - 1] + gaps[n // 2]) / 2)
        if median_gap <= 0:
            return fallback_h
        return median_gap * 0.6

    def _run_iterations(self, data, params, parameters, ctx, fb, poly_path):
        sorted_pts = data['sorted_features']
        pool = data['filtered_features']
        res = {
            'generated_files_info': [],
            'nna_results': [],
            'all_rejected_points': [],
            'temp_files_to_delete': [],
            'hilbert_index_map': {
                p.id: i for i,
                p in enumerate(sorted_pts)},
            'strata': {},
            'kmeans_cluster_map': {},
            'kmeans_validation_data': [],
            'hilbert_group_validation_data': [],
            'json_strata_validation_data': [],
            'k_efectivo': 0,
            'rowcol_pts': None}

        strata = {}
        if params.method == SamplingMethod.STRATIFIED_HILBERT:
            k_efectivo = params.num_groups or self._suggest_groups(len(sorted_pts), params.sample_size)
            strata = self._create_strata(len(sorted_pts), k_efectivo, fb)
            res['k_efectivo'] = k_efectivo
            res['strata'] = strata
        elif params.method == SamplingMethod.KMEANS_GROUPS:
            k_km = params.num_groups or self._suggest_groups(len(pool), params.sample_size)
            res['k_efectivo'] = k_km
            if HAS_SKLEARN:
                res['kmeans_cluster_map'] = self._get_kmeans_cluster_map(pool, k_km, fb)
        elif params.method == SamplingMethod.GROUPS_ROW_COL:
            rowcol_pts = self._sort_points_row_col(pool, data['area_for_filter'])
            k_gf = params.num_groups or self._suggest_groups(len(rowcol_pts), params.sample_size)
            strata_gf = self._create_strata(len(rowcol_pts), k_gf, fb)
            res['k_efectivo'] = k_gf
            res['strata'] = strata_gf
            res['rowcol_pts'] = rowcol_pts

        _n_iter = max(params.num_iterations, 1)
        for i in range(params.num_iterations):
            if fb.isCanceled():
                break
            it = i + 1
            fb.setProgress(50 + int(40 * i / _n_iter))

            if params.method == SamplingMethod.SIMPLE_RANDOM:
                samp = self._select_random_sample(pool, params.sample_size, fb)
            elif params.method == SamplingMethod.STRATIFIED_HILBERT:
                samp, grp_pop_h, grp_sam_h = self._select_stratified_sample(
                    sorted_pts, params.sample_size, strata, fb, it)
                prop_h = self._compute_group_proportionality_report(
                    grp_pop_h, grp_sam_h, len(sorted_pts), len(samp), fb)
                res['hilbert_group_validation_data'].append(
                    HilbertGroupValidationData(
                        iteration=it, proportionality=prop_h))
            elif params.method == SamplingMethod.GROUPS_ROW_COL:
                samp, grp_pop_gf, grp_sam_gf = self._select_stratified_sample(
                    rowcol_pts, params.sample_size, strata_gf, fb, it)
                prop_gf = self._compute_group_proportionality_report(
                    grp_pop_gf, grp_sam_gf, len(rowcol_pts), len(samp), fb)
                res['hilbert_group_validation_data'].append(
                    HilbertGroupValidationData(
                        iteration=it, proportionality=prop_gf))
            elif params.method == SamplingMethod.KMEANS_GROUPS:
                samp, grp_pop_km, grp_sam_km, _km = self._select_kmeans_sample(pool, params.sample_size, k_km, fb, it)
                if HAS_SKLEARN:
                    prop_km = self._compute_group_proportionality_report(
                        grp_pop_km, grp_sam_km, len(pool), len(samp), fb)
                    res['kmeans_validation_data'].append(
                        KMeansValidationData(
                            it,
                            self._get_kmeans_cluster_stats(
                                pool,
                                params.sample_size,
                                k_km,
                                fb,
                                it,
                                km_model=_km),
                            proportionality=prop_km))
            elif params.method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
                samp, grp_pop, grp_req, grp_sam = self._select_stratified_by_json_sample(
                    pool, params.strata_json, fb, it)
                res['json_strata_validation_data'].append(
                    {'iteration': it, 'population': grp_pop, 'requested': grp_req, 'sampled': grp_sam})
            else:
                samp = self._select_systematic_sample(sorted_pts, params.sample_size, fb, it)

            if not samp:
                continue
            sc = [copy.deepcopy(p) for p in samp] if params.apply_pullback else list(samp)
            if params.apply_pullback:
                self._apply_pullback_correction(sc, data['original_area'], params.plot_radius, fb)
            if params.min_sample_distance > 0:
                sc, rej = self._enforce_minimum_distance(
                    sc, params.min_sample_distance, data['distance_calculator'], params.sample_size, fb, it)
                res['all_rejected_points'].extend(rej)
            if not sc:
                continue

            nn_temp = self._calc_nn_temp(sc)
            contig = sum(
                1 for d in nn_temp.values() if d and (
                    d < params.manual_grid_spacing or math.isclose(
                        d,
                        params.manual_grid_spacing,
                        rel_tol=1e-4))) if params.manual_grid_spacing > 0 else 0
            nna, tmp = self._run_nna_for_sample(sc, data, parameters, ctx, fb, nn_temp, poly_path)
            if tmp:
                res['temp_files_to_delete'].append(tmp)
            if nna:
                nna['contiguous_count'] = contig
                res['nna_results'].append(nna)
            res['generated_files_info'].append((sc, nn_temp, nna))
        return res

    def _select_stratified_by_json_sample(self, pts: List[PointData], json_sizes: Dict[str, int], fb, it: int):
        strata_pool = defaultdict(list)
        for p in pts:
            strata_pool[p.stratum_value].append(p)

        s = []
        grp_pop = {}
        grp_sam = {}

        for sv, target in json_sizes.items():
            sv_str = str(sv)
            avail = strata_pool.get(sv_str, [])
            grp_pop[sv_str] = len(avail)
            actual = min(target, len(avail))
            grp_sam[sv_str] = actual
            if actual > 0:
                s.extend(random.sample(avail, actual))

        for sv, p_list in strata_pool.items():
            if sv not in grp_pop:
                grp_pop[sv] = len(p_list)
                grp_sam[sv] = 0

        return s, grp_pop, json_sizes, grp_sam

    def _run_nna_for_sample(self, sample, data, parameters, ctx, fb, nn_map, poly_path):
        tmp = None
        n = len(sample)
        if n < 2:
            return None, tmp
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
        d_esp = 0.5 / math.sqrt(n / A)
        ivmc_base = d_obs / d_esp if d_esp > 0 else 0
        ab = area_geom.boundingBox().area()
        ivmc_corr = ivmc_base * math.sqrt(ab / A) if ab > 0 else ivmc_base
        try:
            z = (d_obs - d_esp) / (0.26136 / math.sqrt(n * n / A))
        except Exception:
            z = 0.0

        polys_count = 0
        try:
            area_layer = data.get('area_layer')
            if area_layer:
                crs_transform = QgsCoordinateTransform(
                    data['source_crs'], area_layer.sourceCrs(), QgsProject.instance())
                poly_feats = [f for f in area_layer.getFeatures() if f.hasGeometry()]
                poly_dict = {f.id(): f for f in poly_feats}
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
        except Exception:
            pass

        return {'indice_nn': ivmc_corr, 'n_puntos': n, 'polygons_hit_count': polys_count,
                'OBSERVED_MD': d_obs, 'EXPECTED_MD': d_esp, 'Z_SCORE': z}, tmp

    def _evaluate_path(self, points, area, params, calc, dups, fb):
        if len(points) < 2:
            return {}
        dists = [calc.calculate_distance(QgsPointXY(points[i - 1].x, points[i - 1].y),
                                         QgsPointXY(points[i].x, points[i].y)) for i in range(1, len(points))]
        tot = sum(dists)
        mean = tot / len(dists)
        h = math.sqrt(area.area() / len(points)) if len(points) > 0 else 0
        var = sum([(d - mean)**2 for d in dists]) / (len(dists) - 1) if len(dists) > 1 else 0
        std = math.sqrt(var)

        h_idxs = set()
        bbox = area.boundingBox()
        mx, my = bbox.xMinimum(), bbox.yMinimum()
        rx, ry = (bbox.width() or 1), (bbox.height() or 1)
        h_size = 1 << params.hilbert_order
        for p in points:
            nx = int((h_size - 1) * (p.x - mx) / rx)
            ny = int((h_size - 1) * (p.y - my) / ry)
            h_idxs.add(HilbertCurveCalculator.xy_to_hilbert(nx, ny, params.hilbert_order))
        coll = len(points) - len(h_idxs)
        g = params.manual_grid_spacing if params.manual_grid_spacing > 0 else 0
        perim = area.length()
        compac = (4 * math.pi * area.area()) / (perim**2) if perim > 0 else 0
        return {"efficiency_index_R": mean / h if h > 0 else 0,
                "cv": std / mean if mean > 0 else 0,
                "order": params.hilbert_order,
                "collisions": coll,
                "collision_ratio": (coll / len(points)) * 100,
                "area_km2": area.area() / 1e6,
                "mesh_distance_h": h,
                "grid_spacing_g": g,
                "num_points": len(points),
                "total_length_km": tot / 1000,
                "std_dev_norm": std / h if h > 0 else 0,
                "length_norm": tot / (len(points) * h) if h > 0 else 0,
                "indice_compacidad": compac,
                "coverage_index": h / g if g > 0 else 0,
                "duplicate_count": dups}

    def _select_systematic_sample(self, pts, k, fb, it):
        n = len(pts)
        if n == 0:
            return []
        if k >= n:
            return list(pts)
        step = n / k
        start = random.uniform(0, step)
        indices = []
        for i in range(k):
            val = int(round(start + i * step))
            if val >= n:
                val = n - 1
            if val < 0:
                val = 0
            indices.append(val)
        return [pts[i] for i in sorted(list(set(indices)))]

    def _select_random_sample(self, pts, k, fb): return random.sample(pts, min(len(pts), k))

    def _select_stratified_sample(self, pts: List, k: int, strata: Dict,
                                  fb: QgsProcessingFeedback, it: int) -> Tuple[List, Dict[int, int], Dict[int, int]]:
        s, group_population, group_sample = [], {}, {}
        base, rem = divmod(min(k, len(pts)), len(strata))
        keys = list(strata.keys())
        random.shuffle(keys)
        for i, key in enumerate(keys):
            info = strata[key]
            sub = pts[info.start: info.end]
            sz = base + (1 if i < rem else 0)
            group_population[key] = len(sub)
            assigned = min(sz, len(sub))
            group_sample[key] = assigned
            if sub:
                s.extend(random.sample(sub, assigned))
        return s, group_population, group_sample

    def _select_kmeans_sample(self, pts: List, k: int, g: int, fb: QgsProcessingFeedback,
                              it: int) -> Tuple[List, Dict[int, int], Dict[int, int]]:
        if not HAS_SKLEARN:
            return [], {}, {}
        n_real = min(g, len(pts))
        coords = [[p.x, p.y] for p in pts]
        km = KMeans(n_clusters=n_real, random_state=it).fit(coords)
        clus = {i: [] for i in range(n_real)}
        for p, lbl in zip(pts, km.labels_):
            clus[lbl].append(p)
        clus = {c_id: c_pts for c_id, c_pts in clus.items() if c_pts}
        total = len(pts)
        group_population = {c_id: len(c_pts) for c_id, c_pts in clus.items()}
        group_sample = {}
        alloc = []
        for c_id, c_pts in clus.items():
            ideal = (len(c_pts) / total) * k
            base = int(ideal)
            alloc.append((c_id, c_pts, base, ideal - base))
        alloc.sort(key=lambda x: x[3], reverse=True)
        rem = k - sum(x[2] for x in alloc)
        s = []
        for c_id, c_pts, base, _ in alloc:
            sz = base + (1 if rem > 0 else 0)
            rem -= 1
            assigned = min(len(c_pts), sz) if sz > 0 else 0
            group_sample[c_id] = assigned
            if assigned > 0:
                s.extend(random.sample(c_pts, assigned))
        for c_id in group_population:
            group_sample.setdefault(c_id, 0)
        return s, group_population, group_sample, km

    def _create_strata(self, tot, num, fb):
        s = {}
        base, rem = divmod(tot, num)
        st = 0
        for i in range(num):
            sz = base + 1 if i < rem else base
            if sz > 0:
                s[i] = StratumInfo(st, st + sz, sz)
                st += sz
        return s

    def _get_kmeans_cluster_map(self, pts: List, g: int, fb: QgsProcessingFeedback) -> Dict:
        if not HAS_SKLEARN or not pts:
            return {}
        try:
            km = KMeans(n_clusters=min(g, len(pts)), random_state=42, n_init='auto').fit([[p.x, p.y] for p in pts])
            return {p.id: int(lbl) for p, lbl in zip(pts, km.labels_)}
        except Exception:
            return {}

    def _get_kmeans_cluster_stats(
            self,
            pts: List,
            k: int,
            g: int,
            fb: QgsProcessingFeedback,
            it: int,
            km_model=None) -> List[Dict]:
        if not HAS_SKLEARN or not pts:
            return []
        try:
            n_real = min(g, len(pts))
            km = km_model if km_model is not None else KMeans(
                n_clusters=n_real, random_state=it, n_init='auto').fit([[p.x, p.y] for p in pts])
            stats = []
            for cid in range(n_real):
                members = [pts[i] for i, lbl in enumerate(km.labels_) if lbl == cid]
                if not members:
                    continue
                cx, cy = km.cluster_centers_[cid]
                dists = [math.hypot(p.x - cx, p.y - cy) for p in members]
                stats.append({'cluster_id': cid, 'size': len(members), 'center_x': cx, 'center_y': cy,
                             'mean_dist_to_center': sum(dists) / len(dists), 'max_dist_to_center': max(dists)})
            return stats
        except Exception:
            return []

    def _compute_group_proportionality_report(self,
                                              group_population: Dict[int,
                                                                     int],
                                              group_sample: Dict[int,
                                                                 int],
                                              total_pop: int,
                                              total_sample: int,
                                              fb: QgsProcessingFeedback) -> Dict:
        rows, abs_diffs, groups_zero = [], [], []
        for gid, pop_n in sorted(group_population.items()):
            sam_n = group_sample.get(gid, 0)
            pct_pop = (pop_n / total_pop) * 100.0 if total_pop > 0 else 0.0
            pct_sam = (sam_n / total_sample) * 100.0 if total_sample > 0 else 0.0
            diff = pct_pop - pct_sam
            abs_diffs.append(abs(diff))
            if sam_n == 0:
                groups_zero.append(gid)
            rows.append({'group_id': gid, 'pop_n': pop_n, 'sam_n': sam_n, 'pct_pop': round(pct_pop, 4),
                        'pct_sam': round(pct_sam, 4), 'diff': round(diff, 4), 'abs_diff': round(abs(diff), 4)})
        return {
            'rows': rows,
            'max_diff': round(
                max(abs_diffs),
                4) if abs_diffs else 0.0,
            'mean_diff': round(
                sum(abs_diffs)
                / len(abs_diffs),  # noqa: W503
                4) if abs_diffs else 0.0,
            'groups_zero': groups_zero}

    def _apply_pullback_correction(self, pts: List, area: QgsGeometry, rad: float, fb: QgsProcessingFeedback) -> None:
        if rad <= 0:
            return
        try:
            safe = area.buffer(-rad, 12)
        except Exception:
            return
        if safe is None or safe.isEmpty():
            return
        for p in pts:
            if not safe.contains(p.get_geometry()):
                try:
                    p.set_geometry(safe.nearestPoint(p.get_geometry()))
                    p.was_corrected = True
                except Exception:
                    pass

    def _enforce_minimum_distance(self, pts, d, calc, k, fb, it):
        random.shuffle(pts)
        k_list, r_list = [], []
        idx = QgsSpatialIndex()
        accepted_geoms = {}
        for p in pts:
            pt_geom = p.get_geometry()
            neighbors = idx.nearestNeighbor(pt_geom, 1)
            too_close = False
            if neighbors:
                nn_geom = accepted_geoms.get(neighbors[0])
                if nn_geom is not None:
                    if nn_geom.distance(pt_geom) < d:
                        too_close = True
            if not too_close:
                f = QgsFeature()
                f.setId(p.id)
                f.setGeometry(pt_geom)
                idx.addFeature(f)
                accepted_geoms[p.id] = pt_geom
                k_list.append(p)
            else:
                p.rejected_in_iteration = it
                r_list.append(p)
        return k_list, r_list

    def _calc_nn_temp(self, pts):
        if not pts:
            return {}
        if len(pts) == 1:
            return {pts[0].id: None}
        idx = QgsSpatialIndex()
        feats_tmp = {}
        for p in pts:
            f = QgsFeature()
            f.setId(p.id)
            f.setGeometry(p.get_geometry())
            idx.addFeature(f)
            feats_tmp[p.id] = p.get_geometry()
        m = {}
        for p in pts:
            neighbors = idx.nearestNeighbor(p.get_geometry(), 2)
            nn_id = next((nid for nid in neighbors if nid != p.id), None)
            m[p.id] = feats_tmp[nn_id].distance(p.get_geometry()) if (
                nn_id is not None and nn_id in feats_tmp) else None
        return m

    def _suggest_groups(self, tot: int, samp: int) -> int:
        base = max(2, math.ceil(math.sqrt(samp)))
        return min(base, max(2, int(math.sqrt(tot))))

    def _create_outputs(self, res, data, params, paths, parameters, ctx, fb, path_eval, t0: float = 0.0):
        source = data['source_layer_for_hydration']
        fields = self._fields(source.fields(), method=params.method)
        idx_map = {p.id: p.original_index for p in res['rowcol_pts']} if params.method == SamplingMethod.GROUPS_ROW_COL and res.get(  
            'rowcol_pts') else res['hilbert_index_map']
        boundary = GeometryUtils.get_polygon_boundary_as_line(data['area_for_filter'], fb)

        all_ids = set()
        for s, *_ in res['generated_files_info']:
            all_ids.update(p.id for p in s)
        if res['all_rejected_points']:
            all_ids.update(p.id for p in res['all_rejected_points'])

        # --- VALIDACIÓN DE MUESTRA VACÍA ---
        if not all_ids and not res.get('rowcol_pts'):
            fb.reportError("❌ ERROR CRÍTICO: La muestra final tiene 0 puntos. No se generarán capas ni reportes HTML.")
            if params.method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
                fb.reportError(
                    "💡 Causa probable: Las claves introducidas en tu diccionario JSON NO existen en la tabla de atributos.")  
                fb.reportError(f"   -> Buscaste en JSON: {list(params.strata_json.keys())}")
                fb.reportError("   -> Revisa el Log azul más arriba para ver qué valores encontró realmente el algoritmo.")  
            return {}
        # -----------------------------------

        req = QgsFeatureRequest().setFilterFids(list(all_ids))
        fmap = {f.id(): QgsFeature(f) for f in source.getFeatures(req)}
        fb.setProgress(91)

        self._sample_layers_pending = []
        self._current_method = params.method
        self._has_rejected = bool(res.get('all_rejected_points'))

        _n_files = max(len(res['generated_files_info']), 1)  # noqa: F841
        for i, entry in enumerate(res['generated_files_info']):
            samp, nn = entry[0], entry[1]
            nna = entry[2] if len(entry) > 2 else None
            it = i + 1
            ivmc = nna.get('indice_nn') if nna else None
            cont = nna.get('contiguous_count') if nna else None
            nm = f"{self.METHOD_ABBREVIATIONS[params.method]}_{paths['sample_word']}_{it:02d}"
            if params.apply_pullback:
                nm += "_RETR"
            if cont is not None:
                nm += f"_Cont_{cont}"
            if ivmc:
                nm += f"_IVMC_{f'{ivmc:.3f}'.replace('.', ',')}"
            if nna is not None:
                nna['muestra'] = nm

            hyd_samp = []
            for p in samp:
                if p.id in fmap:
                    f = QgsFeature(fmap[p.id])
                    f.setGeometry(p.get_geometry())
                    hyd_samp.append(
                        SamplePoint(
                            f,
                            QgsGeometry(
                                f.geometry()),
                            p.was_corrected,
                            p.correction_distance,
                            0,
                            p.stratum_value))

            tmp = QgsProcessingUtils.generateTempFilename(f"{nm}.gpkg")
            ok = self._write_sample_file(
                hyd_samp,
                Path(tmp),
                fields,
                idx_map,
                it,
                params.method,
                boundary,
                data['source_crs'],
                ctx,
                fb,
                params,
                ivmc,
                nn_map=nn)
            if ok:
                self._sample_layers_pending.append((tmp, nm))

        fb.setProgress(92)
        if paths['html_report_path']:
            gen = ReportGenerator(
                res['nna_results'],
                params.project_name,
                params,
                self.METHOD_NAMES,
                path_eval,
                data['initial_point_count'],
                data['filtered_point_count'],
                res['kmeans_validation_data'],
                hilbert_group_val=res['hilbert_group_validation_data'],
                json_strata_val=res.get('json_strata_validation_data'),
                total_polygon_count=data.get(
                    'total_polygon_count',
                    0),
                t0=t0)
            gen.write(paths['html_report_path'], fb, paths['open_report'])
            self._report_gen = gen
        fb.setProgress(95)

        result = {}
        if res['all_rejected_points']:
            hyd_rej = []
            for p in res['all_rejected_points']:
                if p.id in fmap:
                    f = QgsFeature(fmap[p.id])
                    f.setGeometry(p.get_geometry())
                    hyd_rej.append(
                        SamplePoint(
                            f,
                            QgsGeometry(
                                f.geometry()),
                            rejected_in_iteration=p.rejected_in_iteration))
            self._write_rejected(hyd_rej, source, parameters, ctx, result)

        _pts_out = res['rowcol_pts'] if params.method == SamplingMethod.GROUPS_ROW_COL and res.get(
            'rowcol_pts') else data['sorted_features']
        self._create_all_points(
            parameters,
            ctx,
            _pts_out,
            source,
            boundary,
            params.calc_dist_all_points,
            data['source_crs'],
            result,
            res['strata'],
            res['kmeans_cluster_map'],
            fb=fb,
            method=params.method)
        fb.setProgress(97)

        if params.method in (
                SamplingMethod.SYSTEMATIC_HILBERT,
                SamplingMethod.STRATIFIED_HILBERT,
                SamplingMethod.GROUPS_ROW_COL):
            _pts_path = res['rowcol_pts'] if params.method == SamplingMethod.GROUPS_ROW_COL and res.get(
                'rowcol_pts') else data['sorted_features']
            self._create_hilbert_path(
                parameters,
                params,
                ctx,
                _pts_path,
                data['source_crs'],
                path_eval,
                result,
                fb=fb,
                clip_geom=data['area_for_filter'])

        self._assign_post_processors(ctx)
        return result

    def _fields(self, src: QgsFields, method=None) -> QgsFields:
        f = QgsFields(src)
        _idx_field = self.FIELD_NWSE_IDX if method == SamplingMethod.GROUPS_ROW_COL else self.FIELD_HILBERT_IDX
        for n in [self.FIELD_SAMPLE_ID, _idx_field, self.FIELD_ITER_NUM]:
            f.append(QgsField(n, _INT_TYPE))
        for n in [self.FIELD_METHOD, self.FIELD_CORRECTED, self.FIELD_NEIGHBOR_CHECK]:
            f.append(QgsField(n, _STR_TYPE))
        for n in [self.FIELD_NEAREST_NEIGHBOR, self.FIELD_SAMPLE_IVMC, self.FIELD_DIST_BORDER]:
            f.append(QgsField(n, _DBL_TYPE))
        if method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
            f.append(QgsField(self.FIELD_STRATUM, _STR_TYPE))
        return f

    def _create_output_fields(self, src, dist, method=None):
        f = QgsFields(src)
        f.append(QgsField(self.FIELD_NWSE_IDX if method ==  # noqa: W504
                 SamplingMethod.GROUPS_ROW_COL else self.FIELD_HILBERT_IDX, _INT_TYPE))
        if dist:
            f.append(QgsField(self.FIELD_DIST_BORDER, _DBL_TYPE))
        if method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
            f.append(QgsField(self.FIELD_STRATUM, _STR_TYPE))
        return f

    def _write_sample_file(self, pts, path, flds, imap, it, method, boundary, crs, ctx, fb, params, ivmc, nn_map=None):
        opts = QgsVectorFileWriter.SaveVectorOptions()
        opts.driverName = "GPKG"
        opts.layerName = path.stem
        w = QgsVectorFileWriter.create(str(path), flds, QgsWkbTypes.Point, crs, ctx.transformContext(), opts)
        if w.hasError() != QgsVectorFileWriter.NoError:
            return False

        pts.sort(key=lambda sp: imap.get(sp.feature.id(), float('inf')))
        geoms = [p.feature.geometry() for p in pts]
        col_map = {field.name(): i for i, field in enumerate(flds)}

        for i, p in enumerate(pts):
            if nn_map is not None:
                final_nn = nn_map.get(p.feature.id()) or 0
            else:
                curr = geoms[i]
                min_d = float('inf')
                for j, other in enumerate(geoms):
                    if i == j:
                        continue
                    d = curr.distance(other)
                    if d < min_d:
                        min_d = d
                final_nn = min_d if min_d != float('inf') else 0

            n_check = "NULL"
            if params.manual_grid_spacing > 0:
                is_neighbor = final_nn < params.manual_grid_spacing or math.isclose(
                    final_nn, params.manual_grid_spacing, rel_tol=1e-4)
                n_check = "Verdadero" if is_neighbor else "Falso"
            elif params.manual_grid_spacing == 0:
                n_check = "N/A (g=0)"

            dist_border = 0
            if boundary and p.feature.geometry():
                try:
                    dist_border = float(p.feature.geometry().distance(boundary))
                except Exception:
                    pass

            attrs = [None] * len(flds)
            for field in p.feature.fields():
                idx = col_map.get(field.name())
                if idx is not None:
                    attrs[idx] = p.feature.attribute(field.name())

            if self.FIELD_SAMPLE_ID in col_map:
                attrs[col_map[self.FIELD_SAMPLE_ID]] = i + 1
            if self.FIELD_HILBERT_IDX in col_map:
                attrs[col_map[self.FIELD_HILBERT_IDX]] = imap.get(p.feature.id(), -1)
            if self.FIELD_NWSE_IDX in col_map:
                attrs[col_map[self.FIELD_NWSE_IDX]] = imap.get(p.feature.id(), -1)
            if self.FIELD_ITER_NUM in col_map:
                attrs[col_map[self.FIELD_ITER_NUM]] = it
            if self.FIELD_METHOD in col_map:
                attrs[col_map[self.FIELD_METHOD]] = self.METHOD_ABBREVIATIONS[method]
            if self.FIELD_CORRECTED in col_map:
                attrs[col_map[self.FIELD_CORRECTED]] = "Verdadero" if p.was_corrected else "Falso"
            if self.FIELD_NEAREST_NEIGHBOR in col_map:
                attrs[col_map[self.FIELD_NEAREST_NEIGHBOR]] = round(final_nn, 2) if final_nn else 0.0
            if self.FIELD_NEIGHBOR_CHECK in col_map:
                attrs[col_map[self.FIELD_NEIGHBOR_CHECK]] = n_check
            if self.FIELD_SAMPLE_IVMC in col_map and ivmc:
                attrs[col_map[self.FIELD_SAMPLE_IVMC]] = round(float(ivmc), 4)
            if self.FIELD_DIST_BORDER in col_map:
                attrs[col_map[self.FIELD_DIST_BORDER]] = round(float(dist_border), 2)
            if self.FIELD_STRATUM in col_map:
                attrs[col_map[self.FIELD_STRATUM]] = p.stratum_value

            f = QgsFeature()
            f.setGeometry(p.feature.geometry())
            f.setAttributes(attrs)
            w.addFeature(f)
        del w
        return True

    def _write_rejected(self, pts, src, params, ctx, res):
        flds = QgsFields(src.fields())
        flds.append(QgsField(self.FIELD_REJECTED_ITER, _INT_TYPE))
        sink, did = self.parameterAsSink(params, self.OUTPUT_REJECTED, ctx, flds, src.wkbType(), src.sourceCrs())
        if sink:
            for p in pts:
                f = QgsFeature(flds)
                f.setGeometry(p.feature.geometry())
                for a in src.fields():
                    f.setAttribute(a.name(), p.feature.attribute(a.name()))
                f.setAttribute(self.FIELD_REJECTED_ITER, p.rejected_in_iteration)
                sink.addFeature(f, QgsFeatureSink.FastInsert)
            res[self.OUTPUT_REJECTED] = did

    def _assign_post_processors(self, ctx: QgsProcessingContext) -> None:
        framework_names = {self.tr('Salida: Puntos filtrados'), 'Salida: Puntos filtrados'}
        hilbert_names = set()
        if self._current_method in (
                SamplingMethod.SYSTEMATIC_HILBERT,
                SamplingMethod.STRATIFIED_HILBERT,
                SamplingMethod.GROUPS_ROW_COL):
            hilbert_names = {
                self.tr('Salida: Ruta Hilbert'),
                'Salida: Ruta Hilbert',
                self.tr('Salida: Orden NO→SE'),
                'Salida: Orden NO→SE'}
            framework_names |= hilbert_names
        framework_names |= {self.tr('Salida: Rechazados'), 'Salida: Rechazados'}
        self._framework_layer_names = framework_names
        expected = 1
        if self._current_method in (
                SamplingMethod.SYSTEMATIC_HILBERT,
                SamplingMethod.STRATIFIED_HILBERT,
                SamplingMethod.GROUPS_ROW_COL) and self._hilbert_path_created:
            expected += 1
        if self._has_rejected:
            expected += 1
        try:
            self._connect_layers_added_signal(framework_names, list(self._sample_layers_pending), expected)
        except Exception:
            pass

    def _create_all_points(self, params, ctx, pts, src, bnd, dist, crs, res, strata, kmap, fb=None, method=None):
        flds = self._create_output_fields(src.fields(), dist, method=method)
        if strata or kmap:
            flds.append(QgsField(self.FIELD_GROUP_ID, _INT_TYPE))
        sink, did = self.parameterAsSink(params, self.OUTPUT_ALL, ctx, flds, QgsWkbTypes.Point, crs)
        if sink:
            imap = {p.id: i for i, p in enumerate(pts)}
            imap_nwse = {p.id: p.original_index for p in pts}
            fid_to_group = {}
            if strata:
                pos_to_stratum = {}
                for sid, info in strata.items():
                    for pos in range(info.start, info.end):
                        pos_to_stratum[pos] = sid + 1
                for pt in pts:
                    pos = imap.get(pt.id, -1)
                    if pos >= 0 and pos in pos_to_stratum:
                        fid_to_group[pt.id] = pos_to_stratum[pos]
            if kmap:
                for fid, cluster_id in kmap.items():
                    fid_to_group[fid] = cluster_id + 1

            req = QgsFeatureRequest().setFilterFids([p.id for p in pts])
            col_map = {field.name(): idx for idx, field in enumerate(flds)}
            has_border, has_group_id, has_stratum = self.FIELD_DIST_BORDER in col_map, self.FIELD_GROUP_ID in col_map, self.FIELD_STRATUM in col_map  

            stratum_map = {p.id: p.stratum_value for p in pts} if method in (
                SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON) else {}

            for f in src.getFeatures(req):
                nf = QgsFeature(flds)
                nf.initAttributes(len(flds))
                nf.setGeometry(f.geometry())
                for a in src.fields():
                    if a.name() in col_map:
                        nf.setAttribute(a.name(), f.attribute(a.name()))
                if method == SamplingMethod.GROUPS_ROW_COL:
                    nf.setAttribute(self.FIELD_NWSE_IDX, imap_nwse.get(f.id(), -1))
                else:
                    nf.setAttribute(self.FIELD_HILBERT_IDX, imap.get(f.id(), -1))
                if has_border and bnd and f.hasGeometry():
                    try:
                        nf.setAttribute(self.FIELD_DIST_BORDER, round(float(f.geometry().distance(bnd)), 2))
                    except Exception:
                        pass
                if has_group_id and fid_to_group.get(f.id()) is not None:
                    nf.setAttribute(self.FIELD_GROUP_ID, fid_to_group.get(f.id()))
                if has_stratum:
                    nf.setAttribute(self.FIELD_STRATUM, stratum_map.get(f.id(), ""))
                sink.addFeature(nf, QgsFeatureSink.FastInsert)
            res[self.OUTPUT_ALL] = did

    def _create_hilbert_path(self, parameters, params, ctx, pts, crs, path_eval, res, fb=None, clip_geom=None):
        flds = QgsFields()
        flds.append(QgsField("id", _INT_TYPE))
        flds.append(QgsField("descripcion", _STR_TYPE))
        flds.append(QgsField("fila", _INT_TYPE))
        flds.append(QgsField("n_puntos", _INT_TYPE))
        flds.append(QgsField("método", _STR_TYPE))
        sink, did = self.parameterAsSink(parameters, self.OUTPUT_HILBERT_PATH, ctx, flds, QgsWkbTypes.LineString, crs)
        if not sink:
            return
        _is_gf = params.method == SamplingMethod.GROUPS_ROW_COL
        _mname = self.METHOD_NAMES.get(params.method, "")

        if _is_gf:
            if not pts:
                res[self.OUTPUT_HILBERT_PATH] = did
                return
            A_m2 = (path_eval.get("area_km2") or 0) * 1e6
            N_pts = max(len(pts), 1)
            cell_h = self._calc_row_cell_h(pts, fallback_h=max(math.sqrt(A_m2 / N_pts), 1e-6) if A_m2 > 0 else 1e-6)
            y_max = max(p.y for p in pts)
            row_map = defaultdict(list)
            for p in pts:
                row_map[int((y_max - p.y) / cell_h)].append(p)

            for feat_id, (row_id, row_pts) in enumerate(sorted(row_map.items()), start=1):
                row_pts.sort(key=lambda p: (-p.y, p.x))
                if len(row_pts) < 2:
                    p0 = row_pts[0]
                    coords = [QgsPointXY(p0.x, p0.y), QgsPointXY(p0.x + (cell_h * 0.1 if cell_h > 0 else 1.0), p0.y)]
                else:
                    coords = [QgsPointXY(p.x, p.y) for p in row_pts]
                line = QgsGeometry.fromPolylineXY(coords)
                if line and not line.isEmpty():
                    if clip_geom is not None:
                        try:
                            line = line.intersection(clip_geom)
                        except Exception:
                            pass
                    if not line or line.isEmpty():
                        continue
                    parts = line.asGeometryCollection() if not line.isEmpty() else []
                    if not parts:
                        parts = [line]
                    for part in parts:
                        if part.isEmpty() or part.type() != QgsWkbTypes.LineGeometry:
                            continue
                        f = QgsFeature(flds)
                        f.setGeometry(part)
                        f.setAttribute(0, feat_id)
                        f.setAttribute(1, f"Fila {row_id}")
                        f.setAttribute(2, row_id)
                        f.setAttribute(3, len(row_pts))
                        f.setAttribute(4, _mname)
                        sink.addFeature(f, QgsFeatureSink.FastInsert)
        else:
            line = QgsGeometry.fromPolylineXY([QgsPointXY(p.x, p.y) for p in pts])
            if clip_geom is not None:
                try:
                    line = line.intersection(clip_geom)
                except Exception:
                    pass
            if line and not line.isEmpty():
                parts = line.asGeometryCollection() if not line.isEmpty() else [line]
                if not parts:
                    parts = [line]
                for part_id, part in enumerate(parts, start=1):
                    if part.isEmpty() or part.type() != QgsWkbTypes.LineGeometry:
                        continue
                    f = QgsFeature(flds)
                    f.setGeometry(part)
                    f.setAttribute(0, part_id)
                    f.setAttribute(1, "Curva de Hilbert")
                    f.setAttribute(2, params.hilbert_order)
                    f.setAttribute(3, len(pts))
                    f.setAttribute(4, _mname)
                    sink.addFeature(f, QgsFeatureSink.FastInsert)
        res[self.OUTPUT_HILBERT_PATH] = did
        self._hilbert_path_created = True

# --- Reporte HTML ---


class ReportGenerator:
    def __init__(
            self,
            results,
            project,
            params,
            methods,
            path_eval,
            init,
            filt,
            k_val=None,
            hilbert_group_val=None,
            json_strata_val=None,
            total_polygon_count=0,
            t0: float = 0.0):
        self.res = results
        self.proj = project
        self.params = params
        self.meth = methods
        self.eval = path_eval
        self.init = init
        self.filt = filt
        self.k_val = k_val
        self.hilbert_group_val: List = hilbert_group_val or []
        self.json_strata_val: List = json_strata_val or []
        self.total_polygon_count: int = total_polygon_count
        self._t0 = t0

    def _fmt_elapsed(self, secs: float) -> str:
        s = int(round(secs))
        if s < 60:
            return f"{s} s"
        m, s = divmod(s, 60)
        if m < 60:
            return f"{m} min {s:02d} s"
        h, m = divmod(m, 60)
        return f"{h} h {m:02d} min {s:02d} s"

    def write(self, path, fb, open_it):
        if not path:
            return
        if path == 'TEMPORARY_OUTPUT':
            path = QgsProcessingUtils.generateTempFilename('reporte.html')
        self._html_path = path
        self._html_fb = fb
        self._html_open = open_it
        try:
            self._elapsed_at_write = 0.0
            with open(path, 'w', encoding='utf-8') as f:
                f.write(self._html())
        except Exception as e:
            if fb:
                fb.pushWarning(f"[ReportGenerator.write] No se pudo escribir el reporte HTML en '{path}': {e}")

    def rewrite_with_time(self, elapsed: float) -> None:
        self._elapsed_at_write = elapsed
        try:
            with open(self._html_path, 'w', encoding='utf-8') as f:
                f.write(self._html())
            if self._html_open:
                webbrowser.open(Path(self._html_path).as_uri())
        except Exception as e:
            # Sin feedback disponible aquí; usar logger del módulo.
            import logging
            logging.getLogger('MuestreoEspacial').warning(
                f"[ReportGenerator.rewrite_with_time] Fallo al reescribir/abrir HTML '{self._html_path}': {e}"
            )

    def _html(self):
        _r_idx = self.eval.get('efficiency_index_R', 0.0) if self.eval else 0.0
        _coll_h = self.eval.get('collision_ratio', 100.0) if self.eval else 100.0
        _malla_reg = (0.85 <= _r_idx <= 1.35 and _coll_h < 10.0)
        _malla_frag = (not _malla_reg and _coll_h < 10.0 and _r_idx > 1.35)
        _is_disp = (_malla_reg or _malla_frag or self.params.min_sample_distance > 0 or self.params.method in (
            SamplingMethod.SYSTEMATIC_HILBERT, SamplingMethod.STRATIFIED_HILBERT, SamplingMethod.GROUPS_ROW_COL))
        if self.params.method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
            _is_disp = False

        if _malla_reg:
            grid_label = "Regular (rectangular / hexagonal)"
        elif _malla_frag:
            grid_label = "Regular fragmentada (malla sobre multipolígono)"
        else:
            grid_label = "Irregular (nube de puntos)"
        _any_regular = _malla_reg or _malla_frag

        import math as _mth

        def _rej(r):
            o_ = r.get('n_puntos', 0)
            s_ = self.params.sample_size
            return ((s_ - o_) / s_ * 100) if s_ > 0 else 0.0

        def _ivn(r): return (r.get('indice_nn') or 0) / _mth.sqrt(max(r.get('n_puntos') or 1, 1))

        _valid_top = [r for r in self.res if r.get('indice_nn') is not None]
        _n_med = sum(r.get('n_puntos', 0) for r in _valid_top) / len(_valid_top) if _valid_top else 0
        _alto = _n_med > 5000

        if _valid_top:
            if _alto:
                _pvars = len(set(r.get('polygons_hit_count', 0) for r in _valid_top)) > 1
                if _pvars:
                    _t3 = sorted(
                        _valid_top,
                        key=lambda x: (
                            x.get('polygons_hit_count') or 0,
                            _ivn(x)),
                        reverse=True)[
                        :3]
                else:
                    _t3 = sorted(_valid_top, key=_ivn, reverse=True)[:3]
            elif _is_disp:
                _t3 = sorted(_valid_top, key=lambda x: (-(x.get('polygons_hit_count') or 0), _rej(x),
                             (x.get('contiguous_count') or 0), -(x.get('indice_nn') or 0)))[:3]
            else:
                _opt = [r for r in _valid_top if MuestreoEspacialPuntos.NNI_RANDOM_LO <=  # noqa: W504
                        (r.get('indice_nn') or 0) <= MuestreoEspacialPuntos.NNI_RANDOM_HI]
                _pool = _opt if _opt else sorted(_valid_top, key=lambda x: abs((x.get('indice_nn') or 99) - 1.0))[:3]
                _t3 = sorted(_pool, key=lambda x: (-(x.get('polygons_hit_count') or 0),
                             _rej(x), (x.get('contiguous_count') or 0)))[:3]
        else:
            _t3 = []

        _top3_names = {r.get('muestra') for r in _t3}

        rows = ""
        for r in self.res:
            nni = r.get('indice_nn')
            s = self.params.sample_size
            o = r.get('n_puntos', 0)
            rej = ((s - o) / s) * 100 if s > 0 else 0
            if nni is None:
                pat, cal = "N/A", "N/A"
            elif _any_regular and nni <= MuestreoEspacialPuntos.NNI_RANDOM_HI:
                pat, cal = "Agrupado / Aleatorio", "Problema grave ⚠"
            elif nni < 0.8:
                pat, cal = "Agrupado Fuerte", "No Aceptable"
            elif nni < MuestreoEspacialPuntos.NNI_RANDOM_LO:
                pat, cal = "Agrupado Moderado", "No Aceptable"
            elif nni <= MuestreoEspacialPuntos.NNI_RANDOM_HI:
                pat, cal = "Aleatorio", "Óptimo"  # noqa: F841
            elif nni <= MuestreoEspacialPuntos.NNI_DISPERSED_OK:
                pat, cal = "Disperso Moderado", "Bueno"  # noqa: F841
            elif nni <= 1.5:
                pat, cal = "Disperso Fuerte", "Óptimo" if _is_disp else "Aceptable"  # noqa: F841
            elif nni <= 2.8:
                pat, cal = "Muy Disperso", "Óptimo" if _is_disp else "Precaución ⚠"  # noqa: F841
            else:
                pat, cal = "Muy Disperso (área fragmentada) — Ver nota ⓘ", "Ver nota ⓘ"  # noqa: F841
            nni_txt = f"{nni:.4f}".replace('.', ',') if nni is not None else "NULL"
            rej_txt = f"{rej:.2f}".replace('.', ',')
            polys = r.get('polygons_hit_count', 'N/A')
            polys_pct = f" ({
                polys
                / self.total_polygon_count  # noqa: W503
                * 100:.1f}%)".replace(  # noqa: W503
                '.',
                ',') if self.total_polygon_count > 0 and isinstance(
                polys,
                int) else ""
            polys_cell = f"{polys}{polys_pct}" if polys != 'N/A' else 'N/A'
            contig = r.get('contiguous_count', 'N/A')
            contig_pct = f" ({
                contig
                / o  # noqa: W503
                * 100:.1f} %)".replace(  # noqa: W503
                '.',
                ',') if isinstance(
                contig,
                int) and isinstance(
                o,
                int) and o > 0 else ""
            contig_cell = f"{contig}{contig_pct}" if contig != 'N/A' else 'N/A'
            name = r.get('muestra')
            name_cell = f"<strong>{name}</strong>" if name in _top3_names else name
            nni_cell = f'<td title="IVMC = {nni_txt}">{pat}</td>' if nni is not None else '<td>N/A</td>'
            rows += f"<tr><td>{name_cell}</td><td>{s}</td><td>{o}</td><td>{rej_txt}%</td><td>{polys_cell}</td><td>{contig_cell}</td>{nni_cell}</tr>"  

        max_cont = max([r.get('contiguous_count', 0) or 0 for r in self.res]) if self.res else 0
        max_poly = max([r.get('polygons_hit_count', 0) or 0 for r in self.res]) if self.res else 0
        js_labels = [r.get('muestra', '') for r in self.res]
        js_nni = [r.get('indice_nn') or 'null' for r in self.res]
        js_eff = [r.get('n_puntos', 0) for r in self.res]
        js_poly = [r.get('polygons_hit_count') or 'null' for r in self.res]
        js_cont = [r.get('contiguous_count') or 'null' for r in self.res]
        js_rej = [((self.params.sample_size - r.get('n_puntos', 0)) / self.params.sample_size)
                  * 100 if self.params.sample_size > 0 else 0 for r in self.res]  # noqa: W503

        def fc(v): return f"{v:.2f}".replace('.', ',') if v is not None else "0,00"
        def fi(v): return f"{v:,}".replace(',', ' ') if v is not None else "0"

        hg_val = fc(self.eval.get('coverage_index')) if self.params.manual_grid_spacing > 0 else None
        hg_alert = ""
        if hg_val:
            hg_f = self.eval.get('coverage_index', 1.0)
            if hg_f < 0.7:
                hg_alert = f"<span style='color:#c62828; font-weight:bold;'>⚠ h/g = {hg_val} — g sobreestimado (malla más densa que g indicado).</span>"  
            elif hg_f > 1.5:
                hg_alert = f"<span style='color:#c62828; font-weight:bold;'>⚠ h/g = {hg_val} — g subestimado (malla menos densa que g indicado).</span>"  
            else:
                hg_alert = f"<span style='color:#2e7d32; font-weight:bold;'>✓ h/g = {hg_val} — coherente con el espaciado real.</span>"  

        dup_count = self.eval.get('duplicate_count', 0) if self.eval else 0
        dup_alert = f"<span style='color:#c62828; font-weight:bold;'>⚠ {dup_count} puntos duplicados detectados — se recomienda limpiar la capa de entrada.</span>" if dup_count and dup_count > 0 else "<span style='color:#2e7d32;'>✓ Sin puntos duplicados.</span>"  

        coll_pct = self.eval.get('collision_ratio', 0.0) if self.eval else 0.0
        coll_alert = f"<span style='color:#c62828; font-weight:bold;'>⚠ Colisiones: {
            fc(coll_pct)}% — supera el 10%: aumente el Orden Hilbert.</span>" if coll_pct >= 10.0 else (
            f"<span style='color:#2e7d32;'>✓ Colisiones: {
                fc(coll_pct)}% — Orden Hilbert adecuado.</span>{
                '<br><small style=\"color:#555;\">Marco detectado: <strong>Regular fragmentada</strong> — malla regular sobre área multipolígono. El índice R elevado refleja los saltos entre polígonos, no irregularidad de la malla.</small>' if _malla_frag else ''}")  

        return f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Reporte — Muestreo Espacial de Puntos</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/chartjs-plugin-annotation/2.2.1/chartjs-plugin-annotation.min.js"></script>
    <style>
        *, *::before, *::after {{ box-sizing: border-box; }}
        body {{ font-family: Arial, sans-serif; margin: 0; background: #f0f2f5; color: #222; line-height: 1.6; font-size: 15px; }}
        .container {{ max-width: 1280px; margin: 0 auto; padding: 24px 32px; }}
        .section-header {{ display: flex; align-items: center; gap: 12px; background: #1a3a5c; color: white; padding: 14px 22px; border-radius: 8px 8px 0 0; margin-top: 2.2em; margin-bottom: 0; }}
        .section-header .badge {{ background: #f0a500; color: #1a3a5c; font-weight: 800; font-size: 0.85em; padding: 2px 10px; border-radius: 20px; white-space: nowrap; }}
        .section-header h2 {{ margin: 0; font-size: 1.15em; color: white; border: none; }}
        .section-body {{ background: white; border: 1px solid #d0d7e3; border-top: none; border-radius: 0 0 8px 8px; padding: 24px 28px; margin-bottom: 0.5em; }}
        .card {{ background: white; border-radius: 8px; border: 1px solid #d0d7e3; padding: 20px 24px; margin-bottom: 16px; box-shadow: 0 1px 4px rgba(0,0,0,0.07); }}
        .card-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; }}
        .rec-card {{ border-radius: 10px; padding: 20px 22px; border: 2px solid; position: relative; }}
        .rec-card.rank-1 {{ border-color: #f0a500; background: #fffdf0; }}
        .rec-card.rank-2 {{ border-color: #78909c; background: #f7f9fa; }}
        .rec-card.rank-3 {{ border-color: #8d6e63; background: #fdf6f3; }}
        .rec-badge {{ position: absolute; top: -14px; left: 18px; font-weight: 800; font-size: 0.8em; padding: 2px 12px; border-radius: 20px; color: white; }}
        .rec-card.rank-1 .rec-badge {{ background: #f0a500; }}
        .rec-card.rank-2 .rec-badge {{ background: #78909c; }}
        .rec-card.rank-3 .rec-badge {{ background: #8d6e63; }}
        .rec-name {{ font-size: 1.05em; font-weight: 700; margin: 8px 0 14px 0; color: #1a3a5c; }}
        .metric-row {{ display: flex; justify-content: space-between; align-items: center; padding: 5px 0; border-bottom: 1px solid #eee; font-size: 0.92em; }}
        .metric-row:last-child {{ border-bottom: none; }}
        .metric-label {{ color: #555; }}
        .metric-val {{ font-weight: 700; }}
        .metric-val.good {{ color: #2e7d32; }} .metric-val.warn {{ color: #e65100; }} .metric-val.bad  {{ color: #c62828; }}
        .criteria-note {{ background: #e8eaf6; border-left: 4px solid #3949ab; padding: 10px 16px; border-radius: 4px; font-size: 0.88em; margin-bottom: 18px; color: #333; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 1.2em; font-size: 0.93em; }}
        th, td {{ text-align: left; padding: 10px 13px; border-bottom: 1px solid #e0e0e0; }}
        th {{ background: #1a3a5c; color: white; font-weight: 600; }}
        tr:nth-child(even) {{ background: #f8f9fb; }} tr:hover {{ background: #eef3fb; }}
        .config-table td {{ padding: 8px 12px; border: 1px solid #ddd; vertical-align: top; }}
        .config-table td:first-child, .config-table td:nth-child(3) {{ color: #555; font-size: 0.88em; white-space: nowrap; }}
        .chart-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 8px 0 16px 0; }}
        .chart-wrap {{ background: white; border-radius: 8px; border: 1px solid #d0d7e3; padding: 16px; position: relative; height: 280px; }}
        .chart-wrap canvas {{ max-height: 240px; }}
        .alert {{ padding: 10px 16px; border-radius: 6px; margin: 8px 0; font-size: 0.91em; }}
        .alert-green  {{ background: #e8f5e9; border-left: 4px solid #2e7d32; }}
        .alert-orange {{ background: #fff3e0; border-left: 4px solid #e65100; }}
        .alert-blue   {{ background: #e8eaf6; border-left: 4px solid #3949ab; }}
        details {{ border: 1px solid #d0d7e3; border-radius: 6px; margin: 10px 0; background: white; }}
        details > summary {{ padding: 12px 18px; font-weight: 600; cursor: pointer; color: #1a3a5c; list-style: none; user-select: none; }}
        details > summary::before {{ content: "▶ "; font-size: 0.8em; }} details[open] > summary::before {{ content: "▼ "; }}
        details > .details-body {{ padding: 16px 22px 20px 22px; border-top: 1px solid #e0e0e0; }}
        .appendix-label {{ display: inline-block; background: #78909c; color: white; font-size: 0.75em; font-weight: 700; padding: 1px 8px; border-radius: 10px; margin-right: 8px; vertical-align: middle; }}
        .report-title {{ background: #1a3a5c; color: white; border-radius: 8px; padding: 20px 28px; margin-bottom: 4px; }}
        .report-title h1 {{ margin: 0 0 6px 0; font-size: 1.4em; color: white; border: none; }}
        .report-title p  {{ margin: 0; font-size: 0.9em; opacity: 0.85; }}
        .grid-badge {{ display: inline-block; padding: 2px 12px; border-radius: 14px; font-size: 0.85em; font-weight: 700; }}
        .grid-badge.regular  {{ background: #e8f5e9; color: #2e7d32; border: 1px solid #a5d6a7; }}
        .grid-badge.irregular {{ background: #e3f2fd; color: #1565c0; border: 1px solid #90caf9; }}
        .interpretation-table th {{ background: #4CAF50; }} .interpretation-table {{ width:100%; margin:1em 0; }}
        pre {{ background:#f4f4f9; padding:10px; border-radius:4px; font-size:0.88em; overflow-x:auto; }}
        h3 {{ color: #1a3a5c; border-bottom: 1px solid #e0e0e0; padding-bottom: 6px; margin-top: 1.4em; }}
        h4 {{ color: #2e5480; margin-top: 1.2em; }}
        .footer {{ margin-top: 3em; padding-top: 1.2em; border-top: 1px solid #ddd; font-size: 0.85em; color: #777; text-align: center; }}
    </style>
</head>
<body>
<div class="container">
    <div class="report-title">
        <h1>Reporte de Muestreo Espacial de Puntos</h1>
        <p><strong>Proyecto:</strong> {self.proj} &nbsp;|&nbsp; <strong>Generado:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M")} &nbsp;|&nbsp; <strong>Tiempo de ejecución:</strong> {self._fmt_elapsed(getattr(self, "_elapsed_at_write", 0.0))} &nbsp;|&nbsp; <strong>Versión:</strong> {MuestreoEspacialPuntos.VERSION}</p>
    </div>

    <div class="section-header"><span class="badge">SECCIÓN 1</span><h2>Resumen de Ejecución</h2></div>
    <div class="section-body">
        <h3>Configuración</h3>
        <table class="config-table">
            <tbody>
            <tr>
                <td><strong>Método de muestreo:</strong></td>
                <td>{self.meth.get(self.params.method)}</td>
                <td><strong>Tipo de marco muestral detectado:</strong></td>
                <td><span class="grid-badge {'regular' if _any_regular else 'irregular'}">{grid_label}</span></td>
            </tr>
            <tr>
                <td><strong>Tamaño de muestra (n):</strong></td>
                <td>{self.params.sample_size} {str("(Automático JSON)" if self.params.method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON) else "")}</td>
                <td><strong>Puntos en el marco (inicial):</strong></td>
                <td>{fi(self.init)}</td>
            </tr>
            <tr>
                <td><strong>Intensidad de muestreo:</strong></td>
                <td><strong>{f"{(self.params.sample_size / self.filt * 100):.2f}".replace('.', ',') if self.filt > 0 else "N/A"} %</strong> <span style="color:#555; font-size:0.88em;">&nbsp;(n / N filtrado)</span></td>
                <td><strong>Puntos en el marco (filtrado):</strong></td>
                <td>{fi(self.filt)}</td>
            </tr>
            <tr>
                <td><strong>Iteraciones:</strong></td>
                <td>{self.params.num_iterations}</td>
                <td><strong>Distancia mínima entre puntos:</strong></td>
                <td>{self.params.min_sample_distance if self.params.min_sample_distance > 0 else "No aplica"}</td>
            </tr>
            <tr>
                <td><strong>Distancia al perímetro:</strong></td>
                <td>{self.params.perimeter_distance if self.params.perimeter_distance > 0 else "No aplica"}</td>
                <td><strong>Orden Hilbert &nbsp;/&nbsp; Grupos (k):</strong></td>
                <td>{self.params.hilbert_order} &nbsp;&nbsp;/&nbsp;&nbsp; {self.params.num_groups if self.params.num_groups > 0 else f'Auto ({self.res[0].get("k_efectivo", "N/A") if self.res else "N/A"})'}</td>
            </tr>
            <tr>
                <td><strong>Corrección de borde (Retracción):</strong></td>
                <td>{'Sí — radio ' + str(self.params.plot_radius) if self.params.apply_pullback else 'No aplica'}</td>
                <td><strong>Área de estudio:</strong></td>
                <td>{fc(self.eval.get('area_km2'))} km² &nbsp;|&nbsp; Compacidad: {f"{self.eval.get('indice_compacidad', 0):.3f}".replace('.', ',')} ({self._get_compacidad_desc(self.eval.get('indice_compacidad', 0))}) &nbsp;|&nbsp; <strong>Polígonos: {fi(self.total_polygon_count)}</strong> {'<span style="color:#555; font-size:0.88em;"> (área multipolígono)</span>' if self.total_polygon_count > 1 else ''}</td>  
            </tr>
            <tr>
                <td><strong>Distancia de malla digitada (g):</strong></td>
                <td>{'<strong>' + fc(self.params.manual_grid_spacing) + ' m</strong>' if self.params.manual_grid_spacing > 0 else '<span style="color:#999;">No ingresada</span>'}</td>  
                <td><strong>Distancia de malla teórica (h = √A/N):</strong></td>
                <td><strong>{fc(self.eval.get('mesh_distance_h'))} m</strong> &nbsp;—&nbsp; {'<span style="color:#555;">Índice de Cobertura h/g = <strong>' + fc(self.eval.get('coverage_index')) + '</strong> &nbsp;(' + ('g sobreestimado ⚠' if self.eval.get('coverage_index', 1) < 0.7 else ('g subestimado ⚠' if self.eval.get('coverage_index', 1) > 1.5 else 'coherente ✓')) + ')</span>' if self.params.manual_grid_spacing > 0 else '<span style="color:#999;">h/g no calculable — ingrese g para comparar</span>'}</td>  
            </tr>
            </tbody>
        </table>

        <h3>Diagnóstico de Pre-ejecución</h3>
        <p style="font-size:0.88em; color:#555;">Estos indicadores validan la configuración del algoritmo. Valores fuera de rango sugieren ajustes antes de usar los resultados.</p>  
        <div class="card-grid">
            <div class="card"><strong>Índice de Cobertura (h/g)</strong><br><small style="color:#555;">Coherencia entre el espaciado g configurado y el espaciado real h de la malla.</small><br><br>{hg_alert if hg_alert else "<span style='color:#777;'>No aplica — no se ingresó el parámetro g.</span>"}</div>  
            <div class="card"><strong>Puntos Duplicados</strong><br><small style="color:#555;">Puntos con coordenadas idénticas en la capa de entrada.</small><br><br>{dup_alert}</div>  
            <div class="card"><strong>Colisiones Hilbert ({fc(self.eval.get('collision_ratio'))}%)</strong><br><small style="color:#555;">Fracción de puntos que comparten la misma celda de Hilbert.</small><br><br>{coll_alert}</div>  
        </div>
    </div>

    <div class="section-header"><span class="badge">SECCIÓN 2</span><h2>Panel de Selección de Muestras</h2></div>
    <div class="section-body">
        <h3>★ Muestras Recomendadas</h3>
        <div class="criteria-note">
            <strong>Criterios de selección aplicados (en orden de prioridad):</strong>
            {'<ol style="margin:6px 0 0 0; padding-left:20px;"><li><strong>Polígonos muestreados</strong> — maximizar la cobertura de sub-polígonos del área de estudio (criterio eliminatorio).</li><li><strong>% Rechazo</strong> — minimizar la pérdida de puntos por distancia mínima.</li><li><strong>Puntos contiguos</strong> — minimizar pares de puntos adyacentes en la malla.</li><li><strong>IVMC</strong> — maximizar la dispersión espacial global.</li><li><strong>Evaluación visual</strong> — verificar en QGIS que los puntos cubren homogéneamente el polígono.</li></ol>' if _is_disp else '<ol style="margin:6px 0 0 0; padding-left:20px;"><li><strong>IVMC en banda óptima (0,95–1,05)</strong> — filtro de aceptación para malla irregular.</li><li><strong>Polígonos muestreados</strong> — maximizar cobertura entre candidatos aceptados.</li><li><strong>% Rechazo</strong> — minimizar pérdida de puntos.</li><li><strong>Puntos contiguos</strong> — minimizar pares adyacentes.</li><li><strong>Evaluación visual</strong> — confirmar en QGIS que no existen zonas sin cobertura.</li></ol>'}  
            <p style="margin:8px 0 0 0; font-size:0.87em; color:#555;"><strong>Cómo hacer la evaluación visual en QGIS:</strong> active la capa de la muestra candidata junto al polígono del área de estudio y verifique que no existan sub-zonas del polígono sin ningún punto muestreado.</p>  
        </div>
        {self._get_best()}

        <h3>Visualización Comparativa</h3>
        <div class="chart-grid">
            <div class="chart-wrap"><canvas id="polygonsChart"></canvas></div>
            <div class="chart-wrap"><canvas id="rejectionChart"></canvas></div>
            <div class="chart-wrap"><canvas id="contiguousChart"></canvas></div>
            <div class="chart-wrap"><canvas id="nniChart"></canvas></div>
            <div class="chart-wrap"><canvas id="effectiveChart"></canvas></div>
        </div>

        <h3>Resultados por Iteración</h3>
        <table>
            <thead><tr><th>Muestra</th><th title="Número de puntos solicitados en la configuración">n solicitado</th><th title="Puntos efectivamente seleccionados tras aplicar filtros">n obtenido</th><th title="Menor es mejor — ↓ min">% Rechazo ↓</th><th title="Mayor es mejor — ↑ max">Polígonos ↑ (de {fi(self.total_polygon_count)})</th><th title="Menor es mejor — ↓ min">Contiguos ↓</th><th title="Índice de Vecino Más Cercano Corregido">Patrón espacial (IVMC)</th></tr></thead>  
            <tbody>{rows}</tbody>
        </table>
        <details style="margin-top:-0.8em; margin-bottom:1.2em;">
            <summary style="font-size:0.88em; color:#2e5480; padding:8px 4px;">ⓘ Explicación de cada columna</summary>
            <div style="padding:12px 16px; background:#f8f9fb; border:1px solid #e0e0e0; border-top:none; border-radius:0 0 6px 6px; font-size:0.87em; line-height:1.7;">  
                <p><strong>n solicitado:</strong> tamaño de muestra configurado antes de la ejecución.</p>
                <p><strong>n obtenido:</strong> puntos efectivamente incluidos en la muestra tras aplicar todos los filtros.</p>  
                <p><strong>% Rechazo:</strong> porcentaje de puntos descartados = (n solicitado − n obtenido) / n solicitado × 100.</p>  
                <p><strong>Polígonos muestreados:</strong> número de sub-polígonos del área de estudio que contienen al menos un punto.</p>  
                <p><strong>Contiguos:</strong> número de puntos de la muestra cuyo vecino más cercano está a una distancia ≤ g.</p>  
                <p><strong>Patrón espacial (IVMC):</strong> descripción textual del patrón detectado. El valor numérico se muestra al pasar el cursor.</p>  
                <p style="margin-top:10px; padding:8px 12px; background:#e8eaf6; border-left:3px solid #3949ab; border-radius:3px;"><strong>Regla de decisión rápida:</strong> seleccione la muestra con mayor número de polígonos muestreados → entre ellas, la de menor % rechazo → entre ellas, la de menor número de contiguos → el IVMC solo decide si todo lo anterior está empatado.</p>  
            </div>
        </details>
    </div>

    <div class="section-header" style="margin-top:2.2em;"><span class="badge">SECCIÓN 3</span><h2>Apéndices Técnicos</h2></div>  
    <div class="section-body">
        <p style="font-size:0.9em; color:#555; margin-top:0;">Información de referencia para el administrador del algoritmo. No es necesaria para seleccionar una muestra.</p>  
        <details>
            <summary><span class="appendix-label">A1</span> Marco de Validación Cuantitativo de la Curva de Hilbert</summary>  
            <div class="details-body">
                <div class="alert alert-blue"><strong>Nota:</strong> Las métricas R, σ/h, L/Nh y CV son descriptores de la geometría del área de estudio, no indicadores de calidad.</div>  
                <table class="config-table">
                    <tbody>
                    <tr><th colspan="4" style="text-align:center; background:#4CAF50;">Métricas Normalizadas y Diagnóstico</th></tr>  
                    <tr><th>Métrica</th><th>Rango Excelente</th><th>Valor</th><th>Diagnóstico</th></tr>
                    <tr style="background:{self._get_color(self.eval.get('efficiency_index_R', 0), 0.8, 1.5)}"><td>Índice de eficiencia (R = d/h)</td><td>0,8 – 1,5</td><td>{fc(self.eval.get('efficiency_index_R'))}</td><td>{self._get_diag(self.eval.get('efficiency_index_R', 0), 0.8, 1.5)}</td></tr>  
                    <tr style="background:{self._get_color(self.eval.get('std_dev_norm', 0), 0.5, 1.0)}"><td>Desv. estándar normalizada (σ/h)</td><td>0,5 – 1,0</td><td>{fc(self.eval.get('std_dev_norm'))}</td><td>{self._get_diag(self.eval.get('std_dev_norm', 0), 0.5, 1.0)}</td></tr>  
                    <tr style="background:{self._get_color(self.eval.get('length_norm', 0), 0.9, 1.6)}"><td>Longitud total normalizada (L/Nh)</td><td>0,9 – 1,6</td><td>{fc(self.eval.get('length_norm'))}</td><td>{self._get_diag(self.eval.get('length_norm', 0), 0.9, 1.6)}</td></tr>  
                    <tr><td>Coeficiente de variación (CV)</td><td>&lt; 0,5</td><td>{fc(self.eval.get('cv'))}</td><td>{self._get_cv_desc(self.eval.get('cv', 0))}</td></tr>  
                    </tbody>
                </table>
            </div>
        </details>

        <details>
            <summary><span class="appendix-label">A2</span> Guía de Interpretación del IVMC</summary>
            <div class="details-body">
                <p>El <strong>IVMC Corregido</strong> mide el patrón espacial ajustado por la forma del polígono:</p>
                <pre>IVMC_corr = (d_obs / d_esp) × √(A_bbox / A_real)</pre>
                <div class="alert alert-orange"><strong>Malla regular compacta:</strong> IVMC tiene muy baja varianza. Usar como criterio de desempate únicamente.</div>  
                <div class="alert alert-orange" style="background:#fff8e1; border-color:#f57f17; margin-top:8px;"><strong>Malla regular fragmentada:</strong> valores altos (4–7) son normales por huecos. Use ratio IVMC_método / IVMC_AL para comparar.</div>  
                <h4>Malla irregular</h4>
                <table class="interpretation-table">
                    <thead><tr><th>Rango IVMC</th><th>Patrón</th><th>AL / KM (Obj: Aleatoriedad)</th><th>SH / GH (Obj: Dispersión)</th></tr></thead>  
                    <tbody>
                        <tr style="background:#fdecea"><td>IVMC &lt; 0,8</td><td>Agrupado Fuerte</td><td>No Aceptable</td><td>No Aceptable</td></tr>  
                        <tr style="background:#fdecea"><td>0,8 ≤ IVMC &lt; 0,95</td><td>Agrupado Moderado</td><td>No Aceptable</td><td>No Aceptable</td></tr>  
                        <tr style="background:#e8f5e9"><td>0,95 ≤ IVMC ≤ 1,05</td><td>Aleatorio</td><td><strong>Óptimo</strong></td><td>Bajo el objetivo</td></tr>  
                        <tr style="background:#e3f2fd"><td>1,05 &lt; IVMC ≤ 1,2</td><td>Disperso Moderado</td><td>Bueno</td><td>Bueno</td></tr>  
                        <tr style="background:#e8f5e9"><td>1,2 &lt; IVMC ≤ 1,5</td><td>Disperso Fuerte</td><td>Aceptable</td><td><strong>Óptimo</strong></td></tr>  
                        <tr style="background:#fff8e1"><td>IVMC &gt; 1,5</td><td>Muy Disperso</td><td>Precaución ⚠</td><td>Precaución ⚠</td></tr>  
                    </tbody>
                </table>
            </div>
        </details>

        <details>
            <summary><span class="appendix-label">A3</span> Metodología de Cálculo del IVMC Corregido</summary>
            <div class="details-body">
                <p>Cálculo nativo en Python usando <code>QgsSpatialIndex</code> (O(n log n)).</p>
                <pre>IVMC_corr = d_obs × 2 × √(n × A_bbox / A_real²)</pre>
            </div>
        </details>

        {self._get_proportionality_section()}
        {self._get_json_strata_section()}
    </div>
    <div class="footer"><p>Complemento QGIS: <strong>Muestreo Espacial de Puntos</strong> &nbsp;|&nbsp; Versión {MuestreoEspacialPuntos.VERSION} &nbsp;|&nbsp; Jorge Fallas</p></div>  
</div>
<script>
    const labels = {js_labels};
    const maxPoly = {max_poly}; const totalPolys = {self.total_polygon_count};
    const polyData = {js_poly}; const polyCovPct = polyData.map(v => (v !== null && v !== 'null' && totalPolys > 0) ? +((v / totalPolys) * 100).toFixed(2) : null);  
    new Chart(document.getElementById('polygonsChart'), {{ type: 'bar', data: {{ labels: labels, datasets: [{{ label: 'Polígonos muestreados', data: polyData, backgroundColor: 'rgba(46,125,50,0.65)', yAxisID: 'y' }}, {{ label: '% Cobertura', data: polyCovPct, type: 'line', borderColor: 'rgba(183,28,28,0.85)', borderWidth: 2, fill: false, yAxisID: 'y2' }}] }}, options: {{ plugins: {{ title: {{ display: true, text: '① Polígonos Muestreados ({fi(self.total_polygon_count)} total)' }} }}, scales: {{ y: {{ type: 'linear', position: 'left', min: 0, suggestedMax: totalPolys > 0 ? totalPolys : 5 }}, y2: {{ type: 'linear', position: 'right', min: 0, max: 100, grid: {{ drawOnChartArea: false }} }} }} }} }});  
    const rejVals = {js_rej}.filter(v => v > 0);
    new Chart(document.getElementById('rejectionChart'), {{ type: 'bar', data: {{ labels, datasets: [{{ label: '% Rechazo', data: {js_rej}, backgroundColor: 'rgba(198,40,40,0.60)' }}] }}, options: {{ plugins: {{ title: {{ display: true, text: '② % Rechazo (↓ menor es mejor)' }} }}, scales: {{ y: {{ min: 0, suggestedMax: rejVals.length ? Math.ceil(Math.max(...rejVals)*1.2) : 10 }} }} }} }});  
    const maxCont = {max_cont};
    new Chart(document.getElementById('contiguousChart'), {{ type: 'bar', data: {{ labels, datasets: [{{ label: 'Puntos Contiguos', data: {js_cont}, backgroundColor: 'rgba(153,102,255,0.65)' }}] }}, options: {{ plugins: {{ title: {{ display: true, text: '③ Puntos Contiguos (↓ menor es mejor)' }} }}, scales: {{ y: {{ beginAtZero: true, suggestedMax: maxCont > 0 ? Math.ceil(maxCont*1.2) : 5 }} }} }} }});  
    const nniVals = {js_nni}.filter(v => v !== null && v !== 'null');
    new Chart(document.getElementById('nniChart'), {{ type: 'bar', data: {{ labels, datasets: [{{ label: 'IVMC', data: {js_nni}, backgroundColor: 'rgba(54,162,235,0.60)' }}] }}, options: {{ plugins: {{ title: {{ display: true, text: '④ IVMC (↑ desempate final)' }} }}, scales: {{ y: {{ min: 0, suggestedMax: nniVals.length ? +(Math.max(...nniVals)*1.15).toFixed(2) : 3.0 }} }} }} }});  
    const effVals = {js_eff}.filter(v => v > 0);
    new Chart(document.getElementById('effectiveChart'), {{ type: 'bar', data: {{ labels, datasets: [{{ label: 'Puntos Obtenidos', data: {js_eff}, backgroundColor: 'rgba(75,192,192,0.65)' }}] }}, options: {{ plugins: {{ title: {{ display: true, text: 'Puntos Obtenidos' }} }}, scales: {{ y: {{ min: effVals.length ? Math.floor(Math.min(...effVals)*0.95) : 0, suggestedMax: effVals.length ? Math.ceil(Math.max(...effVals)*1.05) : {self.params.sample_size} }} }} }} }});  
</script>
</body></html>"""

    def _prop_quality_label(self, max_diff: float) -> Tuple[str, str]:
        if max_diff <= 1.0:
            return "Excelente", "#e8f5e9"
        if max_diff <= 3.0:
            return "Adecuado", "#fff9c4"
        if max_diff <= 6.0:
            return "Revisión ⚠", "#fff3e0"
        return "Deficiente ⛔", "#ffebee"

    def _render_proportionality_block(self, it_label: str, prop: Dict, method_note: str = "") -> str:
        max_diff = prop.get('max_diff', 0.0)
        mean_diff = prop.get('mean_diff', 0.0)
        zeros = prop.get('groups_zero', [])
        rows_data = prop.get('rows', [])
        label, bg = self._prop_quality_label(max_diff)
        def fp(v): return f"{v:.4f}".replace('.', ',')
        def fi(v): return f"{v:,}".replace(',', ' ')
        MAX_ROWS_HTML = 500
        detail_rows = ""
        for r in rows_data[:MAX_ROWS_HTML]:
            diff_v = r['diff']
            if abs(diff_v) <= 1.0:
                diff_color = ""
            elif diff_v > 0:
                diff_color = ' style="color:#c62828; font-weight:bold;"'
            else:
                diff_color = ' style="color:#1565c0; font-weight:bold;"'
            detail_rows += f"<tr><td>{r['group_id']}</td><td>{fi(r['pop_n'])}</td><td>{fi(r['sam_n'])}</td><td>{fp(r['pct_pop'])} %</td><td>{fp(r['pct_sam'])} %</td><td{diff_color}>{fp(r['diff'])} %</td></tr>\n"  
        truncation_note = f"<p><small>⚠ Tabla truncada: {MAX_ROWS_HTML} de {
            len(rows_data)} grupos.</small></p>" if len(rows_data) > MAX_ROWS_HTML else ""
        zeros_note = f"<p style='color:#b71c1c;'><strong>⚠ Grupos sin representación:</strong> {
            len(zeros)} grupo(s).</p>" if zeros else ""
        method_note_html = f"<p><small>{method_note}</small></p>" if method_note else ""
        return f"<h3>{it_label}</h3>{method_note_html}<table class='config-table'><tbody><tr><td><strong>Máxima Diferencia (|Δ%| máx.):</strong></td><td style='background-color:{bg}; font-weight:bold;'>{  
            fp(max_diff)} % — {label}</td><td><strong>Diferencia Media (|Δ%| prom.):</strong></td><td>{
            fp(mean_diff)} %</td></tr><tr><td><strong>Total grupos:</strong></td><td>{
            fi(
                len(rows_data))}</td><td><strong>Grupos con 0 puntos:</strong></td><td>{
                    len(zeros)}</td></tr></tbody></table>{zeros_note}<details><summary style='cursor:pointer; color:#0056b3; font-weight:bold;'>▶ Ver tabla de detalle por grupo ({  
                        len(rows_data)} filas)</summary>{truncation_note}<table><thead><tr><th>ID Grupo</th><th>N Grupo</th><th>n Asignado</th><th>% Población</th><th>% Muestra</th><th>Diferencia (Δ%)</th></tr></thead><tbody>{detail_rows}</tbody></table></details><hr style='margin:1.5em 0; border:none; border-top:1px solid #ddd;'>"  

    def _get_json_strata_section(self) -> str:
        if self.params.method not in (SamplingMethod.STRATIFIED_BY_FIELD,
                                      SamplingMethod.STRATIFIED_BY_POLYGON) or not self.json_strata_val:
            return ""
        blocks = ""
        for vd in self.json_strata_val:
            it, pop, req, sam = vd['iteration'], vd['population'], vd['requested'], vd['sampled']
            rows = ""
            for strat_val in sorted(pop.keys()):
                p, r, s = pop[strat_val], req.get(strat_val, 0), sam.get(strat_val, 0)
                deficit = r - s
                def_txt = f"<span style='color:#c62828; font-weight:bold;'>{deficit}</span>" if deficit > 0 else "0"
                rows += f"<tr><td>{strat_val}</td><td>{p}</td><td>{r}</td><td>{s}</td><td>{def_txt}</td></tr>"
            blocks += f"<h4>Iteración {it}</h4><table><thead><tr><th>Estrato</th><th>Puntos Disponibles</th><th>Puntos Solicitados (JSON)</th><th>Puntos Obtenidos (Al Azar)</th><th>Déficit (Faltantes)</th></tr></thead><tbody>{rows}</tbody></table>"  
        return f"<details open><summary><span class='appendix-label'>A5</span> Desempeño del Muestreo Estratificado Aleatorio (JSON)</summary><div class='details-body'><p>Verificación de la cuota solicitada por estrato vs la cuota obtenida de manera <b>Aleatoria</b>. Un déficit mayor a 0 significa que no había suficientes puntos disponibles en el área o que fueron rechazados por el filtro de Distancia Mínima.</p>{blocks}</div></details>"  

    def _get_proportionality_section(self) -> str:
        is_gh = self.params.method in (
            SamplingMethod.STRATIFIED_HILBERT,
            SamplingMethod.GROUPS_ROW_COL) and self.hilbert_group_val
        is_km = self.params.method == SamplingMethod.KMEANS_GROUPS and self.k_val
        if not is_gh and not is_km:
            return ""
        _is_gf = self.params.method == SamplingMethod.GROUPS_ROW_COL
        NOTE_GH = "Los grupos siguen el orden fila NO→SE." if _is_gf else "Los estratos Hilbert son de tamaño aproximadamente igual por diseño."  
        NOTE_KM = "La asignación es proporcional al tamaño real del clúster K-Medias."
        blocks = ""
        if is_gh:
            for vd in self.hilbert_group_val:
                label = f"Iteración {
                    vd.iteration} — Grupos Fila-Col (NO→SE)" if _is_gf else f"Iteración {
                    vd.iteration} — Grupos Hilbert"
                blocks += self._render_proportionality_block(label, vd.proportionality, NOTE_GH)
        if is_km:
            for vd in self.k_val:
                if vd.proportionality:
                    blocks += self._render_proportionality_block(
                        f"Iteración {vd.iteration} — K-Medias", vd.proportionality, NOTE_KM)
        return f"<h2>Validación de Proporcionalidad por Grupo</h2><details><summary><span class='appendix-label'>A4</span> Desempeño de Proporcionalidad</summary><div class='details-body'>{blocks}</div></details>"  

    def _get_compacidad_desc(self, v):
        if v > 0.85:
            return "Compacto"
        if v > 0.6:
            return "Moderadamente Compacto"
        if v > 0.4:
            return "Alargado o Irregular"
        return "Muy Alargado o Fragmentado"

    def _get_diag(self, v, min_v, max_v):
        if min_v <= v <= max_v:
            return "Excelente"
        rng = max_v - min_v
        if (min_v - rng * 0.5) <= v <= (max_v + rng * 0.5):
            return "Adecuado"
        if (min_v - rng) <= v <= (max_v + rng):
            return "Revisión"
        return "Deficiente"

    def _get_color(self, v, min_v, max_v):
        if min_v <= v <= max_v:
            return "#e8f5e9"
        rng = max_v - min_v
        if (min_v - rng * 0.5) <= v <= (max_v + rng * 0.5):
            return "#fff9c4"
        if (min_v - rng) <= v <= (max_v + rng):
            return "#fff3e0"
        return "#ffebee"

    def _get_cv_desc(self, v):
        if v < 0.5:
            return "Uniforme / Regular"
        if v < 0.8:
            return "Relativamente Homogéneo"
        if v < 1.2:
            return "Moderadamente Agrupado"
        return "Agrupado / Fragmentado"

    def _get_best(self):
        if not self.res:
            return "<p>No hay datos de iteraciones.</p>"
        valid = [r for r in self.res if r.get('indice_nn') is not None]
        if not valid:
            return "<p>Sin IVMC válido en ninguna iteración.</p>"

        import math as _math
        s_size = self.params.sample_size
        _n_medio = sum(r.get('n_puntos', 0) for r in valid) / len(valid)
        _alto_n = _n_medio > 5000
        _r_e = self.eval.get('efficiency_index_R', 0.0) if self.eval else 0.0
        _co_e = self.eval.get('collision_ratio', 100.0) if self.eval else 100.0
        _malla_r = (0.85 <= _r_e <= 1.35 and _co_e < 10.0)
        _malla_f = (not _malla_r and _co_e < 10.0 and _r_e > 1.35)
        is_disp = (_malla_r or _malla_f or self.params.min_sample_distance > 0 or self.params.method in (
            SamplingMethod.SYSTEMATIC_HILBERT, SamplingMethod.STRATIFIED_HILBERT, SamplingMethod.GROUPS_ROW_COL))
        if self.params.method in (SamplingMethod.STRATIFIED_BY_FIELD, SamplingMethod.STRATIFIED_BY_POLYGON):
            is_disp = False

        def rej_rate(r): return ((s_size - r.get('n_puntos', 0)) / s_size * 100) if s_size > 0 else 0.0
        def ivnorm(r): return (r.get('indice_nn') or 0) / _math.sqrt(max(r.get('n_puntos') or 1, 1))

        if _alto_n:
            pvars = len(set(r.get('polygons_hit_count', 0) for r in valid)) > 1
            if pvars:
                top3 = sorted(valid, key=lambda x: (x.get('polygons_hit_count') or 0, ivnorm(x)), reverse=True)[:3]
                context_note = f"n ≈ {
                    int(_n_medio):,                                      } pts. Criterio principal: cobertura de polígonos. Desempate: IVMC/√n."  
            else:
                top3 = sorted(valid, key=ivnorm, reverse=True)[:3]
                context_note = f"n ≈ {int(_n_medio):,} pts. Un solo polígono. Criterio: IVMC/√n."
        elif is_disp:
            top3 = sorted(valid, key=lambda x: (-(x.get('polygons_hit_count') or 0), rej_rate(x),
                          (x.get('contiguous_count') or 0), -(x.get('indice_nn') or 0)))[:3]
            context_note = "Objetivo: dispersión máxima."
        else:
            opt = [r for r in valid if MuestreoEspacialPuntos.NNI_RANDOM_LO <= (
                r.get('indice_nn') or 0) <= MuestreoEspacialPuntos.NNI_RANDOM_HI]
            if opt:
                top3 = sorted(opt, key=lambda x: (-(x.get('polygons_hit_count') or 0),
                              rej_rate(x), (x.get('contiguous_count') or 0)))[:3]
                context_note = "Objetivo: aleatoriedad. Filtro IVMC 0,95–1,05 aplicado."
            else:
                top3 = sorted(valid, key=lambda x: abs((x.get('indice_nn') or 99) - 1.0))[:3]
                context_note = "Objetivo: aleatoriedad. Ninguna en rango óptimo — más cercanas a IVMC = 1,0."

        rank_labels = ["🥇 1.ª opción", "🥈 2.ª opción", "🥉 3.ª opción"]
        rank_css = ["rank-1", "rank-2", "rank-3"]
        cards_html = ""

        for i, r in enumerate(top3):
            ivmc, n_pts, rej = r.get('indice_nn'), r.get('n_puntos', 0), rej_rate(r)
            polys, contig = r.get('polygons_hit_count', 'N/A'), r.get('contiguous_count', 'N/A')
            name = r.get('muestra', f'Iteración {i + 1}')
            ivmc_s = f"{ivmc:.4f}".replace('.', ',') if ivmc else "N/A"
            rej_s = f"{rej:.1f}".replace('.', ',')
            ivnrm_s = (f"{ivnorm(r):.4f}".replace('.', ',') if _alto_n and ivmc else "")

            def cv(val, good_thresh, bad_thresh, higher_is_better=True):
                try:
                    v = float(val)
                    if higher_is_better:
                        return "good" if v >= good_thresh else ("warn" if v >= bad_thresh else "bad")
                    else:
                        return "good" if v <= good_thresh else ("warn" if v <= bad_thresh else "bad")
                except Exception:
                    return ""

            rej_cls = cv(rej, 0.0, 5.0, False)
            poly_cls = cv(polys, 1, 1, True)
            cont_cls = cv(contig, 0, 3, False) if contig != 'N/A' else ""
            ivmc_cls = cv(
                ivmc,
                MuestreoEspacialPuntos.NNI_DISPERSED_OK,
                MuestreoEspacialPuntos.NNI_RANDOM_HI,
                True) if is_disp else cv(
                ivmc,
                MuestreoEspacialPuntos.NNI_RANDOM_HI,
                MuestreoEspacialPuntos.NNI_RANDOM_LO,
                False)
            extra_row = f'<div class="metric-row"><span class="metric-label">IVMC/√n</span><span class="metric-val">{ivnrm_s}</span></div>' if _alto_n and ivnrm_s else ""  
            contig_note = "" if self.params.manual_grid_spacing > 0 else " <small style='color:#999;'>(requiere g)</small>"  

            cards_html += f"""
            <div class="rec-card {rank_css[i]}">
                <span class="rec-badge">{rank_labels[i]}</span>
                <div class="rec-name">{name}</div>
                <div class="metric-row"><span class="metric-label">① Polígonos muestreados</span><span class="metric-val {poly_cls}">{polys}</span></div>  
                <div class="metric-row"><span class="metric-label">② % Rechazo</span><span class="metric-val {rej_cls}">{rej_s} %</span></div>  
                <div class="metric-row"><span class="metric-label">③ Puntos contiguos{contig_note}</span><span class="metric-val {cont_cls}">{contig}</span></div>  
                <div class="metric-row"><span class="metric-label">④ IVMC (desempate)</span><span class="metric-val {ivmc_cls}">{ivmc_s}</span></div>  
                {extra_row}
                <div class="metric-row" style="margin-top:6px;border-top:2px solid #ddd;padding-top:6px;"><span class="metric-label" style="color:#777;">Puntos obtenidos</span><span class="metric-val">{n_pts} / {s_size}</span></div>  
            </div>"""

        return f"<p style='font-size:0.88em;color:#555;margin-bottom:14px;'><strong>Contexto:</strong> {context_note}</p><div style='display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:24px;margin-top:16px;'>{cards_html}</div>"  
