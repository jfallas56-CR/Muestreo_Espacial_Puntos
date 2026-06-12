# -*- coding: utf-8 -*-
"""
Proveedor de Processing para el complemento Muestreo Espacial de Puntos.

Registra MuestreoEspacialPuntos en la caja de herramientas de QGIS Processing
y expone el ícono del complemento (icon.png) a nivel de proveedor y algoritmo.

Compatibilidad: QGIS 3.28 LTR / 3.44 LTR / 4.0  (Qt5 y Qt6)
Autor:          Jorge Fallas · jfallas56@gmail.com
Versión:        1.0.0
"""

import os

from qgis.PyQt.QtGui import QIcon
from qgis.core import QgsProcessingProvider


class MuestreoEspacialProvider(QgsProcessingProvider):
    """Proveedor Processing que agrupa los algoritmos del complemento."""

    # ID único del proveedor — no cambiar entre versiones
    PROVIDER_ID = 'muestreo_espacial_puntos'

    def __init__(self):
        super().__init__()
        self._plugin_dir = os.path.dirname(__file__)

    # ------------------------------------------------------------------
    # Identidad del proveedor
    # ------------------------------------------------------------------

    def id(self) -> str:
        return self.PROVIDER_ID

    def name(self) -> str:
        return 'Muestreo Espacial de Puntos'

    def longName(self) -> str:
        return 'Muestreo Espacial de Puntos'

    def icon(self) -> QIcon:
        """Ícono del proveedor en la caja de herramientas de Processing."""
        icon_path = os.path.join(self._plugin_dir, 'icon.png')
        if os.path.isfile(icon_path):
            return QIcon(icon_path)
        return super().icon()

    def svgIconPath(self) -> str:
        """Ruta al ícono (PNG aceptado también). Usado por algunos diálogos."""
        icon_path = os.path.join(self._plugin_dir, 'icon.png')
        return icon_path if os.path.isfile(icon_path) else ''

    # ------------------------------------------------------------------
    # Registro de algoritmos
    # ------------------------------------------------------------------

    def loadAlgorithms(self) -> None:
        """Registra todos los algoritmos del proveedor."""
        from .Muestreo_Espacial_Puntos import MuestreoEspacialPuntos
        self.addAlgorithm(MuestreoEspacialPuntos())
