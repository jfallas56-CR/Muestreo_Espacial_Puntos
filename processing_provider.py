# -*- coding: utf-8 -*-
"""
Proveedor de algoritmos de Processing para el complemento
Muestreo_Espacial_Puntos.
"""

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .Muestreo_Espacial_Puntos import MuestreoEspacialPuntos


class MuestreoEspacialProvider(QgsProcessingProvider):
    """Proveedor de algoritmos de Processing para Muestreo Espacial de Puntos."""

    def __init__(self):
        super().__init__()

    def id(self):
        """ID único del proveedor — no debe cambiar entre versiones."""
        return 'muestreo_espacial_puntos'

    def name(self):
        """Nombre visible en el panel de Processing."""
        return 'Muestreo Espacial de Puntos'

    def longName(self):
        """Nombre largo (descripción en el panel de Processing)."""
        return 'Muestreo Espacial de Puntos por Curva de Hilbert'

    def loadAlgorithms(self, *args, **kwargs):
        """Registra todos los algoritmos del proveedor."""
        self.addAlgorithm(MuestreoEspacialPuntos())
        # Aquí se pueden agregar algoritmos adicionales en versiones futuras.

    def icon(self):
        """Ícono del proveedor — QGIS usa el predeterminado si se retorna QIcon()."""
        return QIcon()
