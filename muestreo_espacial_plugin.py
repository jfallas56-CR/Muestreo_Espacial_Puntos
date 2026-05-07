# -*- coding: utf-8 -*-
"""
Complemento QGIS: Muestreo_Espacial_Puntos
Clase principal del complemento — registra el proveedor de Processing.
"""

from qgis.core import QgsApplication
from .processing_provider import MuestreoEspacialProvider


class MuestreoEspacialPlugin:
    """Clase principal del complemento Muestreo Espacial de Puntos."""

    def __init__(self, iface):
        """Constructor.

        :param iface: Instancia de QgisInterface proporcionada por QGIS.
        """
        self.iface = iface
        self.provider = None

    def initGui(self):
        """Llamado cuando el complemento se activa en QGIS.
        Registra el proveedor de algoritmos de Processing.
        """
        self.provider = MuestreoEspacialProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        """Llamado cuando el complemento se desactiva.
        Elimina el proveedor del registro de Processing.
        """
        if self.provider:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
