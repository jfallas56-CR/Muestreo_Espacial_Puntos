# -*- coding: utf-8 -*-
"""
Clase principal del complemento QGIS: Muestreo Espacial de Puntos.

Registra el algoritmo en Processing y agrega un botón en la barra de
herramientas que abre el diálogo del algoritmo directamente.

Arquitectura:
    __init__.py
        └── muestreo_espacial_plugin.py   ← este archivo
                └── muestreo_espacial_provider.py
                        └── Muestreo_Espacial_Puntos.py

Compatibilidad: QGIS 3.28 LTR / 3.44 LTR / 4.0  (Qt5 y Qt6)
Autor:          Jorge Fallas · jfallas56@gmail.com
Versión:        1.0.0
"""

import os

from qgis.PyQt.QtWidgets import QAction
from qgis.PyQt.QtGui import QIcon
from qgis.core import QgsApplication


class MuestreoEspacialPlugin:
    """Complemento QGIS que expone MuestreoEspacialPuntos en Processing."""

    TOOLBAR_NAME = 'Muestreo Espacial de Puntos'

    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self._action = None
        self._toolbar = None
        self._plugin_dir = os.path.dirname(__file__)

    # ------------------------------------------------------------------
    # Ciclo de vida
    # ------------------------------------------------------------------

    def initGui(self):
        """Registra el proveedor y crea el botón en la barra de herramientas."""

        # 1. Registrar proveedor de Processing
        from .muestreo_espacial_provider import MuestreoEspacialProvider
        self.provider = MuestreoEspacialProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

        # 2. Ícono
        icon_path = os.path.join(self._plugin_dir, 'icon.png')
        icon = QIcon(icon_path) if os.path.isfile(icon_path) else QIcon()

        # 3. Acción — abre el diálogo del algoritmo
        self._action = QAction(
            icon,
            'Muestreo Espacial de Puntos',
            self.iface.mainWindow()
        )
        self._action.setToolTip('Muestreo Espacial de Puntos')
        self._action.triggered.connect(self._open_dialog)

        # 4. Barra de herramientas dedicada
        self._toolbar = self.iface.addToolBar(self.TOOLBAR_NAME)
        self._toolbar.setObjectName(self.TOOLBAR_NAME)
        self._toolbar.addAction(self._action)

        # 5. Menú Complementos
        self.iface.addPluginToMenu(self.TOOLBAR_NAME, self._action)

    def unload(self):
        """Elimina proveedor, barra y menú al desactivar el complemento."""
        if self._action:
            self.iface.removePluginMenu(self.TOOLBAR_NAME, self._action)
            self._action = None

        if self._toolbar:
            self._toolbar.clear()
            self.iface.mainWindow().removeToolBar(self._toolbar)
            self._toolbar = None

        if self.provider:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None

    # ------------------------------------------------------------------
    # Slot
    # ------------------------------------------------------------------

    def _open_dialog(self):
        """Abre el diálogo Processing del algoritmo.

        El ID se construye dinámicamente desde la clase del algoritmo para
        evitar que quede desincronizado cuando cambia VERSION.
        Format: '{provider_id}:{algorithm_name}'
        """
        try:
            import processing
            from .Muestreo_Espacial_Puntos import MuestreoEspacialPuntos
            alg = MuestreoEspacialPuntos()
            alg_id = f"{self.provider.id()}:{alg.name()}"
            processing.execAlgorithmDialog(alg_id)
        except Exception as e:
            self.iface.messageBar().pushWarning(
                'Muestreo Espacial',
                f'No se pudo abrir el diálogo: {e}'
            )
