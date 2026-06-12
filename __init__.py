# -*- coding: utf-8 -*-

def classFactory(iface):
    """Carga la clase principal del complemento Muestreo_Espacial_Puntos.

    :param iface: Instancia de QgisInterface proporcionada por QGIS.
    """
    from .muestreo_espacial_plugin import MuestreoEspacialPlugin
    return MuestreoEspacialPlugin(iface)
