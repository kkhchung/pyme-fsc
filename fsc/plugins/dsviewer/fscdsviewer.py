# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 22:06:00 2018

@author: kkc29
"""

import numpy as np
import sys
import fsc

class FourierShellCorrelation:
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.image = dsviewer.image
        
        dsviewer.AddMenuItem('Kenny>FRC', 'Calculate FRC', self.on_calculate_frc_from_images,
                          helpText='')

    def on_calculate_frc_from_images(self, event=None):
        from fsc.plugins.recipes.fsc import CalculateFRCFromImages
        from PYME.recipes.base import ModuleCollection
         
        recipe = ModuleCollection()
        frc = CalculateFRCFromImages(recipe)
        
        if frc.configure_traits(kind='modal'):
            recipe.add_module(frc)
        
            recipe.execute()            
            
            print(recipe.namespace[frc.output_frc_dict])

def Plug(dsviewer):
#    """Plugs this module into the gui"""
#    dsviewer.fsc = FourierShellCorrelation(dsviewer)
#    print("plugging FourierShellCorrelation as visFr.fsc")
    pass