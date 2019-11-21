# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 22:06:00 2018

@author: kkc29
"""

import numpy as np
import sys

class FourierShellCorrelation:
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        
        visFr.AddMenuItem('Kenny>FSC', 'Calculate FSC', self.on_calculate_frc_from_locs,
                          helpText='')

    def on_calculate_frc_from_locs(self, event=None):
        from fsc.plugins.recipes.fsc import CalculateFRCFromLocs
         
        recipe = self.pipeline.recipe
        frc = CalculateFRCFromLocs(recipe,
                           inputName=self.pipeline.selectedDataSourceKey)
        
        if frc.configure_traits(kind='modal'):
            recipe.trait_set(execute_on_invalidation=False)
            recipe.add_module(frc)
        
            recipe.execute()
            
            self.visFr.CreateFoldPanel() #TODO: can we capture this some other way?
            
            print(recipe.namespace[frc.output_frc_dict])
            
            recipe.trait_set(execute_on_invalidation=True)
            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.fsc = FourierShellCorrelation(visFr)
    print("plugging FourierShellCorrelation as visFr.fsc")