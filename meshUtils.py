"""
 this class generates and manipulates archives of meshes
"""

import os
import glob


class meshArxv():
    
    def __init__(self, fName=None):
        self.fName = fName
        self.arxv = None
        
    def setArxv(self, arxv):
        self.arxv = arxv