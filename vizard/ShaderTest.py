import viz
import vizconnect

viz.setMultiSample(4) 

import vizact
import vizshape
import vizfx.postprocess
import sys
import time
import threading
import json
from pprint import pprint

from vizhmd import HMD

class CustomHMD(HMD):
	def __init__(self,*args,**kw):
		HMD.__init__(self,*args,**kw)

json_config = open("HMD.json")
data = json.load(json_config)
                
"""
" Shader to be applied postprocess
"
"""
class DistortionEffect(vizfx.postprocess.BaseShaderEffect):

    def _getFragmentCode(self):
        frag_shader = ""
        with open ("ShaderTest.frag", "r") as frag:
           frag_shader=frag.read()
        return frag_shader

    def _getVertexCode(self):
        with open ("ShaderTest.vert", "r") as vert:
           vert_shader=vert.read()        
        return vert_shader

    def _createUniforms(self):
        
        self.uniforms.addFloat('k1_red', data["hmd"]["distortion"]["k1_red"])
        self.uniforms.addFloat('k1_green', data["hmd"]["distortion"]["k1_green"])
        self.uniforms.addFloat('k1_blue', data["hmd"]["distortion"]["k1_blue"])
        self.uniforms.addFloat('fullscr_center', [data["hmd"]["eyes"][0]["center_proj_x"], data["hmd"]["eyes"][0]["center_proj_y"]])
        self.uniforms.addFloat('left_center', [.25, .5])
        self.uniforms.addFloat('right_center', [.75, .5])
        json_config.close()
        
global effect
effect = DistortionEffect()
vizfx.postprocess.addEffect(effect)

viz.setMultiSample(4)

screenMode=viz.FULLSCREEN

if data["hmd"]["resolutions"][0]["display_mode"] != "full_screen":
    if data["hmd"]["resolutions"][0]["display_mode"]== "vert_side_by_side":
        screenMode=viz.FULLSCREEN | viz.STEREO_VERT
    else:
        screenMode=viz.FULLSCREEN | viz.STEREO_HORZ

hmd=CustomHMD(data["hmd"]["field_of_view"]["monocular_horizontal"],
    data["hmd"]["field_of_view"]["monocular_vertical"],
    overlap=data["hmd"]["field_of_view"]["overlap_percent"],
    leftRollShift=data["hmd"]["rendering"]["left_roll"],
    rightRollShift=data["hmd"]["rendering"]["right_roll"],
    verticalShift=data["hmd"]["field_of_view"]["pitch_tilt"],
    stereo=screenMode);
    



viz.window.setFullscreenMonitor(1)
viz.enable(viz.AUTO_COMPUTE)

viz.go(screenMode)

#Create skylight
sky_light = viz.addLight(euler=(0,0,0))
sky_light.position(0,0,0,0)
sky_light.color(viz.WHITE)
sky_light.ambient([5,5,5]) # necessary to light up images in Vizard 5

#Add the gallery mode
gallery = viz.addChild('gallery.osgb')
