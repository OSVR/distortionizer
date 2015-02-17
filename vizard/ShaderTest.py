import viz
import vizconnect

viz.setMultiSample(4) 
#vizconnect.go('vizconnect_zSight.py') 
#vizconnect.go('vizconnect_config.py') 

import vizact
import vizshape
import sensics
import vizfx.postprocess
import socket
import SocketServer
import sys
import time
import threading
import json
from pprint import pprint


post_gaze = []
pre_gaze = []
recording = False
json_config = open("HMD_Config.json")
data = json.load(json_config)
"""
" Function to update textures based
" off block and trial. 
"""

                
"""
" BubbleZoom effect to be applied postprocess
"
"""
class BubbleEffect(vizfx.postprocess.BaseShaderEffect):

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
effect = BubbleEffect()
vizfx.postprocess.addEffect(effect)


viz.setMultiSample(4)
viz.fov(60)

#viz.window.setFullscreenMonitor(2)
viz.window.setSize(1280, 1024)
viz.enable(viz.AUTO_COMPUTE)
#viz.window.setPosition(1920,0) # Fullscreen wouldn't work without manually setting the position. Could be TeamViewer


if data["hmd"]["fullscreen"] == 1:
    viz.go()
else:
    viz.go(viz.STEREO_HORZ)
viz.MainView.move([0,0,0])


#Create skylight
sky_light = viz.addLight(euler=(0,0,0))
sky_light.position(0,0,0,0)
sky_light.color(viz.WHITE)
sky_light.ambient([5,5,5]) # necessary to light up images in Vizard 5

#Add the gallery mode
gallery = viz.addChild('gallery.osgb')

