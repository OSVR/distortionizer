# Vizard shaders

## Instructions on how to get it running

1. Download and install latest Python version
2. Download and install latest version of Vizard software (free version is ok)
3. You can use the distortionizer tool to generate HMD_Config.json file with 
distortion values and center of projection, and use the values from it to
copy over to HMD.json for vizard. Otherwise, you can use the provided 
sample HDM.json file to edit the values.
4. Depending on configuration settings in HMD.json, the vizard
will display in a fullscreen or split screen mode
For example: in HMD.json 
to set it to vertical split screen edit line as "display_mode": "vert_side_by_side"
to set it to full screen mode edit line as "display_mode": "full_screen"            
5. You can use the mouse to look around in the generated room.