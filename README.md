# Distortionizer

## Instructions on how to get it to build and run in VS

1. Download and install latest version of Qt (with openGL) 
2. Download and unpack latest version of GLEW
3. Open CMake and navigate to `distortionizer` directory
4. Add path for build directory like `/distortionizer/bin`
5. Click "Add Entry", add `CMAKE_PREFIX_PATH` of type PATH with the value containing the location of Qt OpenGL directory and Glew directory like `C:/Qt/Qt5.4.0-gl/5.4/msvc2013_opengl;C:/glew`
6. Click "Configure", and choose the right generator (probably Visual Studio 12 2013).
7. If there is an error, the "Generate" button will not be available. **Just having red entries doesn't mean it failed - you'd see an error message down below.** If it did fail, you'll need to fill in more entries, depending on what it failed to find.  For example:
  - For `Qt5Widgets_DIR` specify like (adjust for your path) `C:/Qt/Qt5.4.0-gl/5.4/msvc2013_opengl/lib/cmake/Qt5Widgets`
  - For `Qt5OpenGL_DIR` specify like (adjust for your path) `C:/Qt/Qt5.4.0-gl/5.4/msvc2013_opengl/lib/cmake/Qt5OpenGL`
  - For `GLEW_DLL` specify like (adjusting for your path) `C:/glew/bin/Release/Win32/glew32.dll`
8. Once you've gotten it to successfully configure, click "Generate" to make the VS solution in the build directory.
9. After you open solution, right click on `distortionizer-calibration` in VS in Solution Explorer widget and select "Set as Startup Project"
10. In menu header, click on Project -> Properties (or Alt+F7) and in "Configuration Properties" -> Debugging set "Environment" to something like `PATH=$(PATH);C:\Qt\Qt5.4.0-gl\5.4\msvc2013_opengl\bin` otherwise it will complain about missing `Qt5OpenGL.dll` file


##Instructions on how to use Distortionizer
Distortion estimation for HMD using K1 (quadratic) term
The program always runs on the last screen, full screen
Click with the mouse in the left eye to locate center
Keyboard controls:
  - r/R: Increase/Decrease distortion in R+G+B
  - g/G: Increase/Decrease distortion in G+B
  - b/B: Increase/Decrease distortion in B only
  - S/L: Save/Load state from JSON config file (HMD_Config.json by default)
  - Left,Right,Up,Down: Move the center of projection by one pixel
  - f/F:    Toggle fullscreen on/off
  - c/C:    Reset center of projection
  - v/V:    Reset distortion values to 0
  - ESC/Q: Quit the application


