/** @file
    @brief Example program that uses the OSVR direct-to-display interface
           and OpenGL to render a scene with low latency.

    @date 2015

    @author
    Russ Taylor working through ReliaSolve.com for Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2015 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Internal Includes
#include <osvr/ClientKit/Context.h>
#include <osvr/ClientKit/Interface.h>
#include <osvr/ClientKit/InterfaceStateC.h>
#include <osvr/Client/RenderManagerConfig.h>
#include "osvr/RenderKit/RenderManager.h"

// Library/third-party includes
#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif
#include <GL/GL.h>

// Standard includes
#include <iostream>
#include <string>
#include <stdlib.h> // For exit()
#include <chrono>

//This must come after we include <GL/GL.h> so its pointer types are defined.
#include "osvr/RenderKit/GraphicsLibraryOpenGL.h"

// Forward declarations of rendering functions defined below.
void draw_cube(double radius);


// Set to true when it is time for the application to quit.
// Handlers below that set it to true when the user causes
// any of a variety of events so that we shut down the system
// cleanly.  This only works on Windows, but so does D3D...
static bool quit = false;

// Get an OSVR client context to use to access the devices
// that we need.
static osvr::clientkit::ClientContext context(
  "com.osvr.renderManager.distortionCorrectRenderManager");

static osvr::renderkit::RenderManager *render = nullptr;
static double scalePower = 0.0;


static std::string osvrRenderManagerGetString(OSVR_ClientContext context, const std::string& path) {
  size_t len;
  if (osvrClientGetStringParameterLength(context, path.c_str(), &len) == OSVR_RETURN_FAILURE) {
    std::string msg = std::string("Couldn't get osvr string length for path ") + path;
    std::cerr << msg << std::endl;
    throw std::runtime_error(msg);
  }

  //struct TempBuffer {
  //    char *buffer;
  //    TempBuffer(size_t len) { buffer = new char[len + 1]; }
  //    ~TempBuffer() { delete[] buffer; }
  //} tempBuffer(len);
  std::vector<char> tempBuffer(len + 1);

  if (osvrClientGetStringParameter(context, path.c_str(), tempBuffer.data(), len + 1) == OSVR_RETURN_FAILURE) {
    std::string msg = std::string("Couldn't get osvr string buffer for path ") + path;
    std::cerr << msg << std::endl;
    throw std::runtime_error(msg);
  }
  return std::string(tempBuffer.data(), len);
}

#ifdef _WIN32
// Note: On Windows, this runs in a different thread from
// the main application.
static BOOL CtrlHandler(DWORD fdwCtrlType)
{
    switch (fdwCtrlType)
    {
        // Handle the CTRL-C signal. 
    case CTRL_C_EVENT:
        // CTRL-CLOSE: confirm that the user wants to exit. 
    case CTRL_CLOSE_EVENT:
    case CTRL_BREAK_EVENT:
    case CTRL_LOGOFF_EVENT:
    case CTRL_SHUTDOWN_EVENT:
        quit = true;
        return TRUE;
    default:
        return FALSE;
    }
}
#endif

void myButtonCallback(void * /*userdata*/, const OSVR_TimeValue * /*timestamp*/,
    const OSVR_ButtonReport *report)
{
  if (report->state == 1) {

    // Get the original distortion correction
    std::string jsonString = osvrRenderManagerGetString(context.get(), "/display");
    OSVRDisplayConfiguration displayConfiguration(jsonString);

    // Create a new set of distortion parameters that is the original one
    // with D parameters scaled by the current factor.
    float scale = static_cast<float>(pow(2.0, scalePower));
    std::cout << "XXX New scale for Ds = " << scale << std::endl;

    osvr::renderkit::RenderManager::DistortionParameters distortion;
    distortion.m_desiredTriangles = 200 * 64;
    std::vector<float> Ds;
    Ds.push_back(
      displayConfiguration.getDistortionDistanceScaleX());
    Ds.push_back(
      displayConfiguration.getDistortionDistanceScaleY());
    distortion.m_distortionD = Ds;
    for (size_t i = 0; i < distortion.m_distortionD.size(); i++) {
      distortion.m_distortionD[i] *= scale;
    }
    // @todo Remove this when COP is always in range 0..1
    for (size_t i = 0; i < distortion.m_distortionCOP.size(); i++) {
      distortion.m_distortionCOP[i] *= scale;
    }
    distortion.m_distortionPolynomialRed =
      displayConfiguration.getDistortionPolynomalRed();
    distortion.m_distortionPolynomialGreen =
      displayConfiguration.getDistortionPolynomalGreen();
    distortion.m_distortionPolynomialBlue =
      displayConfiguration.getDistortionPolynomalBlue();

    // Push the same distortion back for each eye.
    std::vector<osvr::renderkit::RenderManager::DistortionParameters> distortionParams;
    for (size_t i = 0; i < displayConfiguration.getEyes().size(); i++) {
      distortionParams.push_back(distortion);
    }

    // Send a new set of parameters to construct a distortion mesh.
    render->UpdateDistortionMesh(osvr::renderkit::RenderManager::DistortionMeshType::SQUARE,
      distortionParams);
  }
}

bool SetupRendering(osvr::renderkit::GraphicsLibrary library)
{
  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (library.OpenGL == nullptr) {
    std::cerr << "SetupRendering: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return false;
  }

    osvr::renderkit::GraphicsLibraryOpenGL *glLibrary = library.OpenGL;

    // Turn on depth testing, so we get correct ordering.
    glEnable(GL_DEPTH_TEST);

    return true;
}

// Callback to set up a given display, which may have one or more eyes in it
void SetupDisplay(
    void *userData              //< Passed into SetDisplayCallback
    , osvr::renderkit::GraphicsLibrary  library //< Graphics library context to use
    , osvr::renderkit::RenderBuffer     buffers //< Buffers to use
    )
{
  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (library.OpenGL == nullptr) {
    std::cerr << "SetupDisplay: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return;
  }
  if (buffers.OpenGL == nullptr) {
    std::cerr << "SetupDisplay: No OpenGL RenderBuffer, this should not happen" << std::endl;
    return;
  }

  osvr::renderkit::GraphicsLibraryOpenGL *glLibrary = library.OpenGL;

  // Clear the screen to black and clear depth
  glClearColor(0, 0, 0, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

// Callback to set up for rendering into a given eye (viewpoint and projection).
void SetupEye(
    void *userData              //< Passed into SetViewProjectionCallback
    , osvr::renderkit::GraphicsLibrary  library //< Graphics library context to use
    , osvr::renderkit::RenderBuffer     buffers //< Buffers to use
    , osvr::renderkit::OSVR_ViewportDescription viewport  //< Viewport set by RenderManager
    , osvr::renderkit::OSVR_ProjectionMatrix  projection  //< Projection matrix set by RenderManager
    , size_t    whichEye        //< Which eye are we setting up for?
    )
{
  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (library.OpenGL == nullptr) {
    std::cerr << "SetupEye: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return;
  }
  if (buffers.OpenGL == nullptr) {
    std::cerr << "SetupEye: No OpenGL RenderBuffer, this should not happen" << std::endl;
    return;
  }

    // We don't do anything here -- everthing has been configured for us
    // in the RenderManager.
}

// Callbacks to draw things in world space, left-hand space, and right-hand
// space.
void DrawWorld(
    void *userData              //< Passed into AddRenderCallback
    , osvr::renderkit::GraphicsLibrary  library //< Graphics library context to use
    , osvr::renderkit::RenderBuffer     buffers //< Buffers to use
    , osvr::renderkit::OSVR_ViewportDescription viewport  //< Viewport we're rendering into
    , OSVR_PoseState   pose     //< OSVR ModelView matrix set by RenderManager
    , osvr::renderkit::OSVR_ProjectionMatrix  projection  //< Projection matrix set by RenderManager
    , OSVR_TimeValue deadline   //< When the frame should be sent to the screen
    )
{
  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (library.OpenGL == nullptr) {
    std::cerr << "DrawWorld: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return;
  }
  if (buffers.OpenGL == nullptr) {
    std::cerr << "DrawWorld: No OpenGL RenderBuffer, this should not happen" << std::endl;
    return;
  }

  osvr::renderkit::GraphicsLibraryOpenGL *glLibrary = library.OpenGL;

  /// Draw a cube with a 5-meter radius as the room we are floating in.
  draw_cube(5.0);

  glTranslated(2, 0, 0);
  draw_cube(0.3);
}

// This can be used to draw a heads-up display.  Unlike in a non-VR game,
// this can't be drawn in screen space because it has to be at a consistent
// location for stereo viewing through potentially-distorted and offset lenses
// from the HMD.  This example uses a small cube drawn in front of us.
// NOTE: For a fixed-display set-up, you do want to draw in screen space.
void DrawHead(
  void *userData              //< Passed into AddRenderCallback
  , osvr::renderkit::GraphicsLibrary  library //< Graphics library context to use
  , osvr::renderkit::RenderBuffer     buffers //< Buffers to use
  , osvr::renderkit::OSVR_ViewportDescription viewport  //< Viewport we're rendering into
  , OSVR_PoseState   pose     //< OSVR ModelView matrix set by RenderManager
  , osvr::renderkit::OSVR_ProjectionMatrix  projection  //< Projection matrix set by RenderManager
  , OSVR_TimeValue deadline   //< When the frame should be sent to the screen
  )
{
  std::string *stringToPrint = static_cast<std::string *>(userData);

  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (library.OpenGL == nullptr) {
    std::cerr << "DrawHead: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return;
  }
  if (buffers.OpenGL == nullptr) {
    std::cerr << "DrawHead: No OpenGL RenderBuffer, this should not happen" << std::endl;
    return;
  }

  osvr::renderkit::GraphicsLibraryOpenGL *glLibrary = library.OpenGL;

  /// Draw a small cube in front of us.
  glTranslated(0, 0, -0.25);
  draw_cube(0.005);
}

int main(int argc, char *argv[])
{
    // Construct button devices and connect them to a callback
    // that will send new distortion parameters when
    // button "8" on the controller is pressed.
    osvr::clientkit::Interface button8 =
        context.getInterface("/controller/8");
    button8.registerCallback(&myButtonCallback, nullptr);

    // Read the analog trigger, which will let us increase
    // or decrease our D parameters for distortion correction.
    osvr::clientkit::Interface analogTrigger =
      context.getInterface("/controller/trigger");

    // Open Direct3D and set up the context for rendering to
    // an HMD.  Do this using the OSVR RenderManager interface,
    // which maps to the nVidia or other vendor direct mode
    // to reduce the latency.
    render = osvr::renderkit::createRenderManager(context.get(), "OpenGL");

    if ( (render == nullptr) ||
         (!render->doingOkay()) ) {
        std::cerr << "Could not create RenderManager" << std::endl;
        return 1;
    }

    // Set callback to handle setting up rendering in a display
    render->SetDisplayCallback(SetupDisplay);

    // Set callback to handle setting up rendering in an eye
    render->SetViewProjectionCallback(SetupEye);

    // Register callbacks to render things in head, left hand, right
    // hand, and world space.
    render->AddRenderCallback("/", DrawWorld);
    render->AddRenderCallback("/me/head", DrawHead);

    // Set up a handler to cause us to exit cleanly.
#ifdef _WIN32
    SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE);
#endif

    // Open the display and make sure this worked.
    osvr::renderkit::RenderManager::OpenResults ret = render->OpenDisplay();
    if (ret.status == osvr::renderkit::RenderManager::OpenStatus::FAILURE) {
        std::cerr << "Could not open display" << std::endl;
        delete render;
        return 2;
    }

    // Set up the rendering state we need.
    if (!SetupRendering(ret.library)) {
        return 3;
    }

    // Keep track of time so we can scale UI analog adjustments
    // properly.
    std::chrono::time_point<std::chrono::system_clock> lastTime;
    lastTime = std::chrono::system_clock::now();

    // Continue rendering until it is time to quit.
    while (!quit) {
        // Update the context so we get our callbacks called and
        // update analog and button states.
        context.update();

        // Read the current value of the analogs we want
        OSVR_TimeValue  ignore;
        OSVR_AnalogState triggerValue = 0;
        osvrGetAnalogState(analogTrigger.get(), &ignore, &triggerValue);

        // Adjust the distortion D parameters by scaling them by a changing scale
        // that can be controlled by the trigger.  The analog value will be from
        // -1 to 1, which we integrate after scaling, adjusting an initial factor
        // of 0 (2^0 = 1) up and down.  We want the right controller to make
        // the scale factor larger, so we negate it during integration.
        std::chrono::time_point<std::chrono::system_clock> now =
          std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_sec = now - lastTime;
        lastTime = now;
        if (triggerValue != 0) {
          scalePower -= elapsed_sec.count() * triggerValue / 10;
          std::cout << "New D scale: " << pow(2.0, scalePower) << std::endl;
        }

        if (!render->Render()) {
            std::cerr << "Render() returned false, maybe because it was asked to quit" << std::endl;
            quit = true;
        }
    }

    // Close the Renderer interface cleanly.
    delete render;

    return 0;
}

static GLfloat matspec[4] = { 0.5, 0.5, 0.5, 0.0 };
static float red_col[] = { 1.0, 0.0, 0.0 };
static float grn_col[] = { 0.0, 1.0, 0.0 };
static float blu_col[] = { 0.0, 0.0, 1.0 };
static float yel_col[] = { 1.0, 1.0, 0.0 };
static float lightblu_col[] = { 0.0, 1.0, 1.0 };
static float pur_col[] = { 1.0, 0.0, 1.0 };

void draw_cube(double radius)
{
    GLfloat matspec[4] = { 0.5, 0.5, 0.5, 0.0 };
    glPushMatrix();
    glScaled(radius, radius, radius);
    glMaterialfv(GL_FRONT, GL_SPECULAR, matspec);
    glMaterialf(GL_FRONT, GL_SHININESS, 64.0);
    glBegin(GL_POLYGON);
    glColor3fv(lightblu_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, lightblu_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, lightblu_col);
    glNormal3f(0.0, 0.0, -1.0);
    glVertex3f(1.0, 1.0, -1.0);
    glVertex3f(1.0, -1.0, -1.0);
    glVertex3f(-1.0, -1.0, -1.0);
    glVertex3f(-1.0, 1.0, -1.0);
    glEnd();
    glBegin(GL_POLYGON);
    glColor3fv(blu_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, blu_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, blu_col);
    glNormal3f(0.0, 0.0, 1.0);
    glVertex3f(-1.0, 1.0, 1.0);
    glVertex3f(-1.0, -1.0, 1.0);
    glVertex3f(1.0, -1.0, 1.0);
    glVertex3f(1.0, 1.0, 1.0);
    glEnd();
    glBegin(GL_POLYGON);
    glColor3fv(yel_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, yel_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, yel_col);
    glNormal3f(0.0, -1.0, 0.0);
    glVertex3f(1.0, -1.0, 1.0);
    glVertex3f(-1.0, -1.0, 1.0);
    glVertex3f(-1.0, -1.0, -1.0);
    glVertex3f(1.0, -1.0, -1.0);
    glEnd();
    glBegin(GL_POLYGON);
    glColor3fv(grn_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, grn_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, grn_col);
    glNormal3f(0.0, 1.0, 0.0);
    glVertex3f(1.0, 1.0, 1.0);
    glVertex3f(1.0, 1.0, -1.0);
    glVertex3f(-1.0, 1.0, -1.0);
    glVertex3f(-1.0, 1.0, 1.0);
    glEnd();
    glBegin(GL_POLYGON);
    glColor3fv(pur_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, pur_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, pur_col);
    glNormal3f(-1.0, 0.0, 0.0);
    glVertex3f(-1.0, 1.0, 1.0);
    glVertex3f(-1.0, 1.0, -1.0);
    glVertex3f(-1.0, -1.0, -1.0);
    glVertex3f(-1.0, -1.0, 1.0);
    glEnd();
    glBegin(GL_POLYGON);
    glColor3fv(red_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, red_col);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, red_col);
    glNormal3f(1.0, 0.0, 0.0);
    glVertex3f(1.0, -1.0, 1.0);
    glVertex3f(1.0, -1.0, -1.0);
    glVertex3f(1.0, 1.0, -1.0);
    glVertex3f(1.0, 1.0, 1.0);
    glEnd();
    glPopMatrix();
}

