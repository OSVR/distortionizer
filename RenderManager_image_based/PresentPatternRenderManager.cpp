/** @file
    @brief Uses RenderManager to display a calibration pattern on one or
           more screens (defined by the OSVR server configuration file).

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
#include <osvr/ClientKit/InterfaceStateC.h>
#include <osvr/Client/RenderManagerConfig.h>
#include "osvr/RenderKit/RenderManager.h"

// Needed for render buffer calls.  OSVR will have called glewInit() for us
// when we open the display.
#include <GL/glew.h>

// Library/third-party includes
#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif
#include <GL/GL.h>
#include <GL/GLu.h>

// Standard includes
#include <iostream>
#include <string>
#include <stdlib.h> // For exit()
#include <chrono>
#include <thread>

//This must come after we include <GL/GL.h> so its pointer types are defined.
#include "osvr/RenderKit/GraphicsLibraryOpenGL.h"

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

// Quadric and sphere to use for rendering.
static GLUquadric *sphere = nullptr;

// X,Y location
typedef struct {
  double x;
  double y;
} XY;

static std::vector<float> params;  //< Distortion parameters
static int activeParam = 0;  //< Which parameter are we adjusting?

static std::string osvrGetString(OSVR_ClientContext context, const std::string& path)
{
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

static float red_col[] = { 1.0, 0.0, 0.0 };
static float grn_col[] = { 0.0, 1.0, 0.0 };
static float blu_col[] = { 0.0, 0.0, 1.0 };
static float yel_col[] = { 1.0, 1.0, 0.0 };
static float cyn_col[] = { 0.0, 1.0, 1.0 };
static float mag_col[] = { 1.0, 0.0, 1.0 };
static float wht_col[] = { 1.0, 1.0, 1.0 };

bool SetupRendering(osvr::renderkit::GraphicsLibrary library)
{
  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (library.OpenGL == nullptr) {
    std::cerr << "SetupRendering: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return false;
  }

  osvr::renderkit::GraphicsLibraryOpenGL *glLibrary = library.OpenGL;

  GLfloat matspec[4] = { 0.0, 0.0, 0.0, 0.0 };
  glMaterialfv(GL_FRONT, GL_SPECULAR, matspec);
  glMaterialf(GL_FRONT, GL_SHININESS, 64.0);

  // Turn on depth testing, so we get correct ordering.
  glEnable(GL_DEPTH_TEST);

  // Set up lighting so the spheres will be brightest in the center
  // and fall off.
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
  GLfloat specularLight[] = { 0.0f, 0.0f, 0.0f, 1.0f };
  GLfloat position[] = { 0.0f, 0.0f, -2.0e5f, 1.0f };

  // Assign created components to GL_LIGHT0.
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, position);

  glEnable(GL_COLOR_MATERIAL);

  // Construct the quadric for the sphere primitive to draw to show spacing
  sphere = gluNewQuadric();

  return true;
}

// Render the world from the specified point of view.
void RenderView(
  size_t whichEye,  //< Which eye are we rendering?
  const OSVRDisplayConfiguration &displayConfiguration, //< Info we need about display
  const osvr::client::RenderManagerConfigPtr renderManagerConfig, //< Info we need about overfill
  const osvr::renderkit::RenderInfo &renderInfo,  //< Info needed to render
  GLuint frameBuffer, //< Frame buffer object to bind our buffers to
  GLuint colorBuffer, //< Color buffer to render into
  GLuint depthBuffer, //< Depth buffer to render into
  XY const &xSphere,  //< Where to draw the X-axis-marking sphere
  XY const &ySphere,  //< Where to draw the Y-axis-marking sphere
  std::vector<XY> const &spheres,  //< Spheres to draw
  float const *color, //< Color to draw the main spheres
  float radius  //< Radius of the spheres
  )
{
  // Make sure our pointers are filled in correctly.  The config file selects
  // the graphics library to use, and may not match our needs.
  if (renderInfo.library.OpenGL == nullptr) {
    std::cerr << "RenderView: No OpenGL GraphicsLibrary, this should not happen" << std::endl;
    return;
  }

  // Render to our framebuffer
  glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

  // Set color and depth buffers for the frame buffer
  glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
    colorBuffer, 0);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
    GL_RENDERBUFFER, depthBuffer);

  // Set the list of draw buffers.
  GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
  glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers

  // Always check that our framebuffer is ok
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    std::cerr << "RenderView: Incomplete Framebuffer" << std::endl;
    return;
  }

  // Set the viewport to cover our entire render texture.
  glViewport(0, 0,
    static_cast<GLsizei>(renderInfo.viewport.width),
    static_cast<GLsizei>(renderInfo.viewport.height));

  // Set the OpenGL projection matrix 
  GLdouble projection[16];
  osvr::renderkit::OSVR_Projection_to_OpenGL(projection, renderInfo.projection);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMultMatrixd(projection);

  /// Put the transform into the OpenGL ModelView matrix
  GLdouble modelView[16];
  osvr::renderkit::OSVR_PoseState_to_OpenGL(modelView, renderInfo.pose);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixd(modelView);

  // Clear the screen to black and clear depth
  glClearColor(0, 0, 0, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // =================================================================
  // This is where we draw our world and hands and any other objects.
  // We're in World Space.  To find out about where to render objects
  // in OSVR spaces (like left/right hand space) we need to query the
  // interface and handle the coordinate tranforms ourselves.

  // =================================================================
  // Now we want to draw things in screen space, so we construct new
  // projection and modelview matrices to do orthographic projection
  // into the overfill viewport.  We draw spheres at the specified locations
  // with respect to the original (non-overfill) viewport, such that
  // pixels are square and the entire width of the screen is the size
  // of the overfill zone.

  double overFill = renderManagerConfig->getRenderOverfillFactor();
  double width = displayConfiguration.getDisplayWidth();
  double height = displayConfiguration.getDisplayHeight();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-width / 2 * overFill, width / 2 * overFill,
           -height / 2 * overFill, height / 2 * overFill,
           100, -100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Draw a set of spheres at the specified locations in viewport space.
  // They are offset so that (0,0) is at the center of projection for
  // the eye.
  glColor3f(color[0], color[1], color[2]);
  double xCOP = displayConfiguration.getEyes()[whichEye].m_CenterProjX;
  double yCOP = displayConfiguration.getEyes()[whichEye].m_CenterProjY;
  double xOffset = xCOP - 0.5;
  double yOffset = yCOP - 0.5;
  for (size_t i = 0; i < spheres.size(); i++) {
    glPushMatrix();
    glTranslated(spheres[i].x + xOffset, spheres[i].y + yOffset, 0);
    gluSphere(sphere, radius, 40, 40);
    glPopMatrix();
  }

  // Draw the axis spheres at the specified locations in viewport space.
  // They are offset so that (0,0) is at the center of projection for
  // the eye.
  // Draw the x sphere in red and the y in green.
  glPushMatrix();
    glColor3f(red_col[0], red_col[1], red_col[2]);
    glTranslated(xSphere.x + xOffset, xSphere.y + yOffset, 0);
    gluSphere(sphere, radius/2, 40, 40);
  glPopMatrix();
  glPushMatrix();
    glColor3f(grn_col[0], grn_col[1], grn_col[2]);
    glTranslated(ySphere.x + xOffset, ySphere.y + yOffset, 0);
    gluSphere(sphere, radius/2, 40, 40);
  glPopMatrix();
}

void Usage(std::string name)
{
  std::cerr << "Usage: " << name
    << " [color (one of red, green, blue, white, cyan, magenta, yellow)]"
    << std::endl;
  exit(-1);
}

int main(int argc, char *argv[])
{
    // Parse the command line
    std::string colorName = "red";
    const float *color = red_col;
    int realParams = 0;
    for (int i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        Usage(argv[0]);
      }
      else switch (++realParams) {
      case 1:
        colorName = argv[i];
        break;
      default:
        Usage(argv[0]);
      }
    }
    if (realParams > 1) { Usage(argv[0]); }

    if (colorName == "red") {
      color = red_col;
    } else if (colorName == "green") {
      color = grn_col;
    } else if (colorName == "blue") {
      color = blu_col;
    } else if (colorName == "white") {
      color = wht_col;
    } else if (colorName == "cyan") {
      color = cyn_col;
    } else if (colorName == "magenta") {
      color = mag_col;
    } else if (colorName == "yellow") {
      color = yel_col;
    } else {
      std::cerr << "Unrecognized color: " << colorName << std::endl;
      Usage(argv[0]);
    }

    // Open RenderManager and set up the context for rendering to
    // an HMD.  Do this using the OSVR RenderManager interface,
    // which maps to the nVidia or other vendor direct mode
    // to reduce the latency.
    render = osvr::renderkit::createRenderManager(context.get(), "OpenGL");

    if ( (render == nullptr) ||
         (!render->doingOkay()) ) {
        std::cerr << "Could not create RenderManager" << std::endl;
        return 1;
    }

    // Wait until we get a connection to a display object, from which we will
    // read information that we need about display device resolutions and
    // distortion correction parameters.  Once we hear from the display
    // device, we presume that we will also be able to read our
    // RenderManager parameters.  Complain as we don't hear from the
    // display device.
    // @todo Verify that waiting for the display is sufficient to be
    // sure we'll get the RenderManager string.
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string displayInfo;
    do {
      context.update();
      displayInfo = osvrGetString(context.get(), "/display");
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed = end - start;
      if (elapsed.count() >= 1) {
        std::cerr << "Waiting to get Display from server..."
          << std::endl;
        start = end;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
    } while (displayInfo.size() == 0);

    // Get the display information from the context.
    OSVRDisplayConfiguration displayConfiguration(displayInfo);

    // Get the RenderManager pipeline configuration info
    osvr::client::RenderManagerConfigPtr renderManagerConfig = nullptr;
    try {
      renderManagerConfig = osvr::client::RenderManagerConfigFactory::createShared(context.get());
    }
    catch (std::exception& /*e*/) {
      std::cerr << "Could not parse /render_manager_parameters string from server." << std::endl;
      return 100;
    }
    if (renderManagerConfig == nullptr) {
      std::cerr << "Could not parse /render_manager_parameters string from server (NULL pipelineconfig)." << std::endl;
      return 101;
    }

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

    // Do a call to get the information we need to construct our
    // color and depth render-to-texture buffers.
    std::vector<osvr::renderkit::RenderInfo> renderInfo;
    context.update();
    renderInfo = render->GetRenderInfo();
    std::vector<osvr::renderkit::RenderBuffer> colorBuffers;
    std::vector<GLuint> depthBuffers; //< Depth/stencil buffers to render into

    // Use the width of the first eye to figure out the spacing for the
    // spheres based on the number requested across the viewport.
    // Both the radius and position are in pixels.
    const int numSpheresX = 20;
    int width = static_cast<int>(renderInfo[0].viewport.width);
    int height = static_cast<int>(renderInfo[0].viewport.width);
    const float sphereSpace = (width / 2.0f) / (numSpheresX - 0.5f);

    // Make sure we have enough spheres to cover the space in height
    // as well as width.
    int numSpheresY = static_cast<int>(0.5 + (numSpheresX * height) / width);
    std::vector<XY> spheres; //< Where to draw the spheres
    for (int x = -numSpheresX + 1; x < numSpheresX; x++) {
      for (int y = -numSpheresY + 1; y < numSpheresY; y++) {
        XY sphere;
        sphere.x = x * sphereSpace;
        sphere.y = y * sphereSpace;
        spheres.push_back(sphere);
      }
    }
    XY xSphere, ySphere;
    xSphere.x = sphereSpace / 2;
    xSphere.y = 0;
    ySphere.x = 0;
    ySphere.y = sphereSpace / 2;

    // Construct the buffers we're going to need for our render-to-texture
    // code.
    GLuint frameBuffer;               //< Groups a color buffer and a depth buffer
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    for (size_t i = 0; i < renderInfo.size(); i++) {

      // The color buffer for this eye.  We need to put this into
      // a generic structure for the Present function, but we only need
      // to fill in the OpenGL portion.
      //  Note that this must be used to generate a RenderBuffer, not just
      // a texture, if we want to be able to present it to be rendered
      // via Direct3D for DirectMode.  This is selected based on the
      // config file value, so we want to be sure to use the more general
      // case.
      //  Note that this texture format must be RGBA and unsigned byte,
      // so that we can present it to Direct3D for DirectMode
      GLuint colorBufferName = 0;
      glGenRenderbuffers(1, &colorBufferName);
      osvr::renderkit::RenderBuffer rb;
      rb.OpenGL = new osvr::renderkit::RenderBufferOpenGL;
      rb.OpenGL->colorBufferName = colorBufferName;
      colorBuffers.push_back(rb);

      // "Bind" the newly created texture : all future texture
      // functions will modify this texture glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, colorBufferName);

      // Determine the appropriate size for the frame buffer to be used for
      // this eye.
      int width = static_cast<int>(renderInfo[i].viewport.width);
      int height = static_cast<int>(renderInfo[i].viewport.height);

      // Give an empty image to OpenGL ( the last "0" means "empty" )
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
        width,
        height,
        0,
        GL_RGBA, GL_UNSIGNED_BYTE, 0);

      // Bilinear filtering
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

      // The depth buffer
      GLuint depthrenderbuffer;
      glGenRenderbuffers(1, &depthrenderbuffer);
      glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
      glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
        width,
        height);
      depthBuffers.push_back(depthrenderbuffer);
    }

    // Register our constructed buffers so that we can use them for
    // presentation.
    if (!render->RegisterRenderBuffers(colorBuffers)) {
      std::cerr << "RegisterRenderBuffers() returned false, cannot continue" << std::endl;
      quit = true;
    }

    // Continue rendering until it is time to quit.
    while (!quit) {
        // Update the context so we get our callbacks called and
        // update analog and button states.
        context.update();

        renderInfo = render->GetRenderInfo();

        // Render into each buffer using the specified information.
        // @todo Pass the color as a command-line argument
        for (size_t i = 0; i < renderInfo.size(); i++) {
          RenderView(i, displayConfiguration, renderManagerConfig,
            renderInfo[i], frameBuffer,
            colorBuffers[i].OpenGL->colorBufferName,
            depthBuffers[i],
            xSphere, ySphere,
            spheres, color, sphereSpace / 4);
        }

        // Send the rendered results to the screen
        if (!render->PresentRenderBuffers(colorBuffers, renderInfo)) {
          std::cerr << "PresentRenderBuffers() returned false, maybe because it was asked to quit" << std::endl;
          quit = true;
        }
    }

    // Clean up after ourselves.
    glDeleteFramebuffers(1, &frameBuffer);
    for (size_t i = 0; i < renderInfo.size(); i++) {
      glDeleteTextures(1, &colorBuffers[i].OpenGL->colorBufferName);
      delete colorBuffers[i].OpenGL;
      glDeleteRenderbuffers(1, &depthBuffers[i]);
    }

    // Close the Renderer interface cleanly.
    delete render;

    return 0;
}

