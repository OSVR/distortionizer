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
#include <osvr/ClientKit/InterfaceStateC.h>
#include <osvr/Client/RenderManagerConfig.h>
#include "osvr/RenderKit/RenderManager.h"
#include "types.h"
#include "helper.h"

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

// Forward declarations of rendering functions defined below.
void draw_cube(double radius);

// Set this to true to cause the program to exit.
bool quit = false;

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
  GLuint frameBuffer,       //< Frame buffer object to bind our buffers to
  GLuint colorBuffer,       //< Color buffer to render into
  GLuint depthBuffer,       //< Depth buffer to render into
  XY const &leftForward,    //< Normalized screen coordinates of forwards in left eye
  XY const &rightForward    //< Normalized screen coordinates of forwards in right eye
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

  // Draw a cube with a 5-meter radius as the room we are floating in.
  draw_cube(5.0);

  // Draw another cube 1 meter along the -Z axis
  glTranslated(0, 0, -1);
  draw_cube(0.1);

  // =================================================================
  // Now we want to draw things in screen space, so we construct new
  // projection and modelview matrices to do orthographic projection
  // into the overfill viewport.  We draw spheres at the specified locations
  // with respect to the original (non-overfill) viewport, such that
  // pixels are square, the origin is at (0,0) and the width and height
  // match the overfilled width and height.  The origin of the screen
  // is at the middle.

  // @todo Handle overfill.
  double overFill = renderManagerConfig->getRenderOverfillFactor();
  double width = renderInfo.viewport.width;
  double height = renderInfo.viewport.height;
  double overWidth = width * overFill;
  double overHeight = height * overFill;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-overWidth / 2, overWidth / 2,
    -overHeight / 2, overHeight / 2,
    -1000, 1000);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Draw a sphere representing the "forwards" direction in whichever
  // eye this is.  Scale this normalized value to cover the actual
  // screen size.
  double xForward, yForward;
  if (whichEye == 0) {
    xForward = leftForward.x;
    yForward = leftForward.y;
  } else {
    xForward = rightForward.x;
    yForward = rightForward.y;
  }

  // Draw the forwards sphere at the specified locations in viewport space.
  // Change from normalized (0-1) coordinates into absolute viewport coordinates.
  // Draw into the original viewport space, not the oversized viewport.
  glColor3d(1, 1, 1);
  double xOffset = -0.5;
  double yOffset = -0.5;
  double xScale = width;
  double yScale = height;

  glPushMatrix();
    glTranslated((xForward + xOffset) * xScale, (yForward + yOffset) * yScale, 0);
    gluSphere(sphere, width*0.01, 30, 30);
  glPopMatrix();
  
  // Draw a set of horizontal and vertical lines in the right eye
  // Draw into the original viewport space, not the oversized viewport.
  glColor3d(0, 0, 0);
  glBegin(GL_LINES);
  if (whichEye == 1) for (double ofs = -width; ofs <= width; ofs += width/50) {
    glVertex2d(-width, ofs);
    glVertex2d(width, ofs);

    glVertex2d(ofs, -width);
    glVertex2d(ofs,  width);
  }
  glEnd();
}

void Usage(std::string name)
{
  std::cerr << "Usage: " << name
    << " [-eye right|left] (default is right)"
    << " [-depth_meters D] (default is 2.0)"
    << " [-field_angles] (use field angles, default is longitude and latitude)"
    << " [-mm] (screen distance units in the config file, default is meters)"
    << " [-screen screen_left_meters screen_bottom_meters screen_right_meters screen_top_meters]"
    << " (default auto-compute based on ranges seen)"
    << std::endl
    << "  This program reads from standard input a configuration that has a list of" << std::endl
    << "x,y screen coordinates in meters followed by long,lat angles in" << std::endl
    << "degrees where (0,0) is straight ahead from the eye, positive" << std::endl
    << "longitude is left and positive latitude is up." << std::endl
    << "  It displays information helpful to debugging the resulting distortion." << std::endl
    << std::endl;
  exit(1);
}

int main(int argc, char *argv[])
{
  // Parse the command line
  bool testingForwards = true;
  bool useRightEye = true;
  bool computeBounds = true;
  bool useFieldAngles = false;
  double left, right, bottom, top;
  double depth = 2.0;
  double toMeters = 1.0;
  int realParams = 0;
  for (int i = 1; i < argc; i++) {
    if (std::string("-mm") == argv[i]) {
      toMeters = 1e-3;  // Convert input in millimeters to meters
    } else if (std::string("-depth_meters") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      depth = atof(argv[i]);
    } else if (std::string("-field_angles") == argv[i]) {
      useFieldAngles = true;
    } else if (std::string("-screen") == argv[i]) {
      computeBounds = false;
      if (++i >= argc) { Usage(argv[0]); }
      left = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      bottom = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      right = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      top = atof(argv[i]);
    } else if (std::string("-eye") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      std::string eye = argv[i];
      if (eye == "left") {
        useRightEye = false;
      } else if (eye == "right") {
        useRightEye = true;
      } else {
        std::cerr << "Bad value for -eye: " << eye << ", expected left or right" << std::endl;
        Usage(argv[0]);
      }
    } else if ((argv[i][0] == '-') && (atof(argv[i]) == 0.0)) {
      Usage(argv[0]);
    }
    else switch (++realParams) {
    case 1:
    default:
      Usage(argv[0]);
    }
  }
  if (realParams != 0) { Usage(argv[0]); }

  //====================================================================
  // Parse the angle-configuration information from standard input.  Expect white-space
  // separation between numbers and also between entries (which may be on separate
  // lines).
  std::vector<Mapping> mapping;
  while (!std::cin.eof()) {
    // Read the mapping info from the input file.
    Mapping map;
    std::cin >> map.xyLatLong.longitude >> map.xyLatLong.latitude >> map.xyLatLong.x >> map.xyLatLong.y;
    mapping.push_back(map);
  }
  // There will have been one extra added, when running into EOF.
  mapping.pop_back();
  if (mapping.size() == 0) {
    std::cerr << "Error: No input points found" << std::endl;
    return 2;
  }

  //====================================================================
  // If we've been asked to auto-range the screen coordinates, compute
  // them here.
  if (computeBounds) {
    left = right = mapping[0].xyLatLong.x;
    bottom = top = mapping[0].xyLatLong.y;
    for (size_t i = 1; i < mapping.size(); i++) {
      double x = mapping[i].xyLatLong.x;
      double y = mapping[i].xyLatLong.y;
      if (x < left) { left = x; }
      if (x > right) { right = x; }
      if (y < bottom) { bottom = y; }
      if (y > top) { top = y; }
    }
    left *= toMeters;
    right *= toMeters;
    bottom *= toMeters;
    top *= toMeters;
  }

  // Open RenderManager and set up the context for rendering to
  // an HMD.  Do this using the OSVR RenderManager interface,
  // which maps to the nVidia or other vendor direct mode
  // to reduce the latency.
  render = osvr::renderkit::createRenderManager(context.get(), "OpenGL");

  if ((render == nullptr) ||
    (!render->doingOkay())) {
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

  // Read the analog trigger, which will let us increase
  // or decrease our D parameters for distortion correction.
  osvr::clientkit::Interface analogTrigger =
    context.getInterface("/controller/trigger");

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

  // Locate the normalized screen coordinates for the "forwards"
  // direction in the left and right screens; this is the direction
  // where the angles are both 0.
  // @todo if exact sample not found, interpolate between nearest three
  // non-collinear points.
  convert_to_normalized_and_meters(mapping, toMeters, depth,
    left, bottom, right, top, useFieldAngles);
  XY forward;
  forward.x = -1e5;    // Start out off-screen
  forward.y = -1e5;  // Start out off-screen
  for (size_t i = 0; i < mapping.size(); i++) {
    if ((mapping[i].xyLatLong.latitude == 0) &&
        (mapping[i].xyLatLong.longitude == 0)) {
      forward.x = mapping[i].xyLatLong.x;
      forward.y = mapping[i].xyLatLong.y;
    }
  }

  XY leftForward;
  XY rightForward;
  leftForward.y = rightForward.y = forward.y;
  if (useRightEye) {
    rightForward.x = forward.x;
    leftForward.x = 1 - forward.x;
  } else {
    leftForward.x = forward.x;
    rightForward.x = 1 - forward.x;
  }

    // If we're debugging the forwards direction, we need to turn off
    // distortion correction so that we draw the sphere directly into
    // screen coordinates.
    if (testingForwards) {
      osvr::renderkit::RenderManager::DistortionParameters distortionLeft;
      distortionLeft.m_desiredTriangles = 200 * 64;
      distortionLeft.m_type = osvr::renderkit::RenderManager::DistortionParameters
        ::rgb_symmetric_polynomials;
      std::vector<float> params;  //< Distortion parameters
      params.push_back(0);
      params.push_back(1);
      std::vector<float> Ds;
      Ds.push_back(1.0);
      Ds.push_back(1.0);
      distortionLeft.m_distortionD = Ds;
      distortionLeft.m_distortionPolynomialRed = params;
      distortionLeft.m_distortionPolynomialGreen = params;
      distortionLeft.m_distortionPolynomialBlue = params;
      distortionLeft.m_distortionCOP[0] =
        static_cast<float>(displayConfiguration.getEyes()[0].m_CenterProjX);
      distortionLeft.m_distortionCOP[1] =
        static_cast<float>(displayConfiguration.getEyes()[0].m_CenterProjY);

      osvr::renderkit::RenderManager::DistortionParameters distortionRight;
      distortionRight = distortionLeft;
      distortionRight.m_distortionCOP[0] =
        static_cast<float>(displayConfiguration.getEyes()[1].m_CenterProjX);
      distortionRight.m_distortionCOP[1] =
        static_cast<float>(displayConfiguration.getEyes()[1].m_CenterProjY);

      // Push the same distortion back for each eye.
      std::vector<osvr::renderkit::RenderManager::DistortionParameters> distortionParams;
      distortionParams.push_back(distortionLeft);
      distortionParams.push_back(distortionRight);

      render->UpdateDistortionMeshes(
        osvr::renderkit::RenderManager::DistortionMeshType::SQUARE,
        distortionParams);
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

        renderInfo = render->GetRenderInfo();

        // Render into each buffer using the specified information.
        for (size_t i = 0; i < renderInfo.size(); i++) {
          RenderView(i, displayConfiguration, renderManagerConfig,
            renderInfo[i], frameBuffer,
            colorBuffers[i].OpenGL->colorBufferName,
            depthBuffers[i],
            leftForward, rightForward);
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

