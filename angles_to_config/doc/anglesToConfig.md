The source code for the **AnglesToConfig** program can be found in the [https://github.com/OSVR/distortionizer](https://github.com/OSVR/distortionizer) repository.  This program takes a list of mappings from view angle to physical display location and produces a fragment of an OSVR server configuration file that will un-distort the display, along with a description of the canonical screen to be rendered to.

The theory behind the operation of this program is described in the [distortion document](https://github.com/OSVR/OSVR-Docs/blob/master/Configuring/distortion.md).

## Step 1: Determine angle-to-screen mapping

Produce a whitespace-separated file with entries, four entries per line, with no header in the file.  The four columns represent: (1) the longitudinal angle in degrees from the forward gaze direction (positive angles towards the right), (2) the latitudinal angle in degrees from the forward gaze direction (positive angles up), (3) meters (or, with a command-line option, millimeters) from the forward-gaze point on the screen, positive towards the right, (4) meters (or, with a command-line option, millimeters) from the forward-gaze point on the screen, positive up.  The first few lines from an example file follow:

    -35	-45	-21.203		-30.274
    -30	-45	-17.754		-30.901
    -25	-45	-14.506		-31.304
    -20	-45	-11.451		-31.581

The _**forward-gaze direction**_ is the direction forward from the eye, as if looking at an object at infinity right in front of the eye.  The _**forward-gaze point**_ is the location view on the display that is seen at the center of the eye’s fovea when looking in the forward gaze direction.  Unless the corners of the screen are specified as a command-line argument, the program will use the range of X and Y values found within the file to implicitly determine the screen dimensions and the location of the forward-gaze point within the screen.

**Note:** The angles are field angles, which differ from latitude/longitude.  Each specifies the angle of a plane that tilts from the origin to pass through the specified point.  The longitude point is the YZ plane tilted towards the X axis (positive is towards +X) and the latitude point is the XZ plane tilted towards the Y axis (positive is towards +Y).  Earlier versions of the program uses lat/long spherical coordinates to determine the points, which produced extra distortions in the vertical direction.

**Note:** There can be a separate file for red, green and blue to enable chromatic distortion correction or there can be a single file that handles all colors for monochromatic distortion correction.

## Step 2: Compute the unstructured distortion mesh and canonical screen

**AnglesToConfig** reads in the table of unordered mappings from angles to locations on the physical display and produces a distortion map and canonical screen description for use in OSVR server configuration files.  The program is run as a filter, reading in the table on standard input and writing the resulting Json-formatted configuration file to standard output.  It can be run without command-line arguments, but there are optional arguments:
* **–screen screen_left_meters screen_bottom_meters screen_right_meters screen_top_meters** lets you specify the actual boundaries of the screen in case there are not points recorded at the actual boundaries.  This argument takes four parameters, each in meters: the left, bottom, top, and right location of the four corners of the screen with respect to the forward-gaze point.  The default is to determine these boundaries by finding the minimum and maximum x and y coordinate in the file.
* **-mm** specifies that the screen-position entries in the file are in millimeters.  The –screen parameters, if present, are still in meters.  The default is meters.
* **-eye right|left** specifies which eye, and has one parameter (left or right).  The default it right.
* **-depth_meters D** specifies the focal depth of the screen.  Larger distances reduce the impact of changes in IPD.  Note that the distortion is dependent on IPD.  The default is 2m.
* **-latlong** tells that the angles specified are in longitude and latitude, rather than field angles.  Field angles is the default.
* **-verify_angles xx xy yx yy max_degrees** tests each mesh point to ensure that the change between the vectors point to each of its neighbors in angle space, when transformed into screen space, does not differ by more than max_degrees.  The transformation is specified: The vector (xx, xy) points in screen space in the direction of +longitude (left).  The vector (yx, yy) points in screen space in the direction of +latitude (up).
* **-mono infile** takes the name of a file to read from rather than standard input, producing a monochromatic distortion function.
* **-rgb redfile greenfile bluefile** takes three file name arguments, one each for red, green, and blue.

**RGB Example:** One prototype display had a simulation run that produced three colored output files with angles from -40 degrees to 75 degrees in X for a right-eye display (-65 to 65 in Y).  The corresponding range in millimeters varied by color, but was near -27.3 to 73.3 in X and -33.4 to 33.4 in Y.  This meant that the size in X of the screen that was covered by simulation was less than the actual screen size.  The actual screen size was 120.96mm in X and 68.04mm in Y.  **Note:** In this case, the simulated region goes past the edge of the screen in the nasal direction.  This will cause a warning to be printed when the program is run and will produce out-of-bounds grid points that will go unused in the mesh.  The Y area was reported to be symmetric around the forward-gaze point, making the range of the screen -34.02mm to 34.02mm.  The X location of the forward-gaze point was reported to be 25.34mm from the left edge of the display.  This implies that the screen left edge is actually -25.34mm, and the right edge is 120.96mm from there, at 95.62mm.  The command line to support this file is (running with this file will produce a warning message):
* AnglesToConfig –mm –screen -0.02534 -0.03402 0.09562 0.03402 –rgb red_in.dat green_in.dat blue_in.dat > out.json

## Step 3: Constructing configuration files

**Main configuration file:** AnglesToConfig prints out a Json-format file that is a subset of the full configuration file that is required to send to an OSVR server program to support rendering to a display.

Several sections of this example file should be modified based on the information in the _out.json_ file produced in step 2.  The _field_of_view_ and _eyes_ sections should have their data fields changed to match (the new sections can be copy-pasted into a copy of this file, overwriting what is there).

**Distortion correction:** As of OSVR release 0.6, the size of configuration files is limited by a networking transfer unit size, so the distortion section of the file cannot be directly copied into the main configuration file.  Instead, it should remain inside the out.json file (which can be renamed or moved) and the main configuration file should be changed to point to it.  If the file is moved to C:/OSVR/Distortion_client.json, then the following section will describe how to load it on the client side:

    "distortion": {
        "type": "mono_point_samples",
     "mono_point_samples_external_file": "C:/OSVR/Distortion_client.json"
    },

For an RGB configuration file, use:

    "distortion": {
        "type": "rgb_point_samples",
     "rgb_point_samples_external_file": "C:/OSVR/Distortion_client.json"
    },

This can be copy-pasted over the distortion section in the main configuration file.

**Locating displays:** The width and height in the resolutions section of the display configuration should match the landscape size of the display, often 1920x1080.

The actual display orientation determines the value for the _rotation_ field within the _display_ section in the _renderManagerConfig_ part of the configuration file.  If the displays are in landscape mode, this will be 0.  If they are in portrait mode (which is faster to scan out for some displays), the value will be 90.

The actual location of the displays within the Microsoft Windows display layout should be described in the _xPosition_ and _yPosition_ fields within the window section of the _renderManagerConfig_ part of the configuration file.  The upper-left location for the left display should be described in this section.  The upper-left corner of the right display will be just to the right of the left one (either 1920 or 1080 pixels shifted, depending on the rotation).  The display locations within Windows should be configured to line up in this order.

## Step 4: Running programs using the configuration files

The configuration files to be used for distortion correction can be copied into the C:/OSVR directory on the computer to which the display is attached.  If this is done, and if the edited main configuration file is named Distortion_server.json, then the following command-line argument will run the server:

    osvr_server.exe C:/OSVR/Distortion_server.json

Any RenderManager-based OSVR client program can now be run, and it will get the configuration information from the server, reading the distortion-correction information from the _Distortion_client.json_ file locally.



