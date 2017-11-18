# ITK-BSpline-Registration
ITK BSpline registration processing time comparison between V3 and V4

I compared the processing time of ITK BSpline registration using V3 and V4.
Cardiac MRI images are used (20 cropped images of a slice). 
The results are shown below.

|Version| Processing Time [sec]  | CPU Utilization [%] |
|---| ---------------------- | ------------------- |
|V3| 2.4                    | 20                  |
|V4| 7.2                    | 100                 |
    
My PC environment:<br />
Windows 10<br />
Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz<br />
16.0 GB<br />
64-bit Operating System<br />

ITK version:<br />
4.12.2

VTK version:<br />
7.1.0

Installation:
1. Launch CMake
2. Select Configure
3. Set ITK_DIR
4. Set VIS_MODE (w/ or w/o visualization)
5. Select Generate
6. Select Open Project.
7. Build in Release.
8. Launch command prompt.
9. Go to the project directory (RegV3 or RegV4).
10. Run executables. <br />
Release\RegistrationV3.exe <br />
Release\RegistrationV4.exe <br />


