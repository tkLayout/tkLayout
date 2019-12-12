This is to simulate resolution on pixel modules local coordinates X and Y with tkLayout.
2 modes are avaialble : either resolutionLocalX or Y are set as nominal values, either they are calculated from parameterized models.



* Nominal resolutionLocalX or Y :  
Please specify values for resolutionLocalX or Y in cfg file.
Example : nominalResolutionLocalX 0.010

Please note that if no value nor any parameter is specified in the cfg file, a nominal resolutionLocalX or Y is calculated using the following formulae :
nominalResolutionLocalX = RMS over sensors of (meanWidth / (numStripsAcross * sqrt(12)))
nominalResolutionLocalY = RMS over sensors of (length / (maxSegments * sqrt(12)))



* Parameterized resolutionLocalX or Y :  
Please specify the set of parameters for resolutionLocalX or Y in cfg file.
Please note that as soon as any parameter is specified for a given module type and local coordinate, then the parametric mode is chosen to calculate the corresponding resolution.

Example : @include-std CMS_Phase2/Pixel/Resolutions/100x100
Sets of parameters are available for the following module types : 25 um x 100 um, 50 um x 50 um, 50 um x 200 um, 100 um x 100 um.

Initial study was conducted in [1] and [2].
Model being now presently used is from Riccardo del Burgo, and mentioned below.

Resolution on local X axis:
B = Magnetic field intensity projected on local (XY) plane
x = fabs(-cotanAlpha - hallMobility*B);

Resolution on local Y axis:
x = fabs(-cotanBeta);

Model for resolution on local X or local Y coordinate:
resolutionLocalAxis = param0 + param1 * x 
                    + param2 * exp(-param9 * x) * cos(param3 * x + param4)
                    + param5 * exp(-0.5 * pow(((x - param6) / param7), 2))
                    + param8 * pow(x, 0.5);



NB on geometry :
In the plane of any module, local Y axis is always along the strips, and local X axis is the axis orthogonal to it.
As a result : 
Barrel modules : Local X axis is along RPhi, Local Y axis is along Z.
Endcap modules : Local X axis is along RPhi, Local Y axis is along R.

Specific to parametric case : alpha and beta angles are the reference notation angles for a track hitting a pixel module.
They are defined in Figure 7 from [3]. 



[1] Parametrization of the spatial resolution of the reconstructed hits in the Inner Pixel for HL-LHC studies, E. Migliore and M. Musich, 2015/09/18
[2] Summary of available parametrizations, E. Migliore, Torino, 2016/06/13
[3] Commissioning and Performance of the CMS Pixel Tracker with Cosmic Ray Muons, The CMS Collaboration, 2010/02/02
