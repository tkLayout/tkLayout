This is to simulate resolution on pixel modules local coordinates X and Y.
resolutionLocalX or Y are either nominal values, either values calculated from parameterized models.




* Nominal resolutionLocalX or Y : Please specify values for resolutionLocalX or Y in cfg file.
Example : nominalResolutionLocalX 0.010

Please note that if no value nor any parameter is specified in the cfg file, a nominal resolutionLocalX or Y is calculated using the following formulae :
nominalResolutionLocalX = RMS over sensors of (meanWidth / (numStripsAcross * sqrt(12)))
nominalResolutionLocalY = RMS over sensors of (length / (maxSegments * sqrt(12)))






* Parameterized resolutionLocalX or Y : Please specify the set of parameters for resolutionLocalX or Y in cfg file.
Example : @includestd CMS_Phase2/Pixel/Resolutions/Barrel_100x100

Sets of parameters are available for the following module types : 25 um x 100 um, 50 um x 50 um, 50 um x 200 um, 100 um x 100 um.

Initial study was conducted in [1].
Models being now presently used are from [2] and mentioned below.

- Barrel PXB :
If parameters to calculate resolutionLocalXBarrel or resolutionLocalYBarrel are specified, the following models are used :
Model for resolution on local X coordinate (barrel modules):
resolutionLocalXBarrel = resolutionLocalXBarrelParam0 + resolutionLocalXBarrelParam1 * cotg(alpha) + resolutionLocalXBarrelParam2 * cotg(alpha)^2
Model for resolution on local Y coordinate (barrel modules):
resolutionLocalYBarrel = resolutionLocalYBarrelParam0 + resolutionLocalYBarrelParam1 * exp(-resolutionLocalYBarrelParam2 * abs(cotg(beta))) * sin(resolutionLocalYBarrelParam3 * abs(cotg(beta)) + resolutionLocalYBarrelParam4)
With alpha and beta the reference notation angles for a track hitting a pixel module.


- Endcap PXE :
If parameters to calculate resolutionLocalXEndcap or resolutionLocalYEndcap are specified, the following models are used :
Model for resolution on local X coordinate (endcap modules):
resolutionLocalXEndcap = resolutionLocalXEndcapParam0 + resolutionLocalXEndcapParam1 * exp(-pow(cotg(alpha), 2.) / resolutionLocalXEndcapParam3) * cos(resolutionLocalXEndcapParam2 * cotg(alpha))
Model for resolution on local Y coordinate (endcap modules):
resolutionLocalYEndcap = resolutionLocalYBarrelParam0 + resolutionLocalYBarrelParam1 * abs(cotg(beta))
With alpha and beta the reference notation angles for a track hitting a pixel module.


[1] Parametrization of the spatial resolution of the reconstructed hits in the Inner Pixel for HL-LHC studies, E. Migliore and M. Musich, 2015/09/18
[2] Summary of available parametrizations, E. Migliore, Torino, 2016/06/13








NB on geometry :
In the plane of any module, local Y axis is always along the strips, and local X axis is the axis orthogonal to it.
As a result : 
Barrel modules : Local X axis is in RPhi, Local Y axis is along Z.
Endcap modules : Local X axis is in RPhi, Local Y axis is along R.

Alpha and beta angles are defined in Figure 7 from [3] Commissioning and Performance of the CMS Pixel Tracker with Cosmic Ray Muons, The CMS Collaboration, 2010/02/02.

