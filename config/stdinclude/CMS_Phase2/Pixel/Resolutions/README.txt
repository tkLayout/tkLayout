// This is to simulate resolution on pixel modules local coordinates X and Y.
// resolutionLocalX or Y are either nominal values, either values calculated from parameterized models.




// Nominal resolutionLocalX or Y : Please specify values for resolutionLocalX or Y in this file.
// Example : nominalResolutionLocalX 0.010
// Please note that if no value nor any parameter is specified in the cfg file, a nominal resolutionLocalX or Y is calculated using the following formulae :
// nominalResolutionLocalX = RMS over sensors of (meanWidth / (numStripsAcross * sqrt(12)))
// nominalResolutionLocalY = RMS over sensors of (length / (maxSegments * sqrt(12)))




// Parameterized resolutionLocalX or Y : Please specify the set of parameters for resolutionLocalX or Y in this file. Models used are from [1] and mentioned below.
// Example : @includestd CMS_Phase2/Pixel/Resolutions/Barrel_100x100
// [1] Parametrization of the spatial resolution of the reconstructed hits in the Inner Pixel for HL-LHC studies, E. Migliore and M. Musich, 2015/09/18
// Sets of parameters are available for the following module types :
// SR : Small Rectangular (25 um x 100 um)
// SS : Small Square      (50 um x 50 um)
// LR : Large Rectangular (50 um x 200 um)
// LS : Large Square      (100 um x 100 um)
//
//
//
// Barrel PXB :
// If parameters to calculate resolutionLocalXBarrel or resolutionLocalYBarrel are specified, the following models are used :
// Model for resolution on local X coordinate (barrel modules):
// resolutionLocalXBarrel = resolutionLocalXBarrelParam0 + resolutionLocalXBarrelParam1 * cos(alpha) + resolutionLocalXBarrelParam2 * cos(alpha)^2
// Model for resolution on local Y coordinate (barrel modules):
// resolutionLocalYBarrel = resolutionLocalYBarrelParam0 + resolutionLocalYBarrelParam1 * exp(-resolutionLocalYBarrelParam2 * abs(cos(beta))) * sin(resolutionLocalYBarrelParam3 * abs(cos(beta)) + resolutionLocalYBarrelParam4)
// With alpha and beta the reference notation angles for a track hitting a pixel module.
//
//
//
// Endcap PXE :
// If parameters to calculate resolutionLocalXEndcap or resolutionLocalYEndcap are specified, the following models are used :
// Model for resolution on local X coordinate (endcap modules):
// resolutionLocalXEndcap = resolutionLocalXBarrelParam0 + resolutionLocalXBarrelParam1 * cos(alpha)
// Model for resolution on local Y coordinate (endcap modules):
// resolutionLocalYEndcap = resolutionLocalYBarrelParam0 + resolutionLocalYBarrelParam1 * abs(cos(beta))
// With alpha and beta the reference notation angles for a track hitting a pixel module.
