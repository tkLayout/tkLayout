Pixel models

Baseline_tilted_200_Pixel_1_1_1.cfg  Square pixels everywhere (small L1-2 and R1-2)
Baseline_tilted_200_Pixel_1_1_2.cfg  Like 1.1.1, but without DC/DC converters
Baseline_tilted_200_Pixel_1_2_1.cfg  Like 1.1.1, but Rectangular pixels in barrel
Baseline_tilted_200_Pixel_1_3_1.cfg  Like 1.1.1, but Rectangular everywhere
Baseline_tilted_200_Pixel_1_4_1.cfg  Like 1.3.1 But a bit shorter barrel
Baseline_tilted_200_Pixel_1_5_1.cfg  Like 1.3.1 But a lot shorter barrel
Baseline_tilted_200_Pixel_1_5_2.cfg  Like 1.5.1, but without DC/DC converters

Baseline_tilted_200_Pixel_3_5_1.cfg  2-rings step Like, TEDD and FPIX_2 disks moved to allow insertion -- small pixels everywhere

Baseline_tilted_200_Pixel_4_0_0.cfg  Like 3_5_1, but with 22x16.4 mm^2 chips (good for tilt?)
                                     BPIX=25x100 FPIX1=50x50 FPIX2=25x100

Baseline_tilted_200_Pixel_4_0_1.cfg  Like 4_0_0, but BPIX=50x50
                                     i.e. BPIX=50x50 FPIX1=50x50 FPIX2=25x100

Baseline_tilted_200_Pixel_4_0_2.cfg  Like 4_0_0, but FPIX1=25x100
                                     i.e. BPIX=25x100 FPIX1=25x100 FPIX2=25x100

Baseline_tilted_200_Pixel_4_0_3.cfg  Like 4_0_2, but
                                     FPIX2=50x50

Baseline_tilted_200_Pixel_4_0_2_1.cfg Like 4_0_2, but with 7 FPIX1 disks and 4 FPIX2 disks

OT_Tilted_360_200_Pixel_4021.cfg     First OT-numbered version: 3.6.0
                                     Pixel version 4.0.2.1
                                     Contains the new 2S barrel, with 3x multiplicity
                                     Contains new tilted barrel version 2016-07-15 (by kamil Cichy)

OT_Tilted_361_200_Pixel_4021.cfg     OT Version 3.6.1
                                     Pixel version 4.0.2.1
                                     Small adjustment in tilted barrel positions (version 2016-07-16 by Kamil Cichy)

OT_Tilted_362_200_Pixel_4021.cfg     OT Version 3.6.2
                                     Pixel version 4.0.2.1
                                     VERY small adjustment in positioning of TBPS layer 1 ring 8 and 10 (affexts rings 8+)

OT_Tilted_362_200_Pixel_4022.cfg     OT Version 3.6.2
                                     Pixel version 4.0.2.2 <- like 4.0.2.1 but with lower radii for BPIX L3 and L4 (less modules)

OT_Tilted_363_200_Pixel_4021.cfg     OT Version 3.6.3 !! Work in progress !!
                                     This is just an intermediary version towards an optimized TBPS. Main optimization point
                                     should be the reduction of TBPS planks δ and Δ, by 5 mm and 1 mm respectively
                                     Position of first ring (by changing zOverlap) was moved in such a way to keep the distance between the centre of the last flat module and the centre of the first
                                     tilted module (in z)
                          (3.6.2)    z_last_flat z_inner_first_tilted zOverlap zError_outer
                                     121.212     170.363              1.0      3.572
                                     213.069     265.499              1.0      6.262
                                     305.080     355.757              1.5      13.725
                          (.....)    z_last_flat z_inner_first_tilted zOverlap zError_outer
                                     124.029     166.827              18       57.977
                                     216.018     262.683              16       89.562
                                     307.947     352.500              17       138.486

                                     Radii and zOverlap of first rings were then changed so as to obtain zError_outer ~50mm, ~70mm, ~70mm for layers 1, 2, 3 at the transition
                                     Layer  Radius [mm]        zOverlap [mm]  zError_outer
                                     1      227.5   -> 228.0   18 -> 16       57.977   -> 52.231
                                     2      355.175 -> 356.7   16 -> 12.5     89.562   -> 73.053
                                     3      508     -> 511     17 ->  8.5     138.486  -> 75.649

OT_Tilted_364_200_Pixel_4021.cfg     OT Version 3.6.4 !! Work in progress !!
                                     Starting the same process as for 3.6.3, with updated δ
                                     Layer 1: δ = 3.9mm
                                     Layer 2, 3: δ = 3.4mm
                                     Since the first tilted ring -> last flat ring had a Δz = 43.58 mm and was already considered to be "close"
                                     and the δ reduction caused this number to further shrink (down to 43.33) , we adjusted zOverlap (16.0 → 15.5)
                                     to push back Δz to 43.62 (basically the same as 3.6.3). Layer 2 and 3 are not expected to be troublesome with a
                                     Δz(last flat → first tilted) = 47.1 and 46.8 respectively).
                                     Overall the change 3.6.3 → 3.6.4 had a positive effect on the coverage of tracks from the origin, with
                                     zErrorOuter = { 51.8, 76.0, 78.8 }: a slight reduction for layer 1 and an increase for layers 2 and 3


OT_Tilted_462_200_Pixel_4021.cfg     OT Version 4.6.2 <- like 3.6.2 but Avi-style
                                     Pixel version 4.0.2.1

OT463_200_IT4025.cfg                 OT Version 4.6.3 <- like 4.6.2 but with 5 endcap disks
                                     Pixel version 4.0.2.5

OT365_200_IT500.cfg                  OT Version 3.6.5
                                     Pixel version 5.0.0 : created series of tilted pixel layouts. Geometry is first jet and needs to be more optimized :)


OT_Tilted_362_200_Pixel_4023.cfg     OT Version 3.6.2
                                     Pixel version 4.0.2.3 <- like 4.0.2.1 but with 8 FPIX_1 disks and 4 FPIX_2 disks

OT365_200_IT4022.cfg                 OT Version 3.6.5  <- like 3.6.4 but with adjusted tilted ring positions and increased layer radii
                                     Inner tracker version 4.0.2.2

OT365_200_IT4024.cfg                 OT Version 3.6.5  <- like 3.6.4 but with adjusted tilted ring positions and increased layer radii
                                     Inner tracker version 4.0.2.4 <- like 4.0.2.3 but with smaller barrel radii, just like 4.0.2.2

OT365_200_IT4025.cfg                 OT Version 3.6.5
                                     Inner tracker version 4.0.2.5 <- based on 4.0.2.4 but with
                                     BPIX_4.3.0 adjusted radii: #rods reverted to 12,28,24,32
                                                                smaller overlap (0.6mm) in L2
                                                                bigger overlap (2.0mm) in L3,4
                                                                radii 29.000, 70.146, 117.753, 157.388
                                     and FPIX_2 shifted by 10 cm inwards
                                            disk z: 1750.0, 1985.43, 2250.83, 2550.0
                                     * Added the pixel support tube and service cylinders

OT365_200_IT4026.cfg                 OT Version 3.6.5
                                     Inner tracker version 4.0.2.6 <- based one 4.0.2.5 but with 7 FPIX disks

OT366_200_IT4025.cfg                 OT Version 3.6.6 <- like 365, but with adjusted TBPS flat part
                                     bigDelta = 11.9 (unchanged) smallDelta=3.5625 (from CML)
                                     New TEDD following Nick's indications:
                                     Outer Physical envelope = 1125 -- Considering 20.54 mm 2S FEH + 1.1 mm sensor margin => outerRadius 1103.36 (+8.36 mm w.r.t v3.5.1)
                                     Inner Physical envelope = 215/315 -- Considering 10.65 mm PS FEH + 1.45 mm sensor margin => 227/327 mm inner envelope for sensors
                                     Disk    z [mm]
                                     1_1(1)  1326.8
                                     1_2(2)  1550.0
                                     2_1(3)  1853.4
                                     2_2(4)  2216.2
                                     2_3(5)  2650.0
                                        bigDelta =30.7/2=15.35
                                        Module/ring            smallDelta:
                                        2S 4.0mm rings 11-12   17.10/2 =  8.55 mm
                                        2S 4.0mm ring  10      20.10/2 = 10.05 mm
                                        2S 1.8mm rings (11-15) 14.90/2 =  7.45 mm
                                        PS 4.0mm all rings     14.84/2 =  7.42 mm
                                     IT version 4.0.2.5
                                     
                       Geometries with timing layer              
OT400_200_IT4026.cfg   Lindsey Gray version
                       OT Version 4.0.0 <- like 3.6.5 but
                          with the last TEDD disk moved inwards by 10 cm and an additional timing layer
                          active pads are 1×3mm² (x×y)
                       IT Version 4.0.2.6 (standard)
OT466_200_IT4025.cfg   OT Version 4.6.6 <- like 3.6.6 but
                          with the last TEDD disk moved inwards by 10 cm and an additional timing layer
                          covering η 1.8 → 2.8
                          using 5×10 cm² active sensor with vertical orienation (minimal # modules)
                          active pads are 2×2mm² (x×y)
                       IT Version 4.0.2.5 (standard)
OT467_200_IT4025.cfg   OT Version 4.6.7 <- like 4.6.6 but with 1×4mm² (x×y)
OT468_200_IT4025.cfg   OT Version 4.6.8 <- like 4.6.6 but with 4×1mm² (x×y)
OT469_200_IT4025.cfg   Lindsey Gray timing layer plugged on TDR layout
                       OT Version 4.6.9 <- like 6.1.3 but
                          with the last TEDD disk moved inwards by 10 cm and an additional timing layer
                          active pads are 1×3mm² (x×y)
                       IT Version 4.0.2.5 (standard)

OT367_200_IT4025.cfg  OT Version 3.6.7 <- like version 3.6.6, but with small adjustments from Nick (WORK IN PROGRESS!!!)
                                            Outer Radius reduced from 1103 mm to 1100 mm
                                            bigDelta increased from 15.35 mm to 15.65 mm
                                            A suggested inwards shift R7 radius by 2.5mm and that of R5 by 0.5mm is NOT implemented here
                                            the disk stretch was implemented instead in order to obtain
                                            inner radius for ring 1 to be 229.3 mm
                                            inner radius for ring 3 to be 329.3 mm
                                            FIX: rOverlap is set to 0

OT567_200_IT4025.cfg  OT Version 5.6.7 <- like version 3.6.7, but leaving the two TEDD sections free to be different
                                            this allows to remove one ring more from TEDD2
                                            Requirement to reach 229 mm and 329 mm with innermost ring is obtained
                                            by setting different zError: in TEDD1 zError=139.5mm and TEDD2 zError=283mm

OT568_200_IT4025.cfg  OT Version 5.6.8 <- like version 5.6.7, but replacing one 2S ring with 2 PS rings
                                            Requirement to reach 229 mm and 329 mm with innermost ring is obtained
                                            by setting different zError: in TEDD1 zError=139.5mm and TEDD2 zError=176.83mm

OT577_200_IT4025.cfg  OT Version 5.7.7 <- like 5.6.7, but with overlap spread over first rings
OT578_200_IT4025.cfg  OT Version 5.7.8 <- like 5.6.8, but with overlap spread over first rings
OT377_200_IT4025.cfg  OT Version 3.7.7 <- like 3.6.7, but with overlap spread over first rings
OT601_200_IT4025.cfg  OT Version 6.0.1 <- starting from 5.7.8 we solve PS->2S clash by shifting the PS rings instead of adding additional
                                          smallDelta. This allows to reduce bigDelta to 14.15 mm
                                          10,7 PS rings in TEDD1,2
                                             5 2S rings in TEDD1&2
                                          phioverlap = -2
                                          gaps distributed with bigDelta. Plenty of room

OT602_200_IT4025.cfg                       9,7 PS rings in TEDD1,2
                                             6 2S rings in TEDD1&2
                                           overlap mainly through bigDelta (and 2S)
                                           phioverlap = -2

OT603_200_IT4025.cfg                      like 602, but overlap mainly through PS modules

OT611_200_IT4025.cfg                      TDR geometry

OT711_200_IT4025.cfg                      Like 6.1.1 TDR geometry but with paired-up layers:
                                             Z coordinates in TEDD: 1356.8 1440 1854.67 1937.67 2567 2650
                                             while keeping the radial positions as in 6.1.1
                                          TB2D radii taken from 4.6.3
                                          TBPS configuration as in 6.1.1, but reducing L2 radius to have 1 cm clearance (incl. hybrid) to L1


                                          
OT612_200_IT4025.cfg	Like 6.1.1 but with slightly larger PS modules

OT613_200_IT4025.cfg	Like 6.1.2 but fixing bigDelta according to Nick 2017-03-27
                         Zd=29.7mm or bigDelta=14.85
                         New geometry 6.1.3 with slight adjustment in TEDD1 and TEDD2: bigDeelta moved from 14.15mm to 14.85mm
                         this would cause a movement of 1.243mm and 1.351mm in innermost rings of TEDD1 and TEDD2 respectively
                         and add 4 modules in one ring of TEDD2.
                        Any ring movement was instead *avoided* by constraining the ring radii to those of 6.1.2, so that the geometry in the
                        xy plane is exactly the same


