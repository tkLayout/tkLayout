Pixel models

Baseline_tilted_200_Pixel_1_1_1.cfg  Square pixels everywhere (small L1-2 and R1-2)
Baseline_tilted_200_Pixel_1_1_2.cfg  Like 1.1.1, but without DC/DC converters
Baseline_tilted_200_Pixel_1_2_1.cfg  Like 1.1.1, but Rectangular pixels in barrel
Baseline_tilted_200_Pixel_1_3_1.cfg  Like 1.1.1, but Rectangular everywhere
Baseline_tilted_200_Pixel_1_4_1.cfg  Like 1.3.1 But a bit shorter barrel
Baseline_tilted_200_Pixel_1_5_1.cfg  Like 1.3.1 But a lot shorter barrel
Baseline_tilted_200_Pixel_1_5_2.cfg  Like 1.5.1, but without DC/DC converters

Baseline_tilted_200_Pixel_3_5_1.cfg  2-rings step Like, TEDD and TEPX disks moved to allow insertion -- small pixels everywhere

Baseline_tilted_200_Pixel_4_0_0.cfg  Like 3_5_1, but with 22x16.4 mm^2 chips (good for tilt?)
                                     BPIX=25x100 TFPX=50x50 TEPX=25x100

Baseline_tilted_200_Pixel_4_0_1.cfg  Like 4_0_0, but BPIX=50x50
                                     i.e. BPIX=50x50 TFPX=50x50 TEPX=25x100

Baseline_tilted_200_Pixel_4_0_2.cfg  Like 4_0_0, but TFPX=25x100
                                     i.e. BPIX=25x100 TFPX=25x100 TEPX=25x100

Baseline_tilted_200_Pixel_4_0_3.cfg  Like 4_0_2, but
                                     TEPX=50x50

Baseline_tilted_200_Pixel_4_0_2_1.cfg Like 4_0_2, but with 7 TFPX disks and 4 TEPX disks

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

OT362_200_IT4121.cfg		     OT Version 3.6.2
				     Pixel version as in 4.0.2.1, but with 50x50 pixels instead of 25x100

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

OT_Tilted_362_200_Pixel_4023.cfg     OT Version 3.6.2
                                     Pixel version 4.0.2.3 <- like 4.0.2.1 but with 8 TFPX disks and 4 TEPX disks

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
                                     and TEPX shifted by 10 cm inwards
                                            disk z: 1750.0, 1985.43, 2250.83, 2550.0
                                     25x100 everywhere.
                                     Added the pixel support tube and service cylinders.

OT365_200_IT4026.cfg                 OT Version 3.6.5
                                     Inner tracker version 4.0.2.6 <- based one 4.0.2.5 but with 7 FPIX disks

OT366_200_IT4025.cfg                 OT Version 3.6.6 <- like 365, but with adjusted TBPS flat part.
                                     bigDelta = 11.9 (unchanged) smallDelta=3.5625 (from CML).
                                     12-th ring in Layer 1.
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

OT611_200_IT4025.cfg                      Adjusted rings radii in TEDD.

OT711_200_IT4025.cfg                      Like 6.1.1 TDR geometry but with paired-up layers:
                                             Z coordinates in TEDD: 1356.8 1440 1854.67 1937.67 2567 2650
                                             while keeping the radial positions as in 6.1.1
                                          TB2D radii taken from 4.6.3
                                          TBPS configuration as in 6.1.1, but reducing L2 radius to have 1 cm clearance (incl. hybrid) to L1


                                          
OT612_200_IT4025.cfg	Like 6.1.1 but with slightly larger PS modules


======================   TDR LAYOUT   ====================

OT613_200_IT4025.cfg	   Like 6.1.2 but fixing bigDelta according to Nick 2017-03-27
                         Zd=29.7mm or bigDelta=14.85
                         New geometry 6.1.3 with slight adjustment in TEDD1 and TEDD2: bigDelta moved from 14.15mm to 14.85mm
                         this would cause a movement of 1.243mm and 1.351mm in innermost rings of TEDD1 and TEDD2 respectively
                         and add 4 modules in one ring of TEDD2.
                        Any ring movement was instead *avoided* by constraining the ring radii to those of 6.1.2, so that the geometry in the
                        xy plane is exactly the same
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


============   TIMING BARREL LAYER STUDIES   =============
OT613_200_IT4025.cfg  Reference "standard" 78 ladders
OT624_200_IT4025.cfg  -15 mm outer radius  78 ladders
OT625_200_IT4025.cfg  -30 mm outer radius  76 ladders (-2)
OT626_200_IT4025.cfg  -30 mm outer radius  78 ladders
OT627_200_IT4025.cfg  -69 mm outer radius  72 ladders (-6)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


===========   POST-TDR OUTER TRACKER STUDIES   ===========                    
OT614_200_IT4025.cfg	  OT 614: Like 6.1.3 but updates from Nick 2017-11-07 in TEDD.
                        'Identical dees' design (simplify cooling design, huge saving :) ):
                          - Rotate in Phi all modules in all disks in TEDD by + half a module.
                        PS module design has been updated, hence:
                          - Increase bigDelta, from 14.85 mm to 15.075 mm (Zd goes from 29.7 mm to 30.15 mm), in all disks in TEDD.
                          - Decrease smallDelta, from 7.42 mm to 7.375 mm (ZmPS4.0 goes from 14.84mm to 14.75mm), in PS4.0 rings.
                          
OT615_200_IT404.cfg	   Diff with OT614 is in TBPS:
                           L1 small delta: 3.5625 mm -> 3.5475 mm
                           L2 small delta: 3.47 mm -> 3.0475 mm
                           L3 small delta: 3.47 mm -> 3.5844 mm
                           All layers: flipped bigParity.
                           
OT616_200_IT404.cfg	   Diff with OT615:
                       Reduced outermost radius to leave space for BTL. Increased innermost radius for IT insertion.
                         - TEDD:
                           bigDelta: 15.075 mm -> 15.755 mm (needed more disk separation).
                           TEDD 1, inner rings: +7.041 mm (IT insertion) +0.51 mm (margin with dee edge).
                           TEDD 2, inner rings: +2 mm (IT insertion) +0.51 mm (margin with dee edge).
                           TEDD 1 and 2, outer rings: -27 mm (BTL) +0.41 mm (margin with dee edge).
                           TEDD 2: -4 modules in Ring 7, -4 modules in Ring 14.
                           Adjusted intermediate radii.
                         - TB2S:
                             L3: radius -27 mm (OT envelope shrink) +2 mm (smaller no-go zone between TB2S and BTL). numRods: -2 rods.
                         - TBPS:
                             L1: +2 mm in last 5 rings radii.
                             
OT616_IT613.cfg	       Diff with OT616_200:
                       Sensor thickness: 200 um -> 290 um + Add 30 um deep-diffused Si + Review Si in inactive edges and MPA.
                       Add 164 um to s-sensor in PS module.
                       
OT617_IT615.cfg	       Diff with OT616:
                       All TEDD: update smallDeltas to latest info.
                          * PS 4 mm: 7.375 mm -> 7.550 mm.
                          * 2S 4 mm:  8.550 mm -> 8.595 mm.
                          * 2S 1.8 mm: 7.450 mm -> 7.495 mm.
                       Improve transition region between TB2S and TEDD by removing readout hybrid in TEDD1, R15.
                          * Special modules in TEDD1, R15, with halved number of strips / module. WARNING: MB NOT UPDATED!! SHOULD UPDATE MODULE MB.
                          * Outer radius increased by 20.47 mm : R15 sensors centers: 1023.16 mm -> 1043.63 mm.
                          * Radii in TEDD1, R12, R13, R14 adjusted accordingly.
                       
OT618_IT615.cfg	       Diff with OT617:
                       Also improve transition between TB2S/TBPS and TEDD by removing readout hybrid in TEDD2, R15.
                          * Special modules in TEDD2, R15, with halved number of strips / module. WARNING: MB NOT UPDATED!! SHOULD UPDATE MODULE MB.
                          * Outer radius increased by 20.47 mm : R15 sensors centers: 1023.16 mm -> 1043.63 mm.
                          * Radii in TEDD2, R12, R13, R14 adjusted accordingly.
                          
OT800_IT615.cfg	       Based from Outer Tracker version 616.
		       All TEDD: 
		       	  Update smallDeltas to latest info. Identical to the smallDelta updates in OT617. 
		       	  Excluding all the other updates in OT617 though (related to special 2S module in R15 and radii adjustments).
		       TBPS:
		          Adjustments in Layer 1 to enable a safe IT installation, while respecting needed spacing with TBPS Layer 2 (supports and services).
		       	  Layer 1, Ring 12-16: radius +0.5 mm, tiltAngle 74 -> 72 deg.

OT800_IT630.cfg        Based on OT800_IT615. Very small pixels for phase3: TBPX L1, TBPX L2, TFPX R1 have a 15µm × 60µm size

OT800_IT631.cfg        Based on OT800_IT630. But small pixels in Layer 1,2 and Ring 1,2 

OT800_IT632.cfg        Based on OT800_IT631 (small pixels), but with no beam pipe material
		       	  
OT801_IT701.cfg	       Based from Outer Tracker version 800.
                       TB2S: inter-ladder radial spacing increased by 1 mm (smallDelta: 2.25 mm -> 2.75 mm). Z positions recomputed accordingly.		       	  
		                                          
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                          


===========   POST-TDR INNER TRACKER STUDIES   ===========                           
OT613_200_IT4125.cfg	  Like OT613_200_IT4025, but with 50x50 pixels instead of 25x100.

OT613_200_IT404.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.0.4:
                                      - geometry similar to IT4.0.2.5, but with radii of barrel layers reduced (same as flat part of IT 5.0.1).
                                      - module type 25x100 everywhere.
                                      
OT613_200_IT405.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.0.5:
                                      - geometry same as IT4.0.4
                                      - module type 50x50 everywhere.
                                      
OT613_200_IT406.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.0.6:
                                      - geometry same as IT4.0.4
                                      - 25x100 in 1x2 modules, 50x50 in 2x2 modules.                                  

OT613_200_IT407.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.0.7:
                                      - geometry same as IT4.0.4
                                      - 50x50 in 1x2 modules, 25x100 in 2x2 modules.   
                                      
OT613_200_IT408.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.0.8:  IT4.0.2.5 with 1 disk less in FPX1 and 1 disk less in FPX2.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                           
     
                                     
=========   LARGE PIXELS IN 2X2 MODULES STUDIES   =========                                                                
OT613_200_IT420.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.2.0:
                                      - geometry same as IT4.0.4
                                      - 25x100 in 1x2 modules, and 50x200 in 2x2 modules.
                                      
OT613_200_IT421.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.2.1:
                                      - geometry same as IT4.0.4
                                      - 25x100 in 1x2 modules, and 100x100 in 2x2 modules. 
                                      
OT613_200_IT422.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.2.2:
                                      - geometry same as IT4.0.4
                                      - 50x50 in 1x2 modules, and 50x200 in 2x2 modules.
                                      
OT613_200_IT423.cfg                  OT Version 6.1.3
                                     Inner Tracker version 4.2.3:
                                      - geometry same as IT4.0.4
                                      - 50x50 in 1x2 modules, and 100x100 in 2x2 modules.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     


===============    ATLAS CHIP SIZE STUDY   ===============      
OT614_200_IT430.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.3.0: smaller chip: 16.4 mm x 20 mm chip, instead of 16.4 mm x 22 mm chip (same length as Atlas).
                                     Based from Inner Tracker version 4.0.4.
                                      - TBPX: shorter of 18 mm.
                                      - TFPX: First disk moved inwards of 18 mm in Z. Other disks Z positions adjusted accordingly to have the same TFPX max Z. Radii: R1 min and R4 max identical, the radii in between are adjusted (coverage). R2 center: +1 mm, R3 center: -2 mm. 
                                      - TEPX: 4 modules added in R4. Radii: R1 min and R5 max identical, the radii in between are adjusted (coverage). R2 center: -3 mm, R3 center: -2.5 mm, R4 center: -6 mm.  
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


=================   1X1 MODULES STUDIES   ================      
OT614_200_IT440.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.4.0: 1x1 modules in TEPX Ring 1 (Swiss layout).
                                     2x2 modules in TEPX Ring 2.
                                     Based from Inner Tracker version 4.0.4.
                                     TEPX: 
                                     - Ring 1: Rhigh: 108 mm -> 86 mm,  numModules: 40 -> 36.
                                     - Ring 2: Rhigh: 149 mm -> 130 mm, numModules: 56 -> 28.
                                     - Ring 3: Rhigh: 188.5 mm -> 171 mm.
                                     - Ring 4: Rhigh: 232 mm -> 215 mm.
                                     
OT614_200_IT441.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.4.1: Also 1x1 modules in TFPX Ring 1 (Correction by Duccio on Swiss layout).
                                     2x2 modules in TEPX Ring 2.
                                     Based from Inner Tracker version 4.0.4.
                                     TFPX: 
                                     - Ring 1: Rhigh: 73.2 mm -> 51 mm.
                                     TEPX: 
                                     - Ring 1: Rhigh: 108 mm -> 86 mm,  numModules: 40 -> 32.
                                     - Ring 2: Rhigh: 149 mm -> 127 mm, numModules: 56 -> 28.
                                     - Ring 3: Rhigh: 188.5 mm -> 169 mm.
                                     - Ring 4: Rhigh: 232 mm -> 213 mm.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


=================   LUMINOSITY STUDIES   =================     
OT614_200_IT450.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.5.0: Tune radii in TEPX to have radially distributed luminosity measurements.
                                     Based from Inner Tracker version 4.0.4.
                                     TEPX:
                                     - Ring 2: Rhigh: 149 mm -> 145 mm.
                                     - Ring 3: Rhigh: 188.5 mm -> 182 mm.
                                     - Ring 4: Rhigh: 232 mm -> 219 mm.
                                     The resulting radial spacings are, using (Ring (i) Rhigh)  - (Ring (i+1) Rmin):
                                     - i=1: 7.2 mm
                                     - i=2: 7.2 mm
                                     - i=3: 7.2 mm
                                     - i=4: 9.2 mm
                                     
OT614_200_IT451.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.5.1: Tune radii and numModules in TEPX to try to equilibrate luminosity measurements.
                                     Based from Inner Tracker version 4.5.0.
                                     TEPX:
                                     - Ring 2: Rhigh: 145 mm -> 146 mm.
                                     - Ring 3: Rhigh: 182 mm -> 183 mm.
                                     - Ring 4: Rhigh: 219 mm -> 220 mm, numModules: 40 -> 44.
                                     The resulting radial spacings are, using (Ring (i) Rhigh)  - (Ring (i+1) Rmin):
                                     - i=1: 6.2 mm
                                     - i=2: 7.2 mm
                                     - i=3: 7.2 mm
                                     - i=4: 10.2 mm
                                     
OT614_200_IT452.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.5.2: Tune radii in TEPX, perfectly equilibrated among rings thanks to the precise counts of fraction of tracks with at least 3 hits.
                                     Based from Inner Tracker version 4.5.0.
                                     TEPX:                                     
                                     * Ring 2: Rhigh: 145 mm -> 149 mm.
                                     * Ring 3: Rhigh: 182 mm -> 186 mm.
                                     * Ring 4: Rhigh: 219 mm -> 223 mm, numModules: 40 -> 44.
                                     The resulting radial spacings are, using (Ring (i) Rhigh)  - (Ring (i+1) Rmin):
                                     * i=1: 3.2 mm
                                     * i=2: 7.2 mm
                                     * i=3: 7.2 mm
                                     * i=4: 13.2 mm                                     
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


=================   MATERIAL BUDGET STUDIES   =================     
OT614_200_IT460.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.6.0: TBPX conversion station placed on top of TBPX.
                                     
OT614_200_IT461.cfg                  OT Version 6.1.4
                                     Inner Tracker version 4.6.1: TBPX conversion stations placed on services cylinders (Z positions from Yadira).                                                                       
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
                                                                                             

============   TILTED INNER TRACKER STUDIES   ============                 
OT613_200_IT500.cfg                  OT Version 6.1.3
                                     Inner Tracker version 5.0.0 : tilted Inner Tracker. Head of series.
                                     First-jet optimization :)
                                     1-chip modules everywhere in the tilted Barrel.
                                     TFPX : 2 disks were removed with respect to TDR version.
                                     TEPX : same as TDR version.
                                     50x50 in BPIX, 25x100 in FPIX.
                                     
OT613_200_IT501.cfg                  OT Version 6.1.3
                                     Inner Tracker version 5.0.1 : tilted Inner Tracker. 
                                     Barrel : Modules are bigger than in IT500, hence significant reduction of number of modules. Barrel modules: 3748 (IT500) -> 1844 (IT501).
                                         (a) longer modules with two chips in the flat part of layers 1-2 (dimensions 16.4 x 44.2)
                                         (b) still the same 1-chip modules in the tilted part of layers 1-2 (dimensions 16.4 x 22, rotated by 90 degrees)
                                         (c) larger modules with 2x2 chips in the flat part of layers 3-4 (dimensions 33 x 44.2)
                                         (d) 1x2 modules in the tilted part of layers 3-4 (dimensions 16.4 x 44.2, rotated by 90 degrees) 
                                     TFPX : 2 disks were removed with respect to TDR version.
                                     TEPX : same as TDR version.
                                     50x50 in BPIX, 25x100 in FPIX.
                                     
OT613_200_IT502.cfg                  OT Version 6.1.3
                                     Inner Tracker version 5.0.2: 
                                      - geometry same as IT5.0.1
                                      - module type 50x50 everywhere.
                                     
OT613_200_IT503.cfg                  OT Version 6.1.3
                                     Inner Tracker version 5.0.3: 
                                      - geometry same as IT5.0.1
                                      - 50x50 in 1x1 and 1x2 modules, 25x100 in 2x2 modules.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    


============   SKEWED INNER TRACKER STUDIES   ============                 
OT614_200_IT600.cfg                  OT Version 6.1.4
                                     Inner Tracker version 6.0.0 : skewed Inner Tracker. Head of series.
                                     The shift of the edge of a skewed module is set to 5 mm.
                                     The phi positions are automatically computed, so that the angular overlap around the (X=0) plane
                                     is 2.x the angular overlap between 2 standard consecutive rods.
                                     
OT614_200_IT601.cfg                  OT Version 6.1.4                                     
                                     Based from Inner Tracker version 6.0.0.
                                     Mid-radii in TBPX:                                     
                                     * Layer 3: 100 mm -> 102 mm.
                                     * Layer 4: 140 mm -> 143 mm.
                                     
OT614_200_IT602.cfg                  OT Version 6.1.4                                     
                                     Based from Inner Tracker version 6.0.1.
                                     bigDelta in TBPX:      
                                     * Layer 2: 1.5 mm -> 2.5 mm.                                                                    
                                     * Layer 3: 1.5 mm -> 2.5 mm.  
                                     * Layer 4: 1.5 mm -> 2.5 mm.  
                                     Question: so Layer 1 would need to be lifted up?   
                                     
OT614_200_IT610.cfg                  OT Version 6.1.4                                     
                                     Based from TBPX version 6.0.2., and TFPX + TEPX version 4.5.1.
                                     Everywhere:
                                     ROC size: 16.4 mm x 22 mm -> 16.8 mm x 21.6 mm.
                                     Spacing between ROCs: 0.2 mm -> 0.25 mm.
                                     TBPX:      
                                     * Layer 1: layerRho: 29 mm -> 30 mm. bigDelta 1.5 mm -> 2.5 mm.                                                                    
                                     * Layer 2: layerRho: 60 mm -> 61.5 mm.    
                                     * Layer 3: layerRho: 102 mm -> 104.5 mm.  
                                     * Layer 4: layerRho: 143 mm -> 146.5 mm.  
                                     Spacing in Z: 0.2 mm -> 1.3 mm.
                                     TFPX:
                                     * Ring 1: Rmin + 2 mm (pending change).
                                     * Ring 2: Rcenter: +0.2 mm. 
                                     * Ring 4: Rcenter: +0.35 mm (keep same Rhigh). 
                                     TEPX:
                                     * Ring 1: Rcenter: -0.375 mm (keep same Rmin). numModules: 40 -> 36.
                                     * Ring 2: Rcenter: 126.9 mm -> 116.25 mm.
                                     * Ring 3: Rcenter: 166.4 mm -> 155.75 mm.
                                     * Ring 3: Rcenter: 209.9 mm -> 196.25 mm.
                                     * Ring 4: Rcenter: +0.36 mm (keep same Rhigh).       
                                     
OT614_200_IT611.cfg                  OT Version 6.1.4                                     
                                     Based from Inner Tracker version 6.1.0.
                                     TEPX: 2x2 modules in Ring 1 and Ring 2.
                                     Materials: adapted correspondingly.
                                     NB: IT cabling still follows IT404 cabling map, will need slight update (+4 modules in R4).
                                     WARNING: Rings radii to be tuned!! 
                                     Rmin set to 62.9 mm, is that fine? Rmax < 254.52 mm, is that the constrainst?  
                                     What is the Rmax constrainst on TFPX???  
                                     
OT614_200_IT612.cfg                  OT Version 6.1.4                                     
                                     Based from Inner Tracker version 6.1.1.
                                     bigDelta: 4 mm -> 2 mm, smallDelta: 2 mm -> 4 mm. 
                                     
OT614_200_IT613.cfg                  OT Version 6.1.4                                     
                                     Based from Inner Tracker version 6.1.2.
                                     Stretched TEPX in Z (there will be no dedicated lumi device, so can go up to bulkhead).
                                     * Disk 1: same Z.
                                     * Disk 2: 1985.43 mm -> 2009.59 mm   
                                     * Disk 3: 2250.83 mm -> 2307.69 mm       
                                     * Disk 4: 2550.00 mm -> 2650.00 mm (set same Z as last TEDD disk's meanZ for now).                                                                                                                                                                            
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   


=================   INSTALLATION STUDIES   =================     
OT616_IT614.cfg                      OT Version 6.1.6
                                     Based from Inner Tracker version 6.1.3.
                                     TFPX, Disk 6 : -15 mm. Disk 7: + 15mm. All TBPX conversion stations placed between Disk 6 and Disk 7.
                                     
OT616_IT615.cfg                      OT Version 6.1.6
                                     Based from Inner Tracker version 6.1.4.
                                     TFPX: Stretched double-disk system, data from Yadira. 
                                     smallDelta: 2 -> 2.75 mm. bigDelta: 4 -> 6.25 mm.
                                     Sensors on both faces of same dee: inter sensor centers distance (in Z) is now 5.5 mm.
                                     Back sensors of one dee VS front sensors of next dee: inter sensor centers distance (in Z) is now 7 mm.
                                     
OT616_IT616.cfg                      OT Version 6.1.6
                                     Based from Inner Tracker version 6.1.5.
                                     TFPX: Tested shift in Z of +25 mm on all double-disks.                                     
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  


================   LOCAL RESOLUTION STUDIES   ============                 
OT616_200_IT620.cfg                  OT Version 6.1.6
                                     Based from Inner Tracker version 6.1.3.
                                     Zebra layout: Alternation of layers with 25x100 and 50x50 pixel aspect ratios.
                                     Layouts fully with 25x100 or 50x50 were studied, but not one that would include a mix of the 2!
                                     
OT618_200_IT621.cfg                  OT Version 6.1.8
                                     Based from Inner Tracker version 6.1.5.
                                     Same as IT 6.1.5., but with 50x50 everywhere (instead of 25x100).                                   
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


=================   3D SENSORS STUDIES   =================
OT616_IT699.cfg                      OT Version 6.1.6
                                     Based from Inner Tracker version 6.1.5.
                                     Reduced holes in TBPX L1 and L2: 1.3 mm -> 0.6 mm, to account for future 3D sensors spacing.
 
OT616_IT700.cfg                      OT Version 6.1.6
                                     Based from Inner Tracker version 6.9.9.
                                     Same geometry + materials, but 3D sensor type in TBPX L1 + TBPX L2 + TFPX R1 (+ slightly different inactive Si).
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


================   INSTALLATION STUDIES   ================
OT800_IT701.cfg                      OT Version 8.0.0
                                     Based from Inner Tracker version 7.0.0.
                                     TFPX: Change of ring paradigm: modules of the same ring are now on 2 different dees, instead of being on both sides of the same dee.
                                     This model is already used in TEPX. 
                                     This makes the dee model mecanically feasible (cut at X~0), and better equilibrates services distribution among disks.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                     


=============   SENSOR TECHNOLOGY STUDIES   =============
OT800_IT702.cfg                      OT Version 8.0.0
                                     Based on IT version 7.0.0.
                                     TPBX L1 with 3D sensors, planar 25x100 um2 pixels in TBPX L2,L3,L4, TFPX and TEPX

OT800_IT703.cfg                      OT Version 8.0.0
                                     Based on IT version 7.0.0.
                                     TPBX L1 with 3D sensors, planar 25x100 um2 pixels in TBPX L2,L3,L4; planar 50x50um2 pixels in TFPX and TEPX
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                     

=============   BRICKED SENSOR STUDIES   =============
OT800_IT800.cfg                      OT Version 8.0.0
                                     Based on IT version 7.0.0.
                                     TPBX L1 with 3D sensors, TBPX L2 bricked in central rod; TBPX L3+L4 bricked in central 3 rods, planar 25x100 um2 pixels elsewhere in TBPX
                                     Bricked pixels everywhere in TEPX and TFPX

OT800_IT801.cfg                      OT Version 8.0.0
                                     Based on IT version 8.0.0.
                                     TPBX L1 with 3D sensors, TBPX L2 bricked in central rod; TBPX L3+L4 bricked in central 3 rods, planar 25x100 um2 pixels elsewhere in TBPX, TFPX and TEPX

OT800_IT802.cfg                      OT Version 8.0.0
                                     Based on IT version 8.0.0.
                                     TPBX L1 with 3D sensors, TBPX L2,L3 bricked in central rod; TBPX L4 bricked in central 3 rods, planar 25x100 um2 pixels elsewhere in TBPX and disks 1-4 of TFPX. 
                                     Bricked pixels in TFPX disks 5-8 and in all of TEPX. 

OT800_IT803.cfg                      OT Version 8.0.0
                                     Based on IT version 8.0.0.
                                     TPBX L1 with 3D sensors, TBPX L2,L3 bricked in central rod; TBPX L4 bricked in central 3 rods, planar 25x100 um2 pixels elsewhere in TBPX.
                                     50x50 um2 pixels in TFPX and TEPX.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                     


OT806_IT743.cfg                      Just like OT806_IT742, but with moduleSubType assigned to the pixel modules

