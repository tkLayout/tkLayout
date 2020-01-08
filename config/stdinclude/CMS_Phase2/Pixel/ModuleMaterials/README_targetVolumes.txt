      //                                                   (FOR XML EXPORT ONLY)
      //                                           MATERIAL ASSIGNMENT IN TARGET VOLUMES
      //                                                      PIXEL MODULE
      //
      //  Top View 
      //
      //    --------------------       
      //    |        (7)       |
      //    --------------------        y
      //    |   |          |   |        ^
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |  Sensor  |   |        |
      //    |(6)|    (2)   |(5)|        |
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |          |   |        +----> x
      //    --------------------    
      //    |        (8)       |                        
      //    --------------------                        z
      //                                                ^
      //  Side View                                     |
      //         ================          Hybrid  (1)  +----> x
      //     (6) ---------------- (5)      Sensor  (2)
      //     =========================     Chip    (3)
      //
      // Chip(3) volume can contain Bumps and any other material (supports, etc) for simplification.
      //
      // Volumes (5), (6), (7), and (8) are volumes of same thickness as the sensor, located around the sensor.
      // They are used for Si inactive areas.
      
      // Combinations
      // Volume (4) is volume (5) + volume (6) + volume (7) + volume (8).
      // Volume (0) does not target anything.
