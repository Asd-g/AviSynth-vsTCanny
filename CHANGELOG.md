##### 1.1.8:
    Fixed default sigma_vY when the clip has only one plane. (regression from 1.1.7)

##### 1.1.7:
    Changed the behavior of default sigma_vU.

##### 1.1.6:
    Fixed default sigma_U/V/vU/vV for RGB formats.
    Changed default sigma_vU/vV. Now they are inherited from sigmaU/V.

##### 1.1.5:
    Fixed the processing of planes for RGB formats.
    Properly clamped float mask to 0-1 range in mode=1. (VS plugin r14)

##### 1.1.4:
    Fixed the behavior when y/u/v=1.

##### 1.1.3:
    Fixed the uninitialized variables when the clip has only one plane.

##### 1.1.2:
    Fixed a bug when sigma=0 and the plane is not processed.

##### 1.1.1:
    Fixed the processing of clips with one plane.

##### 1.1.0:
    Changed chroma planes range from -0.5..0.5 to 0.0..1.0 (float clips). (VS plugin r13)
    Added AVX512 code. (VS plugin r13)
    Added Kroon, Kirsch and FDoG operatos. (VS plugin r13)
    Renamed `gmmax` parameter to `scale` and changed its default to 1.0. (VS plugin r13)
    Changed default sigma_vY from 1.5 to sigmaY.

##### 1.0.1:
    Fixed sigma for RGB clips.

##### 1.0.0:
    Port of the VapourSynth plugin TCanny r12.
