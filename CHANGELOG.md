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
