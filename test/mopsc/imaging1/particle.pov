#include "colors.inc"
#include "bintree-tem(0.1s, 2).pov"

#declare tem_finish =
    texture {
        finish {ambient 0.25 diffuse 0.7 phong 0.08}
        pigment {color Gray70 transmit 0.38}
    }
background {color White}

blob { 
    MyParticle
    texture {tem_finish}
}

camera {
    look_at <0, 0, 0>
    location <0, ParticleDiameter*1.5, 0>
}
