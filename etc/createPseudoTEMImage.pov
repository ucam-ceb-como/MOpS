//Creates pseudo-TEM image.
//
//Note that you must have installed POV-Ray (http://www.povray.org/download/)
//to use this script. For the binary tree models, if ptrack in the mops.inx
//file is set to true, MOpS will print out the coordinates and radius of each
//primary in an aggregate for each particle which is tracked at each time step.
//To create a more realistic TEM image, the povray image can be further
//processed using the Color to Alpha feature in GIMP
//(https://www.gimp.org/downloads/).
//
//Written by Edward Yapp ekyy2@cam.ac.uk

#include "colors.inc"
#declare tem_finish =
    texture {
        finish {ambient 0.4 diffuse 0.7 phong 0.1}
        pigment {color Gray90 transmit 0.38}
    }
#declare backSphere = 
    texture {
        pigment {color Gray90 transmit 0.95}
    }
background {color White}
#declare RNG = seed (61);

// Should the circles of collision diameter be plotted?
#declare plotCol=0;

// Should the scale bar be drawn? (200 nm)
#declare scaleBar=1;

// Length and width of image in nm
#declare W = 2000;
#declare H = 2000;
#declare dLight = 10*max(W, H);

// Declare base filename
#declare fbase = "calc/0.8/wu-tem";
#declare time = "0.8s";

// Declare number of particles to draw
#declare nump = 25;
#declare i = 1;

// Loop over particles to plot them all
#while (i < nump)
    #declare fname = concat(fbase, "(", time, ", ", str(i,0,0), ").pov")
    #include fname
    #declare d  = ParticleDiameter;
    #declare sp = blob {MyParticle};
    
    blob {
        sp
        texture {tem_finish}
        rotate y*(360*rand(RNG))
        translate <rand(RNG)*W, 0, rand(RNG)*H>
    }
    #if (plotCol)
    sphere {c0, d/2 texture {backSphere}}
    #end
    
    #declare i = i + 1;
#end

// DRAW SCALE BAR
#if (scaleBar)
    box {
    <W/2-100, 0, H/2+10>
    <W/2+100, 0, H/2-10>
    texture {tem_finish}
    }
#end

// DEFINE CAMERA AND REFERENCES
light_source {<W, -dLight, H> White}
#declare Y  = 1.5*W*tan(45.0/2.0);
#declare th = 2*asin(W/(2*Y));
camera {
    look_at <W/2, 0, H/2>
    location <W/2, Y, H/2>
    right x
    up (W/H)*z
}
