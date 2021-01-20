#
# Author:    Eric J. Bringley (eb656) and Casper Lindberg (csl37)
#
# Purpose:   Draws a single TEM frame from a video file
#

import matplotlib.pyplot as plot
import numpy as np
import pandas as pd 
import vtk 
import sys

# Turn on debugging output 
debug=False

def read_coords(file_name):

    frame = list()
    data = pd.read_csv(file_name)

    # counters
    i=0 # primary particle counter
    j=0 # frame counter

    frameDict = {}
    particleList = list()

    # Add frame Dict to frame lists:
    frame.append(frameDict)
    frameDict['Time'] = data['Time (s)'][0]
    frameDict['coords'] = particleList

    # get number of rows
    m = len(data)

    #loop over all rows in file:
    for k in range(m):
        # New particle time:
        if(data['Time (s)'][k] != frameDict['Time']):
            i=0
            j = j+1
            frameDict={}
            particleList = list()
            frameDict['Time'] = data['Time (s)'][k]
            frameDict['coords'] = particleList
            frame.append(frameDict)
        # create particle list
        # coords [ x y z r ]
        coords = [ data['x (m)'][k], data['y (m)'][k], data['z (m)'][k], data['r (m)'][k] ]
        particleList.append(coords)

        # get frame position and orientation 
        # the frame position are the position coordinates of the tracked primary
        # the frame orientation is tracked using two orthogonal unit
        # vectors: the frame x and z axes

        vec_1 = [
            data['orient-x_x (m)'][k],
            data['orient-x_y (m)'][k],
            data['orient-x_z (m)'][k]            
        ]

        vec_2 = [
            data['orient-z_x (m)'][k],
            data['orient-z_y (m)'][k],
            data['orient-z_z (m)'][k]            
        ]
        #  check if tracked primary has non-zero z vector 
        if ((vec_1[0]*vec_1[0]+ vec_1[1]*vec_1[1]+ vec_1[2]*vec_1[2]) != 0):  
           # tracked primary position
           frameDict['position'] = coords
           # tracked frame orientation vectors
           frameDict['orient_x'] = vec_1    # x vector
           frameDict['orient_z'] = vec_2    # z vector
        if(debug):
            print("Read in primary data: ")
            print("Time: {}".format(frameDict['Time']))
            print("Coordinates")
            print(coords)
            print("orientation_x")
            print(vec_1)
            print("orientation_z")
            print(vec_2)
    if(debug):        
        print("END OF READ_COORDS")
    # for f in frame:
    #     print(f.keys())
    return(frame)

def translate_coords(position1, position2):
    new_coords = []
    new_coords.append( position1[0] - position2[0])
    new_coords.append( position1[1] - position2[1])
    new_coords.append( position1[2] - position2[2])
    return new_coords

def rotate_coords(position1, vec_1, vec_2):
    # math is in radians
    #calculate angles for rotation of frame z-axis
    len_z = np.sqrt(vec_2[0]**2+vec_2[1]**2+vec_2[2]**2)#
    if(vec_2[0]!=0):    #no x-component
        phi1 = np.arctan(vec_2[1]/vec_2[0])# 
        if(vec_2[0] > 0 and vec_2[1] < 0):
            phi1 = phi1 + 2*np.pi#
        if(vec_2[0] < 0):
            phi1 = phi1 + np.pi#
    else:
        phi1 = np.pi/2. 
    theta1 = np.arccos(vec_2[2]/len_z)#
    if(debug):
        print(phi1,theta1)
    #construct rotatation matrix
    rot1 = np.empty([3, 3])
    rot1[0,0] = np.cos(theta1)*np.cos(phi1)#
    rot1[0,1] = np.cos(theta1)*np.sin(phi1)#
    rot1[0,2] = -np.sin(theta1)#
    rot1[1,0] = -np.sin(phi1)#
    rot1[1,1] = np.cos(phi1)#
    rot1[1,2] = 0#
    rot1[2,0] = np.sin(theta1)*np.cos(phi1)#
    rot1[2,1] = np.sin(theta1)*np.sin(phi1)#
    rot1[2,2] = np.cos(theta1)#
    
    #rotate position1
    coords1 = np.matmul(rot1,np.transpose(position1))#
    if(debug):
        print(coords1)
    #rotate x orientation vector
    xorient = np.matmul(rot1,np.transpose(vec_1))#
    
    ##calculate rotation about z-axis
    len_x = np.sqrt(xorient[0]**2+xorient[1]**2+xorient[2]**2)#
    theta2 = np.arccos(xorient[0]/len_x)#
    if(xorient[1]<0): #get correct range for  rotation
        theta2 = 2*np.pi - theta2#
    
    #construct rotation matrix#
    rot2 = np.empty([3, 3])

    rot2[0,0] = np.cos(theta2)#
    rot2[0,1] = np.sin(theta2)#
    rot2[0,2] = 0#
    rot2[1,0] = -np.sin(theta2)#
    rot2[1,1] = np.cos(theta2)#
    rot2[1,2] = 0#
    rot2[2,0] = 0#
    rot2[2,1] = 0#
    rot2[2,2] = 1#
    
    ##rotate coords
    coords2 = np.matmul(rot2,coords1)#
    new_coords = np.transpose(coords2)#
    return new_coords

def reorientate(frames): 

    # reorientate all frames in frames
    new_frames = list()
    #check for frame_range (optional)
    #in not supplied reorientate all coordinates

    n_frames = len(frames)
    if(debug):
        print(n_frames)
    frame_range = np.arange(0,n_frames,1)
    for j in frame_range:     #loops over frames
        
        coords = frames[j]['coords']
        new_coords = np.empty_like(coords)
        #loops over primaries
        frameDict = {}
        new_frames.append(frameDict)
        n_prim = len(coords)
        frameDict['Time'] = frames[j]['Time']
        position = frames[j]['position']
        if(debug):
            print("Description of frame: {}".format(j))
            print("Time: {}".format(frames[j]['Time']))
            print("Coordinates")
            print(coords)
            print("orientation_x")
            print(frames[j]['orient_x'])
            print("orientation_z")
            print(frames[j]['orient_z'])
        for i in range(n_prim): 
            if(debug):
                print("Primary {}".format(i))
            #if this is a primary, (radius non-zero)
            if(coords[i][3] != 0 ):
                
                ##translate coordinates
                ##this recentres the frame on the tracked primary
                x,y,z  = translate_coords(coords[i],position)
                if(debug):
                    print(x,y,z)
                ##rotate coordinates after translation
                ##this returns the frame to the correct orientation 
                ##defined by the points on the x- and z-axes.
                x,y,z = rotate_coords([x,y,z], frames[j]['orient_x'], frames[j]['orient_z'])
                if(debug):
                    print(x,y,z)
                new_coords[i]= [x,y,z,coords[i][3]]
        ##copy remaining data
        frameDict['coords'] = new_coords
        ##We haven't copied the frame orientation and position
        ##information. This is no longer needed.
    return new_frames


def draw_particle(ren, coordinates,solid=0.5,ambient = 0.6, diffuse = 0.3, specular = 0.0,opacity=1,color=None):
    # check for solid and opacity range
    if (solid<0 or solid > 1):
        solid = 0.5
    if (opacity<0 or opacity > 1):
        opacity = 1   
    ##loop over primaries
    L = len(coordinates)
    for i in range(L):
        # check if primary (radius > 0)
        if (coordinates[i][3] > 0 ):
            # Draw primary sphere:
            draw_primary(ren, coordinates[i][0], coordinates[i][1], coordinates[i][2], coordinates[i][3], opacity=opacity, color=color,ambient = ambient, diffuse = diffuse, specular = specular)

    return None

def draw_primary(ren, x, y, z, r, ambient = 0.6, diffuse = 0.3, specular = 0.0,opacity = 0.2, color='Seashell'):
    """
    Draw a primary given a rendering window, coordiates, radius, and optical properties

    :param ren: The VTK render window
    :param x: x-coordinage (nm).
    :param y: y-coordinage (nm).
    :param z: z-coordinage (nm).
    :param r: sphere radius (nm).
    :param ambient: Sphere Ambient light, range [0,1], default = 0.6.
    :param diffuse: Sphere Light Diffusivity, range [0,1], default = 0.3.
    :param specular: Sphere Specular... yeah not sure..
    :param opacity: Sphere opacity, ranged [0,1], default = 0.2.
    :param color: String of named color to set sphere color, default 'Seashell'
    :return:
    """
    # Check color name: use seashell by default if not a real color:
    colors = vtk.vtkNamedColors()
    if(colors.ColorExists(color)):
        sphere_color = colors.GetColor3d(color)
    else:
        sphere_color = colors.GetColor3d("Seashell")
    # Check for other argument values, 
    # implement later
    # Sphere source:
    sphere = vtk.vtkSphereSource()
    sphere.SetThetaResolution(100)
    sphere.SetPhiResolution(50)
    sphere.SetRadius(r*1e9) # m
    # The mapper is responsible for pushing the geometry into the graphics
    # library. It may also do color mapping, if scalars or other attributes
    # are defined.
    #
    sphereMapper = vtk.vtkPolyDataMapper()
    sphereMapper.SetInputConnection(sphere.GetOutputPort())
    spheres_actor= vtk.vtkActor()
    spheres_actor.SetMapper(sphereMapper)
    spheres_actor.GetProperty().SetColor(sphere_color)
    spheres_actor.GetProperty().SetAmbient(ambient)
    spheres_actor.GetProperty().SetDiffuse(diffuse)
    spheres_actor.GetProperty().SetSpecular(specular)
    spheres_actor.GetProperty().SetOpacity(opacity)
    spheres_actor.AddPosition([x*1e9,y*1e9,z*1e9])
    ren.AddActor(spheres_actor)
    return None

def WriteImage(fileName, renWin, rgba=True):
    """
    Write the render window view to an image file.

    Image types supported are:
     BMP, JPEG, PNM, PNG, PostScript, TIFF.
    The default parameters are used for all writers, change as needed.

    :param fileName: The file name, if no extension then PNG is assumed.
    :param renWin: The render window.
    :param rgba: Used to set the buffer type.
    :return:
    """

    import os

    if fileName:
        # Select the writer to use.
        path, ext = os.path.splitext(fileName)
        ext = ext.lower()
        if not ext:
            ext = '.png'
            fileName = fileName + ext
        if ext == '.bmp':
            writer = vtk.vtkBMPWriter()
        elif ext == '.jpg':
            writer = vtk.vtkJPEGWriter()
        elif ext == '.pnm':
            writer = vtk.vtkPNMWriter()
        elif ext == '.ps':
            if rgba:
                rgba = False
            writer = vtk.vtkPostScriptWriter()
        elif ext == '.tiff':
            writer = vtk.vtkTIFFWriter()
        else:
            writer = vtk.vtkPNGWriter()

        windowto_image_filter = vtk.vtkWindowToImageFilter()
        windowto_image_filter.SetInput(renWin)
        windowto_image_filter.SetScale(1)  # image quality
        if rgba:
            windowto_image_filter.SetInputBufferTypeToRGBA()
        else:
            windowto_image_filter.SetInputBufferTypeToRGB()
            # Read from the front buffer.
            windowto_image_filter.ReadFrontBufferOff()
            windowto_image_filter.Update()

        writer.SetFileName(fileName)
        writer.SetInputConnection(windowto_image_filter.GetOutputPort())
        writer.Write()
    else:
        raise RuntimeError('Need a filename.')


def TEMframe():

    ##############################################
    ## Input parameters
    ##############################################
    
    ## location of data
    coords_file = './streamline_10_M0_2e13_splits_1_steps_1-video(43).csv'
    
    ## draw scale bar (width metres)?
    scale_bar_width = 10e-9
    
    ## frame number to draw ( < 0 for final frame) 
    n_frame = -1 
    # n_frame = 9
    
    ## half frame sizes in metres
    ## box will be [-x, x] 
    half_size_x =  50e-9
    half_size_y =  50e-9
    half_size_z =  50e-9
    

    ##############################################
    ## VTK set up
    ##############################################

    colors = vtk.vtkNamedColors()
    # Set the background color.
    colors.SetColor("bkg", [255, 255, 255, 255])

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)


    ren.SetBackground(colors.GetColor3d("bkg"))
    renWin.SetSize(500, 500)
    renWin.SetWindowName("Particle")

    light = vtk.vtkLight()
    light.SetFocalPoint(1.875, 0.6125, 0)
    light.SetPosition(0.875, 1.6125, 1)
    ren.AddLight(light)

    camera = vtk.vtkCamera()
    ren.SetActiveCamera(camera)
    ren.GetActiveCamera().SetFocalPoint(0, 0, 0)

    ##############################################
    ##############################################
        
    ## read time, primary coordinates, and frame coordinates
    frames = read_coords(coords_file)
    n = len(frames)
    
    ## frame to draw
    if (n_frame >= n-1 or n_frame < 0): 
        n_frame = n-1
    
    print('Drawing frame: '+str(n_frame)+' of '+str(n-1))
    
    #reorient frame
    new_frames = reorientate(frames, n_frame)
    
    ##draw particle

    draw_particle(ren, frames[n_frame]["coords"],1, color=black, opacity=0.15)
    
    # change axis color
    #ax.Color = [0.7 0.7 0.7]
    
    # ax.set_xlim([-half_size_x, half_size_x])
    # ax.set_ylim([-half_size_y, half_size_y])
    # ax.set_zlim([-half_size_z, half_size_z])
    
    ##light object
    # light('Position',[-1,-1,0])

    #draw scale bar in lower right corner
    # if (scale_bar_width > 0):
    #     sh = 0.8
        # fig.plot([-0.9*half_size_x, -0.9*half_size_x+scale_bar_width],[-sh*half_size_y,-sh*half_size_y],'-k','linewidth',6)
        # text(-0.9*half_size_x,-(sh+0.05)*half_size_y,strcat(num2str(scale_bar_width*1e9),' nm'),'fontsize',20)

    


    # Adjust Position and view before resetting camera to change view:
    ren.GetActiveCamera().SetPosition(10, 15, 100)
    ren.GetActiveCamera().SetViewUp(0, 1, 0)
    ren.GetActiveCamera().ParallelProjectionOn()
    ren.ResetCamera()

    # iren.Initialize()
    renWin.Render()
    WriteImage("test.png", renWin)
    # iren.Start()


def TEMvideo():
    ##############################################
    ## Input parameters
    ##############################################
    
    ## location of data
    coords_file = './streamline_10_M0_2e13_splits_1_steps_1-video(43).csv'
    # coords_file = './test.csv'
    
    ## draw scale bar (width metres)?
    scale_bar_width = 10e-9
    
    ## frame number to draw ( < 0 for final frame) 
    n_frame = -1 
    # n_frame = 9
    
    ## half frame sizes in metres
    ## box will be [-x, x] 
    half_size_x =  50e-9
    half_size_y =  50e-9
    half_size_z =  50e-9
    

    ##############################################
    ## VTK set up
    ##############################################

    colors = vtk.vtkNamedColors()
    # Set the background color.
    colors.SetColor("bkg", [255, 255, 255, 255])

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    ren.SetBackground(colors.GetColor3d("bkg"))
    renWin.SetSize(640, 480)
    renWin.SetWindowName("Ambient Spheres") 

    light = vtk.vtkLight()
    light.SetFocalPoint(1.875, 0.6125, 0)
    light.SetPosition(0.875, 1.6125, 1)
    ren.AddLight(light)

    camera = vtk.vtkCamera()
    ren.SetActiveCamera(camera)
    ren.GetActiveCamera().SetFocalPoint(0, 0, 0)

    ##############################################
    ##############################################
        
    ## read time, primary coordinates, and frame coordinates
    frames = read_coords(coords_file)
    n = len(frames)
    
    ## frame to draw
    if (n_frame >= n-1 or n_frame < 0): 
        n_frame = n-1

    if(debug):
        for j in range(n):
            print(frames[j].keys())
            print("Description of frame: {}".format(j))
            print("Time: {}".format(frames[j]['Time']))
            print("Coordinates")
            print(frames[j]['coords'])
            print("orientation_x")
            print(frames[j]['orient_x'])
            print("orientation_z")
            print(frames[j]['orient_z'])
    
    print('Drawing frame: '+str(n_frame)+' of '+str(n-1))
    
    #reorient frame
    new_frames = reorientate(frames)

    if(debug):
        print("done reorientate")
        for j in range(n):
            print(new_frames[j].keys())
            print("Description of frame: {}".format(j))
            print("Time: {}".format(new_frames[j]['Time']))
            print("Coordinates")
            print(new_frames[j]['coords'])

    ##draw particle
    # print(new_frames[n_frame]["coords"])
    # print(frames[n_frame]["coords"])

    draw_particle(ren, new_frames[n_frame]["coords"],1)
    if(debug):print("Drew particle")
    # change axis color
    #ax.Color = [0.7 0.7 0.7]
    
    # ax.set_xlim([-half_size_x, half_size_x])
    # ax.set_ylim([-half_size_y, half_size_y])
    # ax.set_zlim([-half_size_z, half_size_z])
    
    ##light object
    # light('Position',[-1,-1,0])

    #draw scale bar in lower right corner
    # if (scale_bar_width > 0):
    #     sh = 0.8
        # fig.plot([-0.9*half_size_x, -0.9*half_size_x+scale_bar_width],[-sh*half_size_y,-sh*half_size_y],'-k','linewidth',6)
        # text(-0.9*half_size_x,-(sh+0.05)*half_size_y,strcat(num2str(scale_bar_width*1e9),' nm'),'fontsize',20)

    # Adjust Position and view before resetting camera to change view:
    ren.GetActiveCamera().SetPosition(10, 15, 100)
    ren.GetActiveCamera().SetViewUp(0, 1, 0)
    ren.GetActiveCamera().ParallelProjectionOn()
    ren.ResetCamera()

    iren.Initialize()
    renWin.Render()
    WriteImage("last.png", renWin)
    # iren.Start()
    if(debug):
        print('Wrote last.png')
        for frame in frames:
            print(frame)
    for i in range(len(frames)):
        ren.RemoveAllViewProps()
        if("coords" in frames[i].keys()):
            print('Drawing frame: '+str(i)+' of '+str(n-1))
            draw_particle(ren, new_frames[i]["coords"],1)
            renWin.Render()
        WriteImage("test"+str(i)+".png", renWin)

    


if __name__ == '__main__':
    print("Running")
    sys.exit(TEMvideo())
