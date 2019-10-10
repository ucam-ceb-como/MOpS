%{  
 Author:    Casper Lindberg (csl37)

 Purpose:   Reorientation of particle coordinates to maintain frame
            position and orientation

 In:        array of frame structures, range of frames to reorientate (optional)     
 Out:       frame structure with reorientated coordinates
%}
function new_frames = reorientate(frames, frame_range)
    
    %check for frame_range (optional)
    %in not supplied reorientate all coordinates
    if (~exist('frame_range', 'var'))
        n_frames = length(frames);
        frame_range = [1:n_frames];
    end

    for j = frame_range     %loops over frames
        
        %loops over primaries
        [n_prim,~] = size(frames(j).coords);
        for i = 1:n_prim 
            
            %if this is a primary, (radius non-zero)
            if(frames(j).coords(i,4) ~= 0 ) 
                
                %%translate coordinates
                %%this recentres the frame on the tracked primary
                new_frames(j).coords(i,1:3) = translate_coords(frames(j).coords(i,1:3),frames(j).position(1:3));
                
                %%rotate coordinates after translation
                %%this returns the frame to the correct orientation 
                %%defined by the points on the x- and z-axes.
                new_frames(j).coords(i,1:3) = rotate_coords(new_frames(j).coords(i,1:3),frames(j).orient_x(1:3),frames(j).orient_z(1:3));
                
                %%copy remaining data
                new_frames(j).coords(i,4) = frames(j).coords(i,4);
                new_frames(j).t = frames(j).t;
                %%We haven't copied the frame orientation and position
                %%information. This is no longer needed.
            end
        end
    end
end

%% Translate coordinates to recentre frame
%%IN:   coordinates, new centre 
%%OUT:  new coordinates
function new_coords = translate_coords(coords,centre)
    
    new_coords(1,1) = coords(1) - centre(1);
    new_coords(1,2) = coords(2) - centre(2);
    new_coords(1,3) = coords(3) - centre(3);
end

%% Rotate coordinates to maintain constant frame orientation
%{
Once the frame has been recentred on the tracked primary the particle is
rotated to place the points orient_x and orient_z on the x- and z-axes
respectively.
%}
%%IN:   coordinates to rotate, x-axis orientation, z-axis orientation
%%OUT:  new coordinates
function new_coords = rotate_coords(coords,orient_x,orient_z)
    
    %calculate angles for rotation of frame z-axis
    len_z = sqrt(orient_z(1)^2+orient_z(2)^2+orient_z(3)^2);
    if(orient_z(1)~=0)    %no x-component
        phi1 = atan(orient_z(2)/orient_z(1)); 
        if(orient_z(1) > 0 && orient_z(2) < 0)
            phi1 = phi1 + 2*pi;
        end
        if(orient_z(1) < 0)
            phi1 = phi1 + pi;
        end
    else 
        phi1 = pi/2;
    end
    theta1 = acos(orient_z(3)/len_z);
    
    %construct rotatation matrix
    rot1(1,1) = cos(theta1)*cos(phi1);
    rot1(1,2) = cos(theta1)*sin(phi1);
    rot1(1,3) = -sin(theta1);
    rot1(2,1) = -sin(phi1);
    rot1(2,2) = cos(phi1);
    rot1(2,3) = 0;
    rot1(3,1) = sin(theta1)*cos(phi1);
    rot1(3,2) = sin(theta1)*sin(phi1);
    rot1(3,3) = cos(theta1);
    
    %rotate coords
    coords1 = mtimes(rot1,transpose(coords));
    
    %rotate x orientation vector
    xorient = mtimes(rot1,transpose(orient_x(1:3)));
    
    %%calculate rotation about z-axis
    len_x = sqrt(xorient(1)^2+xorient(2)^2+xorient(3)^2);
    theta2 = acos(xorient(1)/len_x);
    if(xorient(2)<0) %get correct range for  rotation
        theta2 = 2*pi - theta2;
    end
    
    %construct rotation matrix;
    rot2(1,1) = cos(theta2);
    rot2(1,2) = sin(theta2);
    rot2(1,3) = 0;
    rot2(2,1) = -sin(theta2);
    rot2(2,2) = cos(theta2);
    rot2(2,3) = 0;
    rot2(3,1) = 0;
    rot2(3,2) = 0;
    rot2(3,3) = 1;
    
    %%rotate coords
    coords2 = mtimes(rot2,coords1);
    new_coords = transpose(coords2);
    
end