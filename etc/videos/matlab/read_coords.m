%{  
 Author:    Casper Lindberg (csl37)

 Purpose:   Functions to read in particle and frame coordinates for videos
            from MOpS "-videos().csv" file format

 In:        File path
 Out:       Array of frame structures
            frame.t:        time
            frame.coords:   list of primary coordinates (x,y,z,r)
            frame.position: frame position relative to particle
            frame.orient_x: frame x orientation vector
            frame.orient_z: frame z orientation vector
%}
function [frame] = read_coords(file_name)

    %%read data 
    data = csvread(file_name,1,0);
    %data = readmatrix(file_name);
    [m,~] = size(data);
    
    %%initialize counters
    i=1;    %primary particle counter
    j=1;    %frame counter
    
    %%first time coordinate
    frame(j).t = data(1,1);
    
    for k = 1:m
        
        %increment frame counter
        if data(k,1) ~= frame(j).t
            i = 1;
            j = j+1;
            frame(j).t = data(k,1);
        end
        
        %add primary particle coordinates (x,y,z,r)
        frame(j).coords(i,1:4) = data(k,2:5);
        i=i+1;
        
        % get frame position and orientation 
        % the frame position are the position coordinates of the tracked primary
        % the frame orientation is tracked using two orthogonal unit
        % vectors: the frame x and z axes
        if (data(k,6)*data(k,6)+ data(k,7)*data(k,7)+ data(k,8)*data(k,8)) ~= 0   % tracked primary has non-zero z vector    
           %tracked primary position
           frame(j).position(1:3) = data(k,2:4);
           %tracked frame orientation vectors
           frame(j).orient_x(1:3) = data(k,6:8);    %x vector
           frame(j).orient_z(1:3) = data(k,9:11);   %z vector
        end
    end
    
end


