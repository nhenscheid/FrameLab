% FrameLab install script 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY PATHS: CHANGE TO REFLECT YOUR SYSTEM
CUDALIB = '/usr/local/cuda/lib64'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get architecture
arch=computer;
mac=strcmp(arch,'MACI64') || strcmp(arch,'MACI') || strcmp(arch,'MAC');
linux=strcmp(arch,'GLNXA64') || strcmp(arch,'GLNX86');
GPU = sign(gpuDeviceCount) % To compile with gpu or not
if mac
    mexext = '.mexmaci64'; % Assuming 64 bit.
end
if linux
    mexext = '.mexa64'; % Assuming 64 bit.
end

% Which binaries to compile?  (add or remove if you need) 
switch GPU
    case 1
        %bins = {'Ax_fan_mf','Atx_fan_mf','Ax_cone_mf','Atx_cone_mf'};
        %bins = {'Ax_cone_mf','Atx_cone_mf'};
        %bins = {'Ax_fan_mf','Atx_fan_mf'};
        bins = {'Ax_cone_mf'};
    case 0
        %bins = {'Ax_fan_mf_cpu','Atx_fan_mf_cpu','Ax_cone_mf_cpu',...
        %        'Atx_cone_mf_cpu'};
        bins = {'Ax_cone_mf_cpu'};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% *** DO NOT MODIFY BELOW THIS LINE *** %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%!!NOTE!!: Must leave a space at the beginning of each file name!
Ax_fan_mf = [' Ax_fan_mf.cpp',' Ax_fan_mf_cpu_siddon.cpp',...
           ' Ax_fan_mf_cpu_new.cpp',' Ax_fan_mf_cpu_new_fb.cpp',...
           ' Ax_fan_mf_gpu_siddon.cu',' Ax_fan_mf_gpu_new.cu',...
           ' Ax_fan_mf_gpu_new_fb.cu',' find_area.cpp', ' sort_alpha.cpp'];
       
Atx_fan_mf = [' Atx_fan_mf.cpp',' Atx_fan_mf_gpu_new.cu',...
            ' Atx_fan_mf_gpu_new_fb.cu',' Atx_fan_mf_cpu_new.cpp',...
            ' Atx_fan_mf_cpu_new_fb.cpp',' find_l.cpp', ' find_area.cpp'];
        
Ax_cone_mf = [' Ax_cone_mf.cpp',' Ax_cone_mf_cpu_siddon.cpp',...
              ' Ax_cone_mf_cpu_new.cpp',' Ax_cone_mf_gpu_siddon.cu',...
              ' Ax_cone_mf_gpu_new.cu',' sort_alpha.cpp'];
          
Atx_cone_mf = [' Atx_cone_mf.cpp',' Atx_cone_mf_cpu_new.cpp',...
               ' find_l_3d.cpp',' Atx_cone_mf_gpu_new.cu'];
           
Ax_fan_mf_cpu = [' Ax_fan_mf_cpu.cpp',' Ax_fan_mf_cpu_new.cpp',...
                 ' find_area.cpp',' sort_alpha.cpp'];
             
Ax_cone_mf_cpu = [' Ax_cone_mf_cpu.cpp',...
                  ' Ax_cone_mf_cpu_new.cpp',' sort_alpha.cpp'];

              
GPUFLAGS = sprintf('-v -L"%s" -lcudart -I"./"',CUDALIB);
%GPUFLAGS = sprintf('-v -L -lcudart -I"./"');
CPUFLAGS = '-v -I"./"';

switch GPU
    case 1
        % Compile binaries for CPU/GPU
        disp('Compiling for CPU/GPU');
        for i=1:length(bins)
            str = eval(bins{i});
            str = [GPUFLAGS,str];
            args = regexp(str, '\s+', 'split');
            cd mex/src
            mex(args{:})
            movefile(strcat(bins{i},mexext),'../')
            cd ../../
        end
        
    case 0
        % Compile binaries for CPU only 
        for i=1:length(bins)
            disp('Compiling for CPU only');
            str = eval(bins{i});
            str = [CPUFLAGS,str];
            args = regexp(str, '\s+', 'split');
            cd mex/src
            mex(args{:})
            movefile(strcat(bins{i},mexext),'../')
            cd ../../
        end
        
end