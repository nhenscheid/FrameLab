
close all;clear;clc;
if 1
    !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Ax_fan_mf_gpu_siddon.cu
    !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Ax_fan_mf_gpu_new.cu
    !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Ax_fan_mf_gpu_new_fb.cu
    path='/usr/local/cuda-6.5/lib';
    mex(['-L' path],'-lcudart','Ax_fan_mf.cpp','Ax_fan_mf_cpu_siddon.cpp','Ax_fan_mf_cpu_new.cpp','Ax_fan_mf_cpu_new_fb.cpp',...
        'Ax_fan_mf_gpu_siddon.o','Ax_fan_mf_gpu_new.o','Ax_fan_mf_gpu_new_fb.o','sort_alpha.cpp','find_area.cpp');

    !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Atx_fan_mf_gpu_new.cu
    !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Atx_fan_mf_gpu_new_fb.cu
    path='C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.0\lib\x64';
    mex(['-L' path],'-lcudart','Atx_fan_mf.cpp','Atx_fan_mf_cpu_new.cpp','Atx_fan_mf_cpu_new_fb.cpp',...
        'Atx_fan_mf_gpu_new.obj','Atx_fan_mf_gpu_new_fb.obj','find_l.cpp','find_area.cpp');
end

if 1
   !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Ax_cone_mf_gpu_siddon.cu
   !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Ax_cone_mf_gpu_new.cu
    path='C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.0\lib\x64';
    mex(['-L' path],'-lcudart','Ax_cone_mf.cpp','Ax_cone_mf_cpu_siddon.cpp','Ax_cone_mf_cpu_new.cpp',...
        'Ax_cone_mf_gpu_siddon.obj','Ax_cone_mf_gpu_new.obj','sort_alpha.cpp');
    
    !"%VS90COMNTOOLS%vsvars32.bat" & nvcc -c -m64 -arch compute_20 -code sm_20 Atx_cone_mf_gpu_new.cu
    path='C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.0\lib\x64';
    mex(['-L' path],'-lcudart','Atx_cone_mf.cpp','Atx_cone_mf_cpu_new.cpp',...
        'Atx_cone_mf_gpu_new.obj','find_l_3d.cpp');
end