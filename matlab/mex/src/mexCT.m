%% Before starting, read the matlab article 
% http://www.mathworks.com/help/distcomp/run-mex-functions-containing-cuda-code.html
% In particular, you must copy over the mex_CUDA_*.xml to the build
% directory and change any relevant settings.  This assumes you know which
% CUDA compile options you need, e.g. compute_30, etc.  You will also need
% to set any relevant environment variables so that matlab can find your
% nvcc and g++ compilers.


%% Compile the forward fan beam transform: CPU/GPU combo
mex -v -L"/usr/local/cuda/lib64" -lcudart -I"./" Ax_fan_mf.cpp...
    Ax_fan_mf_cpu_siddon.cpp Ax_fan_mf_cpu_new.cpp...
    Ax_fan_mf_cpu_new_fb.cpp Ax_fan_mf_gpu_siddon.cu...
    Ax_fan_mf_gpu_new.cu Ax_fan_mf_gpu_new_fb.cu...
    find_area.cpp sort_alpha.cpp

%% Compile the adjoint fan beam transform: CPU/GPU combo
mex -v -I./ Atx_fan_mf.cpp...
    Atx_fan_mf_gpu_new.cu Atx_fan_mf_gpu_new_fb.cu ...
    Atx_fan_mf_cpu_new.cpp Atx_fan_mf_cpu_new_fb.cpp...
    find_l.cpp find_area.cpp -L/usr/local/cuda/lib64 -lcudart 

%% Compile the forward cone beam transform for the GPU
mex -L/usr/local/cuda/lib64 -lcudart -I./ Ax_cone_mf.cpp Ax_cone_mf_cpu_siddon.cpp...
    Ax_cone_mf_cpu_new.cpp Ax_cone_mf_gpu_siddon.cu Ax_cone_mf_gpu_new.cu...
    sort_alpha.cpp

%% Compile the adjoint cone beam transform for the GPU
mex -v LDFLAGS='-pthread -Wl,--no-undefined -Wl,--verbose'  -L/usr/local/cuda/lib64 -lcudart -I./ Atx_cone_mf.cpp Atx_cone_mf_cpu_new.cpp...
    find_l_3d.cpp Atx_cone_mf_gpu_new.cu 


%% Compile the forward fan beam transform CPU ONLY
mex -v -I"./" Ax_fan_mf_cpu.cpp Ax_fan_mf_cpu_new.cpp find_area.cpp sort_alpha.cpp

%% Compile the adjoint fan beam transform CPU ONLY
cd mex/source' codes'/;
mex -v -I"./" Atx_fan_mf_cpu.cpp Atx_fan_mf_cpu_new.cpp find_l.cpp find_area.cpp
movefile('Atx_fan_mf_cpu.mexa64','../')
cd ../../