CUDALIB = '/usr/local/cuda/lib64'

str = sprintf('-v -L"%s" -lcudart -I"./"',CUDALIB);
str = [str,' Ax_fan_mf.cpp',' Ax_fan_mf_cpu_siddon.cpp',...
    ' Ax_fan_mf_cpu_new.cpp',' Ax_fan_mf_cpu_new_fb.cpp',...
    ' Ax_fan_mf_gpu_siddon.cu',' Ax_fan_mf_gpu_new.cu',...
    ' Ax_fan_mf_gpu_new_fb.cu',' find_area.cpp', ' sort_alpha.cpp'];

args = regexp(str, '\s+', 'split');

mex(args{:})


