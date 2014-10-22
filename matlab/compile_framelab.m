% FrameLab install script 

arch=computer;
mac=strcmp(arch,'MACI64') || strcmp(arch,'MACI') || strcmp(arch,'MAC');
linux=strcmp(arch,'GLNXA64') || strcmp(arch,'GLNX86');
GPU = gpuDeviceCount;

CUDALIB = '/usr/local/cuda/lib64';