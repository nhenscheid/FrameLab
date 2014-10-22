function compileMex(x,arch,gpu)
% Compile mex files from source using various options


end



function mexFanA(arch,gpu)

    switch gpu
        case 1
            disp('Compiling forward fan beam transform for CPU/GPU')
            cd src
            mex -v -L'CUDALIB' -lcudart -I"./" Ax_fan_mf.cpp...
                Ax_fan_mf_cpu_siddon.cpp Ax_fan_mf_cpu_new.cpp...
                Ax_fan_mf_cpu_new_fb.cpp Ax_fan_mf_gpu_siddon.cu...
                Ax_fan_mf_gpu_new.cu Ax_fan_mf_gpu_new_fb.cu...
                find_area.cpp sort_alpha.cpp
        case 0
            disp('Compiling forward fan beam transform for CPU only')
    end

end


function mexFanAt(arch,gpu)


end


function mexConeA(arch,gpu)



end



function mexConeAt(arch,gpu)



end