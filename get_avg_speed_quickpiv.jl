using PythonCall
<<<<<<< HEAD
include("E:/PV/git_repo/multi_quickPIV/src/multi_quickPIV.jl")


function fn(a1, a2)
    IA,SM,ST,TH=50,20,15,100
    Us = []
    Vs = []
    pivparams = multi_quickPIV.setPIVParameters(interSize=IA, searchMargin=SM, step=ST, corr_alg="nsqecc", threshold=TH)
    VF, SN = multi_quickPIV.PIV(a1, a2, pivparams, precision=64)
    push!( Us, VF[1,:,:] );
    push!( Vs, VF[2,:,:] )
    avg_speed = [ sum( sqrt.( Us[t].^2 + Vs[t].^2 ) )/length(Us[t]) for t in 1:length(Us) ]
    return avg_speed
end
