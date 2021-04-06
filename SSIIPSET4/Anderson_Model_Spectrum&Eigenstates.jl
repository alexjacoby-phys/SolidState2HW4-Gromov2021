using LinearAlgebra
using Combinatorics
using Plots

function Proj(L::Int64,X::Int64)
    Proj0 = zeros(Float64,L,L)
    Proj0[X,X] = 1.
    return Proj0
end

function Hop(L::Int64,X::Int64)
    Hop0 = zeros(Float64,L,L)
    Hop0[X,X+1] = 1.
    Hop0[X+1,X] = 1.
    return Hop0
end

function MakeH(L::Int64,w::Number,t::Number)
    Onsite = w*(rand(L)-0.5*ones(L))
    H=zeros(Float64,L,L)
    for i in 1:L-1
        H = H+ t*Hop(L,i)+Onsite[i]*Proj(L,i)
    end
    H=H + Onsite[L]*Proj(L,L)
    return H
end

function SolveH(L,w,t)
    H = MakeH(L,w,t)
    return (H,eigvals(H),eigvecs(H))
end




function MakeDOS(Spectrum::Vector{Float64},FractionalClarity::Int64)
    SpectrumPos = Spectrum - Spectrum[1]*ones(length(Spectrum))
    Increment = SpectrumPos[length(Spectrum)]/FractionalClarity
    EnergyAxis = Increment*Vector(1:FractionalClarity)
    DOS = zeros(FractionalClarity)
    for i in 1:FractionalClarity
        for j in 1:length(Spectrum)
            DOS[i]=length(filter(t -> abs(t-i*Increment) < Increment,Spectrum))
        end
    end
    DOS = DOS*(1/sum(DOS))*(1/Increment)
    return (EnergyAxis, DOS)
end

#This Function is Not Very Useful Unless You Are Only Solving Once#
function SOLVEMODEL(L::Int64,w::Number,t::Number,FractionalClarity::Int64)
    FirstReturn = SolveH(L,w,t)
    SecondReturn = MakeDOS(FirstReturn[2],FractionalClarity)
    return (FirstReturn,SecondReturn)
end


function SpecAVG(L::Int64,w::Number,t::Number,AvgNum::Int64)
    SpecAvg = SolveH(L,w,t)[2]
    for i in 1:AvgNum-1
        SpecAvg = SpecAvg + SolveH(L,w,t)[2]
    end
    SpecAvg = (1/AvgNum)*SpecAvg
    return SpecAvg
end



function InverseParticipation(L::Int64,EigSol::Array{Number,2},q::Number)
    States = Vector{Vector{Float64}}(undef,L)
    for i in 1:L
        States = normalize!(EigSol[i,:])
    end
    AVGINVPR = 0
    for i in 1:L#Indexes states
        INVPR= 0
        for j in 1:L #Indexes positions
            INVPR = INVPR + (norm(Proj(L,j)*States[j]))^q
        end
        AVGINVPR = AVGINVPR + INVPR*(1/L)
    end
    return AVGINVPR
end



L=500
w=2
t=-0.1
FractionalClarity = 50
Trials = 1
DOSGRAPH = MakeDOS(SpecAVG(L#=system size=#,w#=w parm=#,t#=t parm=#,Trials#=number of trials averaged over=#),FractionalClarity#=how many increments there are on Xaxis of DOSGRAPH=#)


plot(DOSGRAPH[1],L*DOSGRAPH[2],title = string("DOS w/ L= ",L,", W=", w, ", and t=",t,"; ",Trials," Trials Averaged"), label =  nothing, #=seriestype = :scatter,=# xlabel = "Energy", ylabel = "Approximate DOS")


savefig(string("DOSL=",L,"W=", w, "t=",t,";",Trials,"Trials Averaged"))
savefig("dos")
#since we assume occupation of 1 per site, we have filled the first L eigenstates#
