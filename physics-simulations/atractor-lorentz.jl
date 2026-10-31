using Plots
using OrdinaryDiffEq
using DifferentialEquations
using GLMakie
using DataStructures: CircularBuffer
σ=Observable(Float64(10))
ρ=Observable(Float64(28))
β=Observable(8/3)
p=[σ[],ρ[],β[]]
fig = GLMakie.Figure(); display(fig)
ax = GLMakie.Axis3(fig[1,1],limits = ((-30, 30), (-30, 30), (0, 60)))
xyz=Observable([1,1,1])
function lorenz(du,u,p,t)
    du[1]=p[1]*(u[2]-u[1])
    du[2]=u[1]*(p[2]-u[3])-u[2]
    du[3]=u[1]*u[2]-p[3]*u[3]
end
xyz[]=[1,1,1]
u0=xyz[]
tspan=(0,75)
prob=OrdinaryDiffEq.ODEProblem(lorenz,u0,tspan,p,reltol=7e-8,abstol=7e-8)
sol=solve(prob,Tsit5())
n=length(sol.t)
Iniciar = Button(fig[2,1]; label = "Iniciar", tellwidth = false)
isrunning = Observable(false)
on(Iniciar.clicks) do clicks; isrunning[] = !isrunning[]; end
on(Iniciar.clicks) do clicks
    @async while isrunning[]
    isopen(fig.scene) || break
    for i in 1:n-1
        p1=Point3(sol.u[i])
        p2=Point3(sol.u[i+1])
        GLMakie.scatter!(ax, p1; marker = :circle, strokewidth = 1, 
        strokecolor = :blue,
        color = :black, markersize = (3,3,3)
        )
        GLMakie.lines!(ax, [p2,p1]; linewidth = 1, color = :red)
        sleep(0.001)
    end 
    end
end

