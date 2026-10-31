ENV["PYTHON"]=""
using Pkg
Pkg.build("PyCall")
Pkg.add("SymPy")
using SymPy
using Pkg, GLMakie, Plots, OrdinaryDiffEq, ModelingToolkit, SpecialFunctions, QuadGK, FastGaussQuadrature
#Resolver la ecuacion diferencial con valores iniciales; airyaiprime(-15), airyai(-15), airybiprime(-15), airybi(-15)
f(du,u,p,t)=u*t
duA=airyaiprime(-15)
uA=airyai(-15)
duB=airybiprime(-15)
uB=airybi(-15)
tspan=(-15,5)
prob1=SecondOrderODEProblem(f,duA,uA,tspan)
prob2=SecondOrderODEProblem(f,duB,uB,tspan)
solAi=solve(prob1, Vern9(), reltol=1e-16, abstol=1e-16)
solBi=solve(prob2, Vern9(), reltol = 1e-16, abstol = 1e-16)
Ai=[tupla[2] for tupla in solAi.u]
Bi=[tupla[2] for tupla in solBi.u]
x_vals = solBi.t
airy_valsAi = [airyai(x)[1] for x in x_vals]
airy_valsBi = [airybi(x)[1] for x in x_vals]

Plots.plot()
Plots.plot(solAi.t,Ai,ylims=(-1,1),label="Ai(x) Numérica")
Plots.plot!(solBi.t,Bi,ylims=(-1,1),label="Bi(x) Numérica")
Plots.plot!(x_vals, airy_valsAi, label="Función Ai(x) ", linestyle=:dash, linewidth=2)
Plots.plot!(x_vals, airy_valsBi, label="Función Bi(x) ", linestyle=:dash, linewidth=2)
Plots.xlabel!("x")
Plots.ylabel!("y")
Plots.title!("Comparación de la Solución Numérica")

#Pkg.add("ApproxFun")
#using ApproxFun
#using SpecialFunctions
#a, b = -15, 10
#d = a..b;
#x = Fun(d);
#D = Derivative(d);
#L = D^2 - x;
#B = Dirichlet(d);
#B_vals = [airyai(a), airyai(b)];
#B_vals2 = [airybi(a), airybi(b)];
#u = [B; L] \ [B_vals, 0];
#u2 = [B;L] \ [B_vals2, 0];
#Plots.plot(u2 , xlabel="x", ylabel="u(x)", xlims=[-15,10],ylims=[-0.5,1] ,legend=false)
#No mirar esto de aquí arriba en verde.

#Resolver la ecuación mediante integración numérica
#Sin cambio de variable
function funAi(x,t)
    return (1/SymPy.pi)*(cos(t^3/3+x*t))
end
function funBi(x,t)
    return (1/SymPy.pi)*((exp(-(1/3)*t^3+x*t)+sin(t^3/3+x*t)))
end
#Con cambio de variable
function funAi(x,t)
    return (1/SymPy.pi)*(cos((t/(1-t))^3/3+x*(t/(1-t))))*(1/(1-t)^2)
end
function funBi(x,t)
    return (1/SymPy.pi)*((exp(-(1/3)*(t/(1-t))^3+x*(t/(1-t)))+sin((t/(1-t))^3/3+x*(t/(1-t)))))*(1/(1-t)^2)
end
function integralAi(x)
    result, _ =quadgk(t -> funAi(x,t),0,0.9)
    return result
end
function integralBi(x)
    result, _ =quadgk(t -> funBi(x,t),0,0.9)
    return result
end
x_values=LinRange(-15,5,1000)
f_valuesAi=[integralAi(x) for x in x_values]
f_valuesBi=[integralBi(x) for x in x_values]
Plots.plot()
Plots.plot!(x_values,f_valuesAi, ylims=(-1,1),label="Ai(x) Integración Numérica")
Plots.plot!(x_values,f_valuesBi, ylims=(-1,1),label="Bi(x) Integración Numérica")
Plots.plot!(x_vals, airy_valsAi, label="Función Ai(x) ", linestyle=:dash, linewidth=2)
Plots.plot!(x_vals, airy_valsBi, label="Función Bi(x) ", linestyle=:dash, linewidth=2)
Plots.xlabel!("x")
Plots.ylabel!("y")
Plots.title!("Solución a la ecuacion de airy mediante integración numérica")

#Resolver mediante expansión en serie de potencias
x_values=LinRange(-15,5,1000)
n=150
Plots.plot()
frames=1:n
anim=@animate for frame in frames
    function AiSerie(x)
        return (1/(3^(2/3)*pi))*sum((1/(factorial(big(k))))*gamma((k+1)/3)*sin((2*pi*(k+1))/3)*(3^(1/3)*x)^k for k in 0:frame)
    end
    function BiSerie(x)
        return (1/(3^(1/6)*pi))*sum((1/(factorial(big(k))))*gamma((k+1)/3)*abs(sin((2*pi*(k+1))/3))*(3^(1/3)*x)^k for k in 0:frame)
    end
    valoresSerieAi=[AiSerie(x) for x in x_values]
    valoresSerieBi=[BiSerie(x) for x in x_values]
    Plots.plot(x_values,valoresSerieAi,ylims=(-1,1),label="Ai(x) Serie",legend=:bottomright)
    Plots.plot!(x_values,valoresSerieBi,ylims=(-1,1),label="Bi(x) Serie",legend=:bottomright)
    Plots.plot!(x_vals, airy_valsAi, label="Función Ai(x) ", linestyle=:dash, linewidth=2)
    Plots.plot!(x_vals, airy_valsBi, label="Función Bi(x) ", linestyle=:dash, linewidth=2)
    Plots.xlabel!("x")
    Plots.ylabel!("y")
    Plots.title!("Ecuación de airy por series de potencias")
end
gif(anim, "AiryPotencias.gif", fps = 30)
