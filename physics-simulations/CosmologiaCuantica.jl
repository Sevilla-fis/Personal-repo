using Pkg, Plots, DifferentialEquations, GLMakie, Distributions, Random, LinearAlgebra, ModelingToolkit

#En primer lugar se definen las constantes

const G = 6.67430e-11  # m^3/kg/s^2
const c = 3.0e8  # m/s
const H0 = 70.4  # km/(s*Mpc)
k = 0  # (0 = universo plano)
Λ= 1.76e-52
# Densidades actuales (en kg/m^3)
ρ_r0 = 4.2e-31  # Radiación
ρ_m0 = 3e-27  # Materia fría (materia oscura + bariónica)
ρ_L0 = 6e-27  # Energía oscura

# Ecuación diferencial de Friedmann

function friedmann!(du, u, p, t)
    a = u[1]
    ρ=ρ_r0 * a^(-4) + ρ_m0 * a^(-3) + ρ_L0
    H = sqrt((8π*G/3) * (ρ) - (k * c^2) / a^2+Λ*c^2/3)
    du[1] = a * H  # da/dt = a * H(a)
end

# Condiciones iniciales

a0 = [0.001]  # Factor de escala inicial
tspan = (0, 5.35e17)  # Segundos

# Resolver la ecuación

prob = ODEProblem(friedmann!, a0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-14)
a_vals = [u[1] for u in sol.u]
tiempo=sol.t/(3.154e16) #aprox 20 M.M de años

# Graficar resultado

Plots.plot(tiempo,a_vals,label="a(t)")
Plots.xlabel!("Tiempo")
Plots.ylabel!("a(t)")
Plots.title!("Tamaño relativo del universo")
n=length(a_vals)

#Funcion para generar un circulo (h,k centro y r radio)

function circle(h,k,r)
    θ=LinRange(0,2*pi,500)
    U=h.+r.*sin.(θ),k.+r.*cos.(θ)
    return U
end

#Funcion para generar cuerpos en funcion de r y θ aleatorios

function body(h,k,r,θ)
    B=[h.+r.*sin.(θ),k.+r.*cos.(θ)]
    return B
end

#Si se generan mas de dos cuerpos, a la hora de animar sus interacciones gravitatorias ocurre un problema, parece que las velocidades de estos cuerpos van ligadas
#Es decir, si dos cuerpos por ejemplo se juntan, se frenan todos los demás cuerpos hasta que estos se separen
#No es necesario tener nada de esto en cuenta para la simulación sin gravedad
#Se recomienda poner 10 cuerpos para ver la atracción gravitatoria y el problema con esta, y 1000 si solo se quiere ver la expansión del universo

ncuerpos=20

#Definimos posiciones aleatorias de los cuerpos

b=[]
r=[]
θ=[]
for j in 1:ncuerpos
    push!(r,rand())
    push!(θ,rand(Uniform(0,2*pi)))
    push!(b,(body(0,0,r[j],θ[j])[1],body(0,0,r[j],θ[j])[2]))
end
frames=1:n

#Creo la animación, dibujando el universo teniendo en cuanto su tamaño en cada iteración.
#Para los cuerpos hago lo mismo teniendo en cuenta su distancia al origen (r).

anim=@animate for i in frames
    limgrafica=2
    valores=sum(a_vals[1:i]) #Sumatorio del factor de escala, el total da el tamaño final
    valorescte=sum(repeat(a0,n)[1:i]) #Valores pero para un universo con expansión constante
    Plots.plot(valores.*circle(0,0,1),xlims=(-limgrafica,limgrafica),ylims=(-limgrafica,limgrafica),seriestype=[:shape,],lw=0.5,c=:black,linecolor=:orange,legend=false,fillalpha=0.2,aspectratio=1)
    Plots.scatter!([valores.*b[k] for k in 1:ncuerpos],xlims=(-limgrafica,limgrafica),ylims=(-limgrafica,limgrafica),ms=:2)
    Plots.plot!(valorescte.*circle(0,0,1),xlims=(-limgrafica,limgrafica),ylims=(-limgrafica,limgrafica),seriestype=[:shape,],lw=0.5,c=:black,linecolor=:blue,legend=false,fillalpha=0.2,aspectratio=1)
    Plots.scatter!([valorescte.*b[k] for k in 1:ncuerpos],xlims=(-limgrafica,limgrafica),ylims=(-limgrafica,limgrafica),ms=:2)
end
gif(anim,"ExpUniverse.gif",fps = 30)

#Definimos las constantes para calcular la tasa de expansión

Ω_m=0.315
Ω_r=9e-5
Ω_Λ=0.685

#Termino que va dentro de la raíz de la fórmula H=H0*√(vec)

ve=Ω_m*[a_vals[i]^-3 for i in 1:n]+Ω_r*[a_vals[i]^-4 for i in 1:n].+Ω_Λ 

#Calculamos la tasa de expansión y ploteamos

H=H0*[sqrt(ve[i]) for i in 1:length(ve)] 
Plots.plot(tiempo,H,label="H(t)",ylims=(0,500),xlims=(0,20))
Plots.xlabel!("Tiempo (M.M. Años)")
Plots.ylabel!("H(t) (km/(s*Mpc))")
Plots.title!("Velocidad de expansión del universo")

#No mirar a partir de aquí-----------------
#A partir de aqui no está terminado, estoy intentando implementar fuerza gravitatoria a los cuerpos

#Las masas son las de las partículas elementales que había al inicio del universo

m0=[1.673*10^-27,3.344*10^-27,5.008*10^-27,6.646*10^-27,1.162*10^-26,1.165*10^-26]
masas=[]
for j in 1:ncuerpos
    push!(masas,m0[rand(1:length(m0))])
end

# Condiciones iniciales (concatenar posiciones y velocidades)

function v0func()
    v0=[0.0,0.0]
    return v0
end
v0=[]
for j in 1:ncuerpos
    push!(v0,v0func())
end
b1=a_vals[1]*[]
for j in 1:ncuerpos
    push!(b1,body(0,0,r[j],θ[j]))
end
u0 = vcat(vcat(b1...), vcat(v0...))

#Escribimos una matriz para las posiciones iniciales(b0) y (v0)

matrixb0=reshape(u0[1:2*ncuerpos],2,ncuerpos)
matrixv0=reshape(u0[2*ncuerpos+1:end],2,ncuerpos)

const ∑ = sum

#Definimos los vectores x e y que vamos a querer calcular

@variables t u(t)[1:2, 1:ncuerpos]
u = collect(u)
D = Differential(t)

#Definimos el potencial para la fuerza gravitatoria en n cuerpos
potential = -G *
            ∑(
    i -> ∑(j -> (masas[i] * masas[j]) / √(∑(k -> (u[k, i] - u[k, j])^2, 1:2)), 1:(i - 1)),
    2:ncuerpos)
    
#Definimos las ecuaciones con el Hamiltoniano para los n-cuerpos, teniendo en cuenta que dp/dt=-∇V y dq/dt=(1/M)*p

    eqs = vec(@. D(D(u))) .~ .-ModelingToolkit.gradient(potential, vec(u)) ./
                         repeat(masas, inner = 2)

#Con nuestro sistema definido, resolvemos la EDO

@named sys = ODESystem(eqs, t)
ss = structural_simplify(sys)
prob = ODEProblem(ss, [vec(u .=> matrixb0); vec(D.(u) .=> matrixv0)], tspan)
sol = solve(prob, Tsit5(), abstol=1e-16, reltol=1e-17);


# Extraer las soluciones

t_values = sol.t

#Igualamos el numero de datos en las posicones de los cuerpos y los datos en a(t), para poder animarlos conjuntamente

datau=sol.u[1:convert(Int64,ceil(length(sol.t)/n)):end]
datat=sol.t[1:convert(Int64,ceil(length(sol.t)/n)):end]
frames=1:length(datau)
anim=@animate for i in frames
    valores=sum(a_vals[1:i]) #Sumatorio del factor de escala, el total da el tamaño final
    valorescte=sum(repeat(a0,n)[1:i]) #Valores pero para un universo con expansión constante
    Plots.plot(valores.*circle(0,0,1),xlims=(-1,1),ylims=(-1,1),seriestype=[:shape,],lw=0.5,c=:black,linecolor=:orange,legend=false,fillalpha=0.2,aspectratio=1)
    Plots.plot!(valorescte.*circle(0,0,1),xlims=(-1,1),ylims=(-1,1),seriestype=[:shape,],lw=0.5,c=:black,linecolor=:blue,legend=false,fillalpha=0.2,aspectratio=1)
#Representamos las posiciones de los cuerpos a lo largo del tiempo, sabiendo que el vector sol.u nos da en los 6 primeros términos las velocidades de cada cuerpo
#Y en los 6 ultimos las posiciones de cada cuerpo (Esto suponiendo que dibujemos 3 cuerpos, para mas o menos cuerpos seria equivalente)
    Plots.scatter!(datau[i][2*ncuerpos+1:2:end],datau[i][2*ncuerpos+2:2:end],xlims=(-1,1),ylims=(-1,1),ms=:2)
end
gif(anim,"ExpUniverse.gif",fps = 30)
#En la simulación los cuerpos se mueven lo que parece atraidos por una fuerza gravitatoria de manera correcta,
#Pero sus velocidades parece que van ligadas
#Todavia no está implementada la expansión del universo a la atracción gravitatoria, simplemente se ve como se atraen los cuerpos
-












#test para sumar cuerpo
r_fusion=1e-5
function detectar_colisiones(posiciones, r_fusion)
    n = length(posiciones)
    colisiones = []
    for i in 1:n
        for j in i+1:n
            r_ij = norm(posiciones[i] - posiciones[j])
            if r_ij < r_fusion
                push!(colisiones, (i, j))
            end
        end
    end
    return colisiones
end
function fusionar_particulas!(posiciones, velocidades, masas, i, j)
    # Calcular la masa total
    masa_total = masas[i] + masas[j]

    # Calcular la posición del centro de masa
    r_cm = (masas[i] * posiciones[i] + masas[j] * posiciones[j]) / masa_total

    # Calcular la velocidad del centro de masa
    v_cm = (masas[i] * velocidades[i] + masas[j] * velocidades[j]) / masa_total

    # Reemplazar la partícula i con la nueva partícula fusionada
    posiciones[i] = r_cm
    velocidades[i] = v_cm
    masas[i] = masa_total

    # Eliminar la partícula j
    deleteat!(posiciones, j)
    deleteat!(velocidades, j)
    deleteat!(masas, j)
end

function test(x)
    if x==1 || typeof(x)==Float64
        return "True"
    end
end
test(1.0)