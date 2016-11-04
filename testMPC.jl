module testMPC
println("loading external Modules.......")
using JuMP
using Ipopt
using PyPlot
println("loaded")

v_ref = 1.0  #1.0 -> these values make infeasible with hard constraint
dt = 0.1    #0.1
N = 10      #10
zInit = zeros(3)# just to initialize the parameter, not initial states
x_obst = 4.0 #4.0 # x coordinate of center of ellipse
y_obst = 0.4 #2.0 # y coordinate of center of ellipse
rx = 0.2      #0.2 # "radius" of the elipses along x axis ## total length along x: 2*a
ry = 0.5     #2.35 #"radius" of the elipses along y axis


println("initialize Model.......")

m = Model(solver = IpoptSolver(print_level=0, max_cpu_time=1.0))
@variable(m, x[1:N+1,1:3] )
@variable(m, u[1:N, 1:2]  )

#for soft constraint
@variable(m, z[1:N+1]  )


@NLparameter(m, x0[j=1:3] == zInit[j] )


@NLconstraint(m, [i=1:3], x[1,i] == x0[i])


@NLconstraint(m, [i=1:N], x[i+1,1] == x[i,1]+ x[i,3]*cos(u[i,2])*dt)
@NLconstraint(m, [i=1:N], x[i+1,2] == x[i,2]+ x[i,3]*sin(u[i,2])*dt)
@constraint(m, [i=1:N], x[i+1,3] == x[i,3]+ u[i,1]*dt)


#hard constraint
#@NLconstraint(m, [i=1:N+1], ((x[i,1]-x_obst)/rx)^2+((x[i,2]-y_obst)/ry)^2 >= 1)

#soft constraint
#@NLconstraint(m, [i=1:N+1], ((x[i,1]-x_obst)/rx)^2+((x[i,2]-y_obst)/ry)^2 == 1+z[i]^2)
#@NLexpression(m, costObstacle, sum{-0.02*log(0.1+ z[i]^2),i=1:N+1})

#@NLexpression(m, costObstacle, sum{0.02*1/(  ((x[i,1]-x_obst)/rx)^2+((x[i,2]-y_obst)/ry)^2 - 1 )^2,i=1:N+1})
@NLexpression(m, costObstacle, sum{-0.02*log((  ((x[i,1]-x_obst)/rx)^2+((x[i,2]-y_obst)/ry)^2 - 1 )^2),i=1:N+1})

#!!
@NLexpression(m, costStatesInputs, sum{x[i,2]^2+ u[i,1]^2+u[i,2]^2+(x[i,3] -v_ref)^2, i = 1:N}+ x[N+1,2]^2+(x[N+1,3] -v_ref)^2)
@NLobjective(m, Min, costStatesInputs + costObstacle)
#print(m)


buffer =300
xCurr = zeros(buffer,3)
xCurr[1,1]= 0.0
xCurr[1,2]= 0.0
xCurr[1,3]= 0.2
uCurr = zeros(buffer,2)
x_log = zeros(N+1,3,buffer)
u_log = zeros(N,2,buffer)
tCurr = zeros(buffer)
costs = zeros(buffer, 2)

i  = 1
finished    = false

println("Solving..........")

while i < buffer && !finished

    setvalue(x0,xCurr[i,:])
    tCurr[i] = (i-1)*dt
    status = solve(m)
    #Simulate
    uCurr[i,:] = getvalue(u[1,:])
    x_log[:,:,i] = getvalue(x)
    u_log[:,:,i] = getvalue(u)
    xCurr[i+1,1]= xCurr[i,1]+xCurr[i,3]*cos(uCurr[i,2])*dt
    xCurr[i+1,2]= xCurr[i,2]+xCurr[i,3]*sin(uCurr[i,2])*dt
    xCurr[i+1,3]= xCurr[i,3]+uCurr[i,1]*dt


    costs[i,:] = [getvalue(costStatesInputs),getvalue( costObstacle)]

    if xCurr[i+1]>=10
        finished= true
        break
    end
   
    
    i=i+1

end

println("Finished Solving")


###post processing
i_final = i
close("all")

figure(2)#x-y plot of "car" and obstacle
subplot(211)
x1 = collect(-rx+x_obst:0.001:rx+x_obst) # calculate the domain of valid x for the ellipse
x1= x1[2:end-1] # ommit the first and last values because rounding errors sometimes give something like -1e-15 except for 0, so we cant calculate sqaure route
y = y_obst + ry*sqrt( 1-( ( ( x1-x_obst )/rx ).*( ( x1-x_obst )/rx ) ) ) #calculate upper of of ellipse
plot(x1,y)
y = y_obst - ry*sqrt( 1-( ( ( x1-x_obst )/rx ).*( ( x1-x_obst )/rx ) ) )
plot(x1,y)
scatter(xCurr[:,1],xCurr[:,2])
xlabel("x in [m]")
ylabel("y in [m]")
ax = gca()
grid()
#ax[:set_xlim]([3,5])
#ax[:set_ylim]([-1,1])
subplot(212, sharex=ax)
plot(xCurr[1:i_final], costs[1:i_final,1])
plot(xCurr[1:i_final], costs[1:i_final,2])
ax2 = gca()
#ax2[:set_xlim]([-2,12])
xlabel("x in [m]")
ylabel("Costs")
grid()


for i = 1:i_final
    figure(3 ,figsize=(10, 8))   
    clf()
    subplot(511)
    plot(tCurr[1:i_final], xCurr[1:i_final,1])
    plot(tCurr[i:i+10], x_log[:,1,i] ,marker="o")
    xlabel("t in [s]")
    ylabel("x in [m]")

    subplot(512)
    plot(tCurr[1:i_final], xCurr[1:i_final,2])
    plot(tCurr[i:i+10], x_log[:,2,i] ,marker="o")
    xlabel("t in [s]")
    ylabel("y in [m]")

    subplot(513)
    plot(tCurr[1:i_final], xCurr[1:i_final,3])
    plot(tCurr[i:i+10], x_log[:,3,i] ,marker="o")
    xlabel("t in [s]")
    ylabel("v in [m/s]")

    subplot(514)
    plot(tCurr[1:i_final], uCurr[1:i_final,1])
    plot(tCurr[i:i+9], u_log[:,1,i] ,marker="o", color ="red")
    xlabel("t in [s]")
    ylabel("a in [m/s^2]")

    subplot(515)
    plot(tCurr[1:i_final], uCurr[1:i_final,2])
    plot(tCurr[i:i+9], u_log[:,2,i] ,marker="o", color ="red")
    xlabel("t in [s]")
    ylabel("delta in [rad]")

    println("Press Enter for next plot step")
    println("Press c to cancel plot")
    a = ' '
    a = readline()
    if a == "c\r\n" 
         break
    end
end

end #end Module