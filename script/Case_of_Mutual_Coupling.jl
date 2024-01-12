## Network J_Mutual coupling/ Sensitivity of voltage and power
import Pkg
Pkg.add("Ipopt")
using Pkg
Pkg.activate("./")
using OpenDSSDirect
import OpenDSSDirect
using Ipopt
using JuMP
using Plots
using Distributions
using StatsPlots
using DataFrames
const _ODSS = OpenDSSDirect

data_path="/Users/raj055/Documents/Julia files/Network_J-main/data/NVLFT_J"
cd(data_path) 
dss_string = open(f->read(f, String), "Master.dss")
dss(dss_string)

function get_bus_voltage_snap()
    bus_dict = Dict() 
    phase_mapping = Dict(1=>"a", 2=>"b", 3=>"c", 4=>"n");
    for bus_name in OpenDSSDirect.Circuit.AllBusNames()
        bus_dict[bus_name]  = Dict() 
        OpenDSSDirect.Circuit.SetActiveBus(bus_name)
        
        bus_phases = OpenDSSDirect.Bus.Nodes() 
        bus_voltages = OpenDSSDirect.Bus.PuVoltage()
        if 4 ∈ bus_phases
            filter!(e->e≠4, bus_phases) 
            for (i, phase) in enumerate(bus_phases) 
                bus_dict[bus_name]["vm"*phase_mapping[phase]] = abs(bus_voltages[i] - bus_voltages[4])
                bus_dict[bus_name]["va"*phase_mapping[phase]] = angle(bus_voltages[i] - bus_voltages[4]) * 180/pi
            end
        else
            for (i, phase) in enumerate(bus_phases) 
                bus_dict[bus_name]["vm"*phase_mapping[phase]] = abs(bus_voltages[i])
                bus_dict[bus_name]["va"*phase_mapping[phase]] = angle(bus_voltages[i]) * 180/pi
            end
        end  
    end
    return bus_dict
end
bus_dict = get_bus_voltage_snap()# give the voltages of all phases of all buses

vma = [bus["vma"] for (i,bus) in bus_dict]
vmb = [bus["vmb"] for (i,bus) in bus_dict]
vmc = [bus["vmc"] for (i,bus) in bus_dict]

function get_branch_flows_snap()
    branch_dict = Dict()
    phase_mapping = Dict(1=>"a", 2=>"b", 3=>"c", 4=>"n");
    for branch_name in OpenDSSDirect.PDElements.AllNames()
        branch_dict[branch_name]  = Dict()
        OpenDSSDirect.PDElements.Name(branch_name)
        
        nphases = OpenDSSDirect.CktElement.NumPhases()
        nwires = Int(length(OpenDSSDirect.CktElement.NodeOrder())/2)
        branch_flows = OpenDSSDirect.CktElement.Powers()
        branch_losses = OpenDSSDirect.CktElement.PhaseLosses()
        branch_dict[branch_name]["pf"] = real.(branch_flows[1:nwires])
        branch_dict[branch_name]["qf"] = imag.(branch_flows[1:nwires])
        branch_dict[branch_name]["pt"] = real.(branch_flows[nwires+1:end])
        branch_dict[branch_name]["qt"] = imag.(branch_flows[nwires+1:end])
        branch_dict[branch_name]["ploss"] = real.(branch_losses)
        branch_dict[branch_name]["qloss"] = imag.(branch_losses)
    end

    return branch_dict
end
branch_dict = get_branch_flows_snap()

all_loads = _ODSS.Loads.AllNames()
loads = Dict(ld => Any[] for ld in all_loads) 

## Function: Get Bus Components: Lines, Transformers, Loads
function get_bus_components(;comp_types=["line", "transformer", "load", "gen"])::Dict{String,Any}
    bus_components = Dict{String,Any}()
    for bus_name in _ODSS.Circuit.AllBusNames()
        _ODSS.Circuit.SetActiveBus(bus_name)

        pde_elements = _ODSS.Bus.AllPDEatBus()
        loads = _ODSS.Bus.LoadList()

        bus_components[bus_name] = Dict(comp_type => Any[] for comp_type in comp_types)
        for pde_element in pde_elements
            if startswith(pde_element, "Line")
                append!(bus_components[bus_name]["line"], split(pde_element,".")[2:end])
            elseif startswith(pde_element, "Transformer")
                append!(bus_components[bus_name]["transformer"], split(pde_element,".")[2:end])
            end
        end
        for load in loads
            append!(bus_components[bus_name]["load"], split(load,".")[2:end])
        end
    end
    return bus_components
end

# Get Bus Components 
bus_components = get_bus_components()
 
_ODSS.Basic.SetActiveClass("Load")    
name_set=_ODSS.ActiveClass.AllNames()
load_name_set2index=Dict(name_set[k]=>k for k=1:length(name_set))

bus_names = _ODSS.Circuit.AllBusNames()  #Rearranging the bus names from 1 t0 32
number_of_buses = length(bus_names) 

Vmag_Phase_a = zeros(number_of_buses,2001); # Phase voltages
Vmag_Phase_b = zeros(number_of_buses,2001);
Vmag_Phase_c = zeros(number_of_buses,2001);

@show Vsources.AllNames()
Vsources.Name("source")
_ODSS.Basic.SetActiveClass("Load")


## Defining distributions of generation and demand

α=0.8; #low-pass filter coefficient

# load_phaseA=[0 1.0 1.0 1 0 1 1 1 0 1 0 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 0 0 1 1 1 1 0  0 1 0 0 1 0 1 0 0  0 1 #=
# =# 1 0 1 1 0 1 0 1 1 1 0 0 1 0 1 0 0 1 0 1 1 0 1 0 0 0 1 1 1 0 1 1  0  1 1 1 0 0 1 0 0 0 0];
# Pgeneration1= -6 *load_phaseA'* ones(1,5501)
# n1=findall(x->x==1,load_phaseA)

# Pgeneration1= -7.1 * ones(87,2001) # when having load at tall phases
# Pgeneration1[16,:] = -8.8 * ones(1,2001);
# Pgeneration1[71,:] = -8.8 * ones(1,2001);

# # when having load at Oseperate diffrenet phases
load_phaseA=[0 1.0 1.0 1 0 1 1 1 0 1 0 0 1 0 0 1.5 0 1 1 0 0 0 1 0 0 1 0 0 1 1 1 1 0  0 1 0 0 1 0 1 0 0  0 1 #=
=# 1 0 1 1 0 1 0 1 1 1 0 0 1 0 1 0 0 1 0 1 1 0 1 0 0 0 1.5 1 1 0 1 1  0  1 1 1 0 0 1 0 0 0 0];
# # load_phaseA[1,59]=3;load_phaseA[1,75]=3;load_phaseA[1,67]=3;
n1=findall(x->x==1,load_phaseA) #41 element
# # Pgeneration1 = -5.5 *load_phaseA'* ones(1,5501)
load_phaseB=[ 1.0 0 0 0 0 0 0 0 0 0 1.0 1 0 1 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0  0 0 0 1 0 1 0 0 0 1 0 #=
=# 0  0 0 0 1 0 0 0 0 0 1 1 0 1 0 0 1 0 1 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 0 0];
# # load_phaseB[1,49]=3;load_phaseB[1,14]=3;load_phaseB[1,66]=3;
n2=findall(x->x==1,load_phaseB) #23 element
load_phaseC=[ 0 0 0 0 1 0 0  0 1.0 0 0 0 0 0 1 0 1 0 0 0 1.0 1 0 0 0  0 1 1 0 0 0 0 1 1  0 1 0 0 0 0 1 1 0 0 #=
=# 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1.0 0 0 1 0 0 0 0 0 0 0 0 1 1];
# # load_phaseC[1,17]=3;
n3=findall(x->x==1,load_phaseC) # 21 element
# Pgeneration1 =-0.9767 *load_phaseA'* ones(1,5501) + -1.8261 * load_phaseB' *ones(1,5501) + -2* load_phaseC' *ones(1,5501);#having equal loads at all phases
Pgeneration1 =-2*load_phaseA'* ones(1,2001) + -2.021 * load_phaseB' *ones(1,2001) + -2.875* load_phaseC' *ones(1,2001);#having equal loads at all phases
# Pgeneration1[18,:] = -9 * ones(1,2001);
# Pgeneration1[80,:] = -9 * ones(1,2001);
# Pgeneration1[79,:] = -9 * ones(1,2001);

Qgeneration1 =  0.0 * ones(87,2001)
# plot(Pgeneration1)


Pdemand1 = 2.1 * ones(87,2001) # when having load at all phases

load_kW = zeros(87,2001); #Pdemand1
PV_kW = zeros(87,2001); #Inverters export
net_load_kW = PV_kW + load_kW ;

load_voltage_dss = zeros(87,2000); # voltage obtained by solving PF equation using _ODSS
load_voltage_filter = zeros(87,2000);# voltage after passing through low-pass filter
load_voltage_filter[:,1]=230.94*ones(87,1); 

node_power = zeros(87,2000); # power at each node
node_power[:,1]=zeros(87,1);


## MAIN LOOP
for i = 2:2000
 # defining the inverter at all load buses
    _ODSS.Basic.SetActiveClass("Load")
    for ld in all_loads
        _ODSS.Loads.Name(ld)
        load_kW[parse(Int64,ld),i] = copy(Pdemand1[parse(Int64,ld),i])
        PV_kW[parse(Int64,ld),i] = copy(Pgeneration1[parse(Int64,ld),i])
        net_load_kW[parse(Int64,ld),i] = PV_kW[parse(Int64,ld),i]+ load_kW[parse(Int64,ld),i]

        #  _ODSS.Properties.Value("kW",string(net_load_kW[parse(Int64,ld),i] ))
        _ODSS.Loads.kW(net_load_kW[parse(Int64,ld),i])
    end
    

    #Rising slack bus voltage
    _ODSS.Basic.SetActiveClass("Load")
    # _ODSS.Loads.AllNames()
    # if  i>=100 && i<= 1000
    #     _ODSS.Vsources.PU(1.1) #changing V0 to be 230v and raising the slack bus voltage  
    # elseif i>=1001 && i<= 1900
    #     _ODSS.Vsources.PU(1.1) 
    # elseif i>=1901 && i<= 5500
    #     _ODSS.Vsources.PU(1.1)
    #  else
    #     _ODSS.Vsources.PU(1.0)
    # end

    if i>=1000
        _ODSS.Vsources.PU(1.01)
    else
        _ODSS.Vsources.PU(1.01)
    end
    Solution.Solve()
    _ODSS.dss("solve") # Solving pf equations by _ODSS

     # The for all 87 loads at different phases
    load_volt =  Dict{String, Any}(ld => Any[] for ld in all_loads)
    for ld in all_loads
        _ODSS.Loads.Name(ld)
        load_bus = [i for (i, comps) in bus_components if ld in comps["load"]][1]
        # v_ld = _ODSS.CktElement.VoltagesMagAng()
        _ODSS.Circuit.SetActiveBus(load_bus)
        bus_voltages = _ODSS.Bus.puVmagAngle()
        vbus_mag = bus_voltages[1:2:6] .* 230.94
        vbus_ang = bus_voltages[2:2:6]

        load_phases = _ODSS.Loads.Phases()
        

        append!(load_volt[ld], vbus_mag[load_phases])
        
        # load_voltage_dss[load_name_set2index[ld], i] = Array{Float64}(load_volt[ld])[1,1]
        load_voltage_dss[load_name_set2index[ld], i] = vbus_mag[load_phases]

        #passing tvoltages throgh low-pass filter
        load_voltage_filter[load_name_set2index[ld], i] = α* load_voltage_filter[load_name_set2index[ld], i-1]+ (1-α)*load_voltage_dss[load_name_set2index[ld], i]
    # node_power[load_name_set2index[ld],i] = _ODSS.CktElement.Powers()
    end

    #Running power flow
    bus_dict = get_bus_voltage_snap()
    branch_dict = get_branch_flows_snap()
    V_source = Vsources.PU()
    # @show V_source

    # Disconnection conditions of inverters 1) if V_mag_average for 5 Minutes > 258v -> disconnect inverter and keep it disconnected for 50 seconds
                                    #       2) if V_mag >265 -> disconnect inverter and keep it disconnected for 50 seconds
                                    # ONLY reconnect inverter if 216v < Vmag <253v for 60 seconds
     
    Vmag_ave= zeros(87,2000);
    for ld in all_loads
        _ODSS.Loads.Name(ld)

       # Inverter Disconnection condition 1
        if i > 601
            Vmag_ave[parse(Int64,ld)]= sum(load_voltage_filter[load_name_set2index[ld],i-600:i])/600
            if Vmag_ave[parse(Int64,ld)] .> 258 # it should be 258V
                # print("true")
                if i<1439
                PV_kW[parse(Int64,ld),i+2:i+60] .= 0
                end
            # elseif all(load_voltage_filter[load_name_set2index[ld],i-60:i] .> 216 )  &&  all(load_voltage_filter[load_name_set2index[ld],i-60:i] .< 253 )
            #        PV_kW[parse(Int64,ld),i+1] = Pgeneration1[parse(Int64,ld),i+1] 
            end
        end

         #  Inverter Disconnection condition 2
        if i> 61 
            if load_voltage_filter[load_name_set2index[ld],i] .> 265 #It should be 265v
                #  print("true")
                 if i<1439
                 PV_kW[parse(Int64,ld),i+1:i+60].= 0
                 end
            elseif all(load_voltage_filter[load_name_set2index[ld],i-60:i] .> 216) &&  all(load_voltage_filter[load_name_set2index[ld],i-60:i] .< 253)
                   PV_kW[parse(Int64,ld),i+1] = Pgeneration1[parse(Int64,ld),i+1]
            end
        end
        
    end

    # voltages at each phase
    Vmag_Phase_a[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(1)#32*5500
    Vmag_Phase_b[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(2)
    Vmag_Phase_c[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(3)
end

## at all buses
#= P1 = plot(230.94*Vmag_Phase_a[:,1900:5500]',legend = false, ylimits=(242,275)) # Phase voltages
# title!("Phase voltages, random loads at all buses, V0=1.1 p.u")
P2 = plot(230.94*Vmag_Phase_b[:,1900:5500]',legend = false,  ylimits=(242,275))
ylabel!("Vmag(Volt)");
P3 = plot(230.94*Vmag_Phase_c[:,1900:5500]',legend = false,  ylimits=(242,280)) =#

## for only bus 15
V1=plot(230.94*Vmag_Phase_a[2,900:end],legend=false, ylims=(240,250),xlims=(0,1100),xticks=00:200:1100) # Phase voltages
title!("Phase voltage distribution, bus 7")
V2=plot(230.94*Vmag_Phase_b[2,900:end] ,ylims=(237,250),legend=false, label ="Phase B",xlims=(0,1100),xticks=0:200:1100)
ylabel!("Vmag(Volt)");
# V3=plot(230.94*Vmag_Phase_c[2,900:end], legend=false, label ="Phase C",xlims=(0,1100),xticks=0:200:1100)
xlabel!("Time(t)");
V_allphases = [V1,V2]
P4= plot(V1,V2, layout=(2,1))
# display(P4)

#phase power for all loads at one bus 7
PKw1 = plot(net_load_kW[18,900:end], legend=false,xticks=0:200:1100)
PKw2 = plot(net_load_kW[80,900:end], legend=false,xticks=0:200:1100)
P_phaseA = plot(PKw1,PKw2)  # Phase powers
P_phaseB = plot(net_load_kW[81,900:end], legend=false,xticks=0:200:1100)
# P_phaseC = plot(net_load_kW[69,900:end], legend=false,xticks=0:200:1100)
plot(P_phaseA,P_phaseB, layout=(2,1))

Vmag16 = plot(load_voltage_filter[18,900:end], legend = false,xticks=0:200:1100)
Vmag71 = plot(load_voltage_filter[80,900:end], legend = false,xticks=0:200:1100)
V_phaseA = plot(Vmag16,Vmag71) 
V_phaseB = plot(load_voltage_filter[81,900:end], legend = false,xticks=0:200:1100)
# V_phaseC = plot(load_voltage_filter[15,900:end], legend = false,xticks=0:200:1100)
plot(V_phaseA,V_phaseB, layout=(2,1))

plot(net_load_kW[:,900:end]' ,legend=false)
plot(net_load_kW[:,900:end]' * load_phaseA' ,legend=false)
Vmag_plot = plot(load_voltage_filter[:,100:2000]', legend = false) # This is the voltage we want to measure
xlabel!("Time(t)");ylabel!("Vmag(Volt)");title!("Voltage magnitude, No generation")
# Vmag_plot = plot(load_voltage_filter[:,2:5]') # This is the voltage we want to measure
# savefig(Vmag_plot,"Vmag2.png")
Vmag_final = plot(load_voltage_filter[:,101:2000]', legend = false) # This is the voltage we want to measure at only last period
xlabel!("Time(t)");ylabel!("Vmag(Volt)");title!("Voltage at bus 15, loads only at phase A")
# plot(load_voltage_dss[:,2:5500]')
P_PV = PV_kW[:,101:2000]
plot( PV_kW[:,101:2000]',legend=false)
xlabel!("Time(t)");ylabel!("PV_kW");
#  for i in eachindex(P_PV)
#     if P_PV[i]- == 0
#         P_PV[i] = NaN
#     end

    
# load_phaseB[81:87]
PV_kW[:,1901:5500]
plot(PV_kW[4,101:2000])
typeof(PV_kW[:,101:2000])
PV_gen = PV_kW[:,901:2000]#Pgeneration at final period
gen_bb= !=(0).(PV_gen)
gen_diff= diff(gen_bb, dims=2)
count1 = sum(diff(gen_bb, dims=2) .== -1, dims=2)#number of switch of per each node
# count2 = sum(diff(gen_bb, dims=2) .== 1, dims=2)
count11 = sum(diff(gen_bb, dims=2) .== -1) #toatal number of on-off in last time period
# count22 = sum(diff(gen_bb, dims=2) .== 1)

# for i in eachindex(gen_diff[:])
#     if gen_diff[i] == -1
#         @show i
#     end
# end
 
# inv_off=findall(x->x==-1,gen_diff[:,:])

inds = Tuple.(findall(x->x==-1,gen_diff[:,:]))
a = first.(inds)
b = last.(inds)
bb=[a b] #those buses that switched off with the time of switching off
# inv_off=findall(x->x!=0,bb[:,:],dims=2)
no_off_node15 = findall(x->x==15,bb[:,1])
no_off_node79 = findall(x->x==16,bb[:,1])
no_off_node81 = findall(x->x==70,bb[:,1])
no_off_node82= findall(x->x==82,bb[:,1])

no_off_node16= findall(x->x==16,bb[:,1])
no_off_node71 = findall(x->x==71,bb[:,1])
c=zeros(87,3600)
for i in bb[:,1] , j in bb[:,2]
    c[i,j]=1



end


c
plot(c',legend=false)

P_agg=-sum(PV_gen, dims=1)
Pagg_mean = mean(P_agg[:,1:3600], dims=2)
Pagg_plot = plot(P_agg[:,1:3600]')
savefig(Pagg_plot,"P_agg1.png")

plot(Pgeneration1[1:2,1901:5500]',legend=false)
n111=findall(x->x==1,Pgeneration1[24,1901:5500])
# savefig("/Users/raj055/Documents/Julia files-PSCC paper/Network_J_EGS/Network_J_EGS-main/data/NVLFT_J" * "plot_43.png")
# img = load("/Users/raj055/Documents/Julia files-PSCC paper/Network_J_EGS/Network_J_EGS-main/data/NVLFT_Jplot_43.png")
ppv=PV_kW[:,1901:5500]+Pgeneration1[:,1901:5500]
plot(ppv[:,:]',legend=false)
x=1900:5500;
y1=plot(PV_kW[71,1900:5500] ,legend=false,ylabel="(PV_kW)",palette= :Dark2_5,linestyle =:dash)
y2 =load_voltage_filter[:,1900:5500];
y3=plot(Vmag_Phase_a[:,1900:5500]',ylabel="Vmag_Phase_a",xlabel=("Time(s)"));
# plot(x,y1',legend=false, color=:blue, label="Net load")
plot(y1',label = "My y1 label", legend = false, ylabel = "PV_kW on phase B",grid=:off, xlabel="Time", box=:on)
plot!(twinx(),y2', label = "my y2 label", legend=false, ylabel = "Voltage magnitude", grid=:off, xlabel="Time", box=:on,color=:blue, linestyle=:dash)
plot!(twinx(),y3', label = "my y2 label", legend=false, ylabel = "Phase B Voltage magnitude", grid=:off, xlabel="Time", box=:on,color=:blue,linestyle=:dash)


plot(y1,y3,layout=(2,1),legend=false,xlims=(0,3600),xticks=0:300:3600)
_ODSS.Bus.puVmagAngle()
typeof(_ODSS.Circuit.AllBusMagPu).name.mt
_ODSS.CktElement.Powers
_ODSS.Circuit.AllBusMagPu
_ODSS.Circuit.AllElementLosses
_ODSS.Circuit.AllNodeNamesByPhase
