## Network J_ EGS
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

data_path="C:/Users/c3264323/Documents/GitHub/Network_J_Generation_Shedding/data/NVLFT_J"
cd(data_path) 
dss_string = open(f->read(f, String), "Master.dss")

const _ODSS = OpenDSSDirect
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
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>100)
# https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/#Getting-started-with-JuMP
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

_ODSS.Basic.SetActiveClass("Load")    
name_set=_ODSS.ActiveClass.AllNames()
load_name_set2index=Dict(name_set[k]=>k for k=1:length(name_set))
bus_names = _ODSS.Circuit.AllBusNames()  #Rearranging the bus names from 1 t0 32
number_of_buses = length(bus_names) 

Vmag_Phase_a = zeros(number_of_buses,1000) #Phase voltages
Vmag_Phase_b = zeros(number_of_buses,1000)
Vmag_Phase_c = zeros(number_of_buses,1000)

@show Vsources.AllNames()
Vsources.Name("source")
_ODSS.Basic.SetActiveClass("Load")

# Defining random distributions for generation and demand
α=0.8;
m1=8
std1=1
γ1 = 1+std1^2/m1^2
μ1 = log(m1/sqrt(γ1))
σ1 = sqrt(log(γ1))
Pgeneration1= -rand(LogNormal(μ1,σ1),87,1001 )
Qgeneration1=0.0*ones(87,1001)
# plot(Pgeneration1)

m2=0.5;m3=0.2
std2=0.2;std3=0.05;
γ2 = 1+std2^2/m2^2
γ3 = 1+std3^2/m3^2
μ2 = log(m2/sqrt(γ2))
μ3 = log(m3/sqrt(γ3))
σ2 = sqrt(log(γ2))
σ3 = sqrt(log(γ3))
Pdemand1= rand(LogNormal(μ2,σ2),87,1001 )
Qdemand1= rand(LogNormal(μ3,σ3),87,1001)
# plot(Pdemand1)


load_profile =Pdemand1 
all_load_profiles =Dict()

for ld in all_loads
    all_load_profiles[ld] = Pdemand1
end
TT = 1000
load_kW =Pdemand1
PV_kW =Pgeneration1 
net_load_kW= PV_kW + load_kW 

load_voltage_dss = zeros(87,1000) # voltage obtained by solving PF equation using _ODSS
load_voltage_filter = zeros(87,1000) # voltage after passing through low-pass filter
load_voltage_filter[:,1]=230*ones(87,1); 

## MAIN LOOP
for i = 2:TT
 # defining the inverter at all load buses
    _ODSS.Basic.SetActiveClass("Load")
    for ld in all_loads
        _ODSS.Loads.Name(ld)
        load_kW[parse(Int64,ld),i] = Pdemand1[parse(Int64,ld),i]
        PV_kW[parse(Int64,ld),i] = Pgeneration1[parse(Int64,ld),i]
        net_load_kW[parse(Int64,ld),i] = PV_kW[parse(Int64,ld),i]+ load_kW[parse(Int64,ld),i]
    end
    
    for ld in all_loads
         _ODSS.Properties.Value("kW",string(load_kW[parse(Int64,ld),i]+PV_kW[parse(Int64,ld),i]))
    end

    #Rising slack bus voltage
    _ODSS.Basic.SetActiveClass("Load")
    _ODSS.Loads.AllNames()
    if  i>=200
       _ODSS.Vsources.PU(0.92*1.07) #changing V0 to be 230v and raising the slack bus voltage  
     else
        _ODSS.Vsources.PU(0.92)
    end
    #Solution.Solve()
    _ODSS.dss("solve") # Solving pf equations by _ODSS

     # The for all 87 loads at different phases
    load_volt =  Dict{String, Any}(ld => Any[] for ld in all_loads)
    for ld in all_loads
        _ODSS.Circuit.SetActiveElement(ld)
        v_ld = _ODSS.CktElement.VoltagesMagAng()
        append!(load_volt[ld], [v_ld[1,1]])
        
        load_voltage_dss[load_name_set2index[ld], i] = Array{Float64}(load_volt[ld])[1,1]
    end
    #passing tvoltages throgh low-pass filter
    for ld in all_loads
        load_voltage_filter[load_name_set2index[ld], i] = α* load_voltage_filter[load_name_set2index[ld], i-1]+ (1-α)*load_voltage_dss[load_name_set2index[ld], i]
    end
    #Running power flow
    bus_dict = get_bus_voltage_snap()
    branch_dict = get_branch_flows_snap()
    V_source = Vsources.PU()
    # @show V_source
    
   # voltages at each phase
    Vmag_Phase_a[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(1)
    Vmag_Phase_b[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(2)
    Vmag_Phase_c[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(3)

    # Disconnection conditions of inverters 1) if V_mag_average for 5 Minutes > 258v -> disconnect inverter and keep it disconnected for 50 seconds
                                    #       2) if V_mag >265 -> disconnect inverter and keep it disconnected for 50 seconds
                                    # ONLY reconnect inverter if 216v < Vmag <253v for 60 seconds
     
    Vmag_ave= zeros(87,1);
    for ld in all_loads
        _ODSS.Loads.Name(ld)

       # Inverter Disconnection condition 1
        if i> 301
            Vmag_ave[parse(Int64,ld)]= sum(load_voltage_filter[load_name_set2index[ld],i-300:i])/300
            if Vmag_ave[parse(Int64,ld)] .> 255 # it should be 258V
                print("true")
                if i<949
                PV_kW[parse(Int64,ld),i+2:i+50] .= 0
                end
            elseif all(load_voltage_filter[load_name_set2index[ld],i-60:i] .> 216 )  &&  all(load_voltage_filter[load_name_set2index[ld],i-60:i] .< 253 )
                   PV_kW[parse(Int64,ld),i+1] = Pgeneration1[parse(Int64,ld),i+1] 
            end
        end

         #  Inverter Disconnection condition 2
        if i> 61 
            if load_voltage_filter[load_name_set2index[ld],i] .> 256 #It should be 265v
                 print("true")
                 if i<949
                 PV_kW[parse(Int64,ld),i+1:i+50].= 0
                 end
            elseif all(load_voltage_filter[load_name_set2index[ld],i-60:i] .> 216) &&  all(load_voltage_filter[load_name_set2index[ld],i-60:i] .< 253)
                   PV_kW[parse(Int64,ld),i+1] = Pgeneration1[parse(Int64,ld),i+1]
            end
        end
        
    end
end

plot(250*Vmag_Phase_a[:,2:1000]') # Phase voltages
plot(250*Vmag_Phase_b[:,2:1000]')
plot(250*Vmag_Phase_c[:,2:1000]')

plot(load_kW[:,2:1000]') 
plot(net_load_kW[:,2:1000]')
Vmag_plot = plot(load_voltage_filter[:,2:1000]') # This is the voltage we want to measure
savefig(Vmag_plot,"Vmag.png")
plot(load_voltage_dss[:,2:1000]')
plot(PV_kW[:,2:1000]')
