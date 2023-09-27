## Running  original network j power flow without considering delay time off for the inverters
using Pkg
Pkg.activate("./")
using OpenDSSDirect
import OpenDSSDirect
using Ipopt
using JuMP
#using PowerModelsDistribution
using Plots
using Distributions
using StatsPlots
using DataFrames

data_path="C:/Users/c3264323/Documents/GitHub/Network_J_Generation_Shedding/data/NVLFT_J"
cd(data_path) # when there is an error here, check the names and address in data_path
dss_string = open(f->read(f, String), "Master.dss")
#OpenDSSDirect.dss(dss_string)
const _ODSS = OpenDSSDirect
dss(dss_string)
function get_bus_voltage_snap()
    bus_dict = Dict() # Creating an Empty dictionary with name bus_dict
    phase_mapping = Dict(1=>"a", 2=>"b", 3=>"c", 4=>"n");#key,value Dictionary_name = Dict{Key_datatype, Value_datatype}("Key1" => value1, ...)
    for bus_name in OpenDSSDirect.Circuit.AllBusNames()
        bus_dict[bus_name]  = Dict() #Dictionary_name[key_name] # define another dict in bus_dict
        OpenDSSDirect.Circuit.SetActiveBus(bus_name)
        
        bus_phases = OpenDSSDirect.Bus.Nodes() #3-element Vector{Int64}: 1 2 3 
        bus_voltages = OpenDSSDirect.Bus.PuVoltage()# 3-element Vector{ComplexF64}: 3 voltage magnitudes with real and im part for each
        if 4 ∈ bus_phases
            filter!(e->e≠4, bus_phases) #filter(function, collection)
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
        # branch_dict[branch_name]["phases"] = branch_phases
    end
    return branch_dict
end
branch_dict = get_branch_flows_snap()# give values of powers for lines
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>100)
# https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/#Getting-started-with-JuMP
all_loads = _ODSS.Loads.AllNames() #what are our loads
loads = Dict(ld => Any[] for ld in all_loads) #generating a dic that contain loads
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

#determining the loads and their buses
_ODSS.Basic.SetActiveClass("Load")    
##determining the loads and their buses
_ODSS.Basic.SetActiveClass("Load")    #
name_set=_ODSS.ActiveClass.AllNames()
load_name_set2index=Dict(name_set[k]=>k for k=1:length(name_set))

## computing the load voltages out of the loop
# load_voltage=zeros(87,1000)
# for TT=1:1000
#     load_volt =  Dict{String, Any}(ld => Any[] for ld in all_loads)
#     for ld in all_loads
#         _ODSS.Circuit.SetActiveElement(ld)
#         v_ld = _ODSS.CktElement.VoltagesMagAng()
#         # @show v_ld
#         append!(load_volt[ld], [v_ld[1,1]])
        
#         load_voltage[load_name_set2index[ld], TT] = Array{Float64}(load_volt[ld])[1,1]
#         #df_load_voltage[df_load_Voltage."Load" ==ld,"t = $TT"] = load_volt[ld]
#     end
# end
##
bus_names = _ODSS.Circuit.AllBusNames()  #Rearranging the bus names from 1 t0 32
number_of_buses = length(bus_names) 

Vmag_Phase_a = zeros(number_of_buses,1000)
Vmag_Phase_b = zeros(number_of_buses,1000)
Vmag_Phase_c = zeros(number_of_buses,1000)

@show Vsources.AllNames()
Vsources.Name("source")
_ODSS.Basic.SetActiveClass("Load")
α=0.8; # Low-pass filter parameter
# Defining Distributions for generation and demand
m1=8
std1=1
γ1 = 1+std1^2/m1^2
μ1 = log(m1/sqrt(γ1))
σ1 = sqrt(log(γ1))
Pgeneration1= -rand(LogNormal(μ1,σ1),87,1001)
Qgeneration1=0.0*ones(87,1001)
# plot(Pgeneration1)
# plot(Qgeneration1)

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
# plot(Qdemand1)

load_profile =Pdemand1 
all_load_profiles =Dict()
#defining a random distribution for loads
for ld in all_loads
    all_load_profiles[ld] = Pdemand1
end
TT = 1000
load_kW = Pdemand1
PV_kW = Pgeneration1 
net_load_kW= PV_kW + load_kW # P=P_generation+P_demand

load_voltage_dss = zeros(87,1000)
load_voltage_filter = zeros(87,1000)
load_voltage_filter[:,1]=230*ones(87,1);

## MAIN LOOP
# Running the simulations for 1000 sec

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

    #Rising the slack bus voltage
    _ODSS.Basic.SetActiveClass("Load")
    _ODSS.Loads.AllNames()
    if  i>=200
       _ODSS.Vsources.PU(0.92*1.1) # 0.06,250. 0.07 255
     else
        _ODSS.Vsources.PU(0.92)
    end
    #Solution.Solve()
    _ODSS.dss("solve") # opendss solves the power flow equations

    # Computing the load voltages 
    load_volt =  Dict{String, Any}(ld => Any[] for ld in all_loads)
    for ld in all_loads
        _ODSS.Circuit.SetActiveElement(ld)
        v_ld = _ODSS.CktElement.VoltagesMagAng()
        append!(load_volt[ld], [v_ld[1,1]])
        load_voltage_dss[load_name_set2index[ld], i] = Array{Float64}(load_volt[ld])[1,1]
    end

    #Passing the voltages through low-pass filter
    for ld in all_loads
        load_voltage_filter[load_name_set2index[ld], i] = α* load_voltage_filter[load_name_set2index[ld], i-1]+ (1-α)*load_voltage_dss[load_name_set2index[ld], i]
    end

    #Running power flow
    bus_dict = get_bus_voltage_snap()
    branch_dict = get_branch_flows_snap()
    V_source = Vsources.PU()
     
    Vmag_Phase_a[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(1)
    Vmag_Phase_b[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(2)
    Vmag_Phase_c[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(3)

# Disconnection conditions of inverters 1) if V_mag_average for 5 Minutes > 258v --> disconnect inverter and keep it disconnected for 50 seconds
                                    #   2) if V_mag >265 --> disconnect inverter and keep it disconnected for 50 seconds
                                    # ONLY reconnect inverter if 216v < Vmag <253v for 60 seconds
    #Condition 1                                
    Vmag_ave= zeros(87,1);
    for ld in all_loads
        _ODSS.Loads.Name(ld)
        if i> 301
            Vmag_ave[parse(Int64,ld)]= sum(load_voltage_filter[load_name_set2index[ld],i-300:i])/300
            if Vmag_ave[parse(Int64,ld)] .> 258 #258V
                print("true")
                if i<949
                PV_kW[parse(Int64,ld),i+1:i+50] .= 0
                end
            elseif all(load_voltage_filter[load_name_set2index[ld],i-60:i] .> 216 )  &&  all(load_voltage_filter[load_name_set2index[ld],i-60:i] .< 253 )
                   PV_kW[parse(Int64,ld),i+1] = Pgeneration1[parse(Int64,ld),i+1] 
            end
        end

        #Condition 2     
        if i> 61 
            if load_voltage_filter[load_name_set2index[ld],i] .> 258 # it should be 265
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

plot(250*Vmag_Phase_a[:,2:1000]') #
plot(250*Vmag_Phase_b[:,2:1000]')
plot(250*Vmag_Phase_c[:,2:1000]')

plot(load_kW[:,2:1000]')

plot(net_load_kW[:,2:1000]')
plot(load_voltage_filter[:,2:1000]') # Reverse buses 5 8 16 18 19 26 28 33 36 42 50 62 65 67 71 72 73 76 77 78 79 80 83
plot(load_voltage_dss[10:40,2:1000]')

plot(PV_kW[:,2:1000]')

##
for i in keys(bus_dict) 
    println(i)
  
end
print(OpenDSSDirect.Vsources.PU())
_ODSS.CktElement.Powers()
_ODSS.PDElements.AllPowers()
ww=ones(1,3)
mat=rand(1)*ones(1,3)
##


vv = [245 248 248 249 250 257 258 259 260 258 259 255 261 259 263 258 259 256 257 259]
pp = 10*ones(20,1)
for i = 10:20
    ave_vv = sum(vv[i-9:i])/10
    # if i > 9
        if  all(ave_vv .> 257)
            pp[i]=0
        end
    # end
end

m=4.5
std=0.5
function myLogNormal(m,std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))
    return LogNormal(μ,σ)
end
plot(myLogNormal(m,std))