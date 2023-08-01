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
# for bus_name in _ODSS.Circuit.AllBusNames()
#     _ODSS.Circuit.SetActiveBus(bus_name)
    
#     loads = _ODSS.Bus.LoadList()
#     for pf in loads
#         if pf<0
#             pf=-pf
#         else
#             pf=pf
#         end
#     end
# end

# loads_value = branch_dict[Lines][pf]
# print(branch_dict"[Lines(Any))]")
# print( _ODSS.PDElements.AllPowers())

m=4.5
std=0.5
function myLogNormal(m,std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))

    return LogNormal(μ,σ)
end
plot(myLogNormal(m,std))



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
for ld in all_loads
    for (bus, components) in bus_components
        if ld ∈ components["load"]
            append!(loads[ld], [bus])
        end
    end
end
# print(_ODSS.Loads.AllNames())
# loads = Dict(ld => Any[] for ld in all_loads)

bus_names = _ODSS.Circuit.AllBusNames()  #Bus names
number_of_buses = length(bus_names) 
Vmag_Phase_a = zeros(number_of_buses,300)
Vmag_Phase_b = zeros(number_of_buses,300)
Vmag_Phase_c = zeros(number_of_buses,300)
#first make all the positive loads negative, then raise the voltage



@show Vsources.AllNames()
Vsources.Name("source")
_ODSS.Basic.SetActiveClass("Load")
# Loads.Name("1")
load_profile = rand(1)*ones(1,300) #[0.7 0.75 0.68 0.7 0.76 0.89 0.8 0.92 0.95 0.9]
all_load_profiles =Dict()

#defining a random distribution for loads
for ld in all_loads
    all_load_profiles[ld] = 1.8*rand(1)
end

load_kW = zeros(87,300)
PV_kW = rand(87,300)

#
# MAIN LOOP
for i = 1:300
 # defining the inverter at all load buses
    
    _ODSS.Basic.SetActiveClass("Load")
    for ld in all_loads
        # _ODSS.ActiveClass.Name(ld)
        _ODSS.Loads.Name(ld)
        load_kW[parse(Int64,ld),i] = _ODSS.Loads.kW()
        # @show load_kW
        if load_kW[parse(Int64,ld),i] >0 #if load>0 make it negative
            _ODSS.Properties.Value("kW",string(-load_profile[i]*all_load_profiles[ld]*load_kW[parse(Int64,ld),i])) #make the load changing with time i
            # load_kW[ld] = _ODSS.Loads.kW()
        else
            _ODSS.Properties.Value("kW",string(load_profile[i]*all_load_profiles[ld]*load_kW[parse(Int64,ld),i]))
        end
        load_kW[parse(Int64,ld),i] = _ODSS.Loads.kW()
    end
    #@show i
    @show load_kW
    Solution.Solve()

    #Rising slack bus voltage
    _ODSS.Basic.SetActiveClass("Load")
    _ODSS.Loads.AllNames()
    if  i>=50
       _ODSS.Vsources.PU(1.05)
     else
        _ODSS.Vsources.PU(1.0)
    end

    # V_source =_ODSS.Vsources.PU()
    # @show V_source
    
    #Running power flow
    bus_dict = get_bus_voltage_snap()
    branch_dict = get_branch_flows_snap()
    V_source = Vsources.PU()
    @show V_source
    # _ODSS.Basic.SetActiveClass("Load")
    
    Loads.Name("1")
    @show CktElement.Name(), real(CktElement.TotalPowers())
    Vmag_Phase_a[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(1)
    Vmag_Phase_b[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(2)
    Vmag_Phase_c[:,i] = _ODSS.Circuit.AllNodeVmagPUByPhase(3)

    for ld in all_loads
        _ODSS.Loads.Name(ld)
        if i> 120
            Vmag_a_ave= sum(Vmag_Phase_a[i-120:i])/120
            if all(Vmag_a_ave .> 1.122) #258v
                load_kW[parse(Int64,ld),i]=0
            end
        end

        if i> 60
            if Vmag_Phase_a[i-60:i] in 0.939:0.00001:1.1  #[216, 253]V
                load_kW[parse(Int64,ld),i] = load_kW[parse(Int64,ld),i]
            end
        end


        if i>2
            Vmag_Phase_a[:,i-2] .> 1.152  #265v
                load_kW[parse(Int64,ld),i]=0
            

        end

        # if Vmag_Phase_a[:,i] .> 1.196 && i>1    #275v
        #     load_kW[parse(Int64,ld),i][i-0.8]=0
        # end
    end






   


end

plot(Vmag_Phase_a[:,:]')
plot(Vmag_Phase_b[:,:]')
plot(Vmag_Phase_c[:,:]')

plot(load_kW[:,:]')

for i in keys(bus_dict) 
    println(i)

  
end

print(OpenDSSDirect.Vsources.PU())
_ODSS.CktElement.Powers()
_ODSS.PDElements.AllPowers()
ww=ones(1,3)
mat=rand(1)*ones(1,3)
#

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
pp