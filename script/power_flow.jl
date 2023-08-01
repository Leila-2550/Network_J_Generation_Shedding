using Pkg
Pkg.activate("./")

using OpenDSSDirect
import OpenDSSDirect
using Ipopt
using JuMP
using PowerModelsDistribution
using Plots
using Distributions
using StatsPlots
data_path="C:/Users/c3264323/Documents/GitHub/Network_J_Generation_Shedding/data/NVLFT_J"

cd(data_path) # when there is an error here, check the names and address in data_path
dss_string = open(f->read(f, String), "Master.dss")
#OpenDSSDirect.dss(dss_string)
const _ODSS = OpenDSSDirect

dss(dss_string)

base_kV1 = _ODSS.Bus.kVBase()
base_kV2 = _ODSS.Bus.kVBase()*1.1

m=0.75
std=0.1
function myLogNormal(m,std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))

    return LogNormal(μ,σ)
end
# plot(myLogNormal(m,std))

# variable_loads_buses = []
# bus_names = _ODSS.Circuit.AllBusNames()  #Bus names
# number_of_buses = length(bus_names) 
# number_of_variable_loads_buses = length(variable_loads_buses)
load_kW = _ODSS.Loads.kW()

bus_names = _ODSS.Circuit.AllBusNames()  #Bus names
number_of_buses = length(bus_names) 
bus_names_to_index = Dict(bus_names[k]=>k for k=1:length(bus_names))   #Indexing the buses names
bus_index_to_names = Dict(k=>[k] for k=1:length(bus_names)) 


#defining the inverter's Network_J_Generation_Shedding

all_loads = _ODSS.Loads.AllNames()
# print(_ODSS.Loads.AllNames())
loads = Dict(ld => Any[] for ld in all_loads)

for ld in all_loads
#    myLogNormal(m,std)  =_ODSS.Loads.kW()
load_kW = _ODSS.Loads.kW()
   @show load_kW
end


for i =1:200

    Vmax, Vmin = 1.05,0.95
    V_source = 1.00

    _ODSS.Bus.kVBase = base_kV1   
        load_profile_kW = append!(loads_profile_kW,_ODSS.Loads.kW()*loadshape_mult[i])
    end
end



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

bus_dict = get_bus_voltage_snap()

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

branch_dict = get_branch_flows_snap()

solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>100)







