using Pkg
Pkg.activate("./")
using OpenDSSDirect
import OpenDSSDirect
using Ipopt
using JuMP
using PowerModelsDistribution
using Plots

data_path="C:/Users/c3264323/Documents/GitHub/Network_J_Generation_Shedding/data/NVLFT_J"

cd(data_path) # when there is an error here, check the names and address in data_path
dss_string = open(f->read(f, String), "Master.dss")
#OpenDSSDirect.dss(dss_string)


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
            filter!(e->e≠4, bus_phases)#it only gives the values that are correct for the first condition
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

    for branch_name in OpenDSSDirect.PDElements.AllNames()#Power Delivery:capacitor, transformer,line
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














## thsi section uses PowerModelsDistribution, which does not support inverter defining

# data_eng = parse_file("MasterAll.dss", transformations=[transform_loops!])

# data_math = transform_data_model(data_eng; kron_reduce=true)

# results_pmd = solve_mc_pf(data_math, ACPUPowerModel, solver)


# ##
# vbase = data_math["bus"]["1"]["vbase"]
# bus_voltages_pmd = []
# bus_voltages_odss = []

# for (i, bus_pmd) in results_pmd["solution"]["bus"]
#     bus_name = data_math["bus"][i]["name"]

#     if !occursin("sourcebus", bus_name) && !occursin("transformer", bus_name) && !occursin("source", bus_name)
#         append!(bus_voltages_pmd, bus_pmd["vm"])

#         @show bus_name
#         bus_odss = bus_dict[bus_name]
#         bus_odss_keys = [k for (k,v) in bus_odss]
#         vm_phases = bus_odss_keys[occursin.("m", bus_odss_keys)]
#         append!(bus_voltages_odss, [bus_odss[v] for v in vm_phases])
#     end
# end

# plot(bus_voltages_pmd, seriestype=:scatter, label="PMD")
# plot!(bus_voltages_odss, seriestype=:scatter, label="DSS")
