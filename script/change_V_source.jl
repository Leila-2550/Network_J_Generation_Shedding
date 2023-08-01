using OpenDSSDirect

# Load a circuit somehow -- for conciseness we create a one-load circuit head.
dss("""
    clear
    new circuit.test_circuit bus1=1 BasekV=10
    new load.test_load phases=3 bus1=1 kV=10 kW=1000
    set mode=snapshot
""")

# Display the list of sources available
@show Vsources.AllNames()

# Activate the target source:
# You can either use the name...
Vsources.Name("source")
# Or the index of the object like
# Vsources.Idx(1)

# And activate the load too. 
# This will also activate it as the current circuit element, be careful when
# manipulating multiple elements and using different interfaces!
Loads.Name("test_load")

for i=1:10
    if i >= 5
        # Note that if you only have one Vsource and left it active out of the loop,
        # you don't need to reactivate it here. When in doubt, it's fine to reactivate.
        Vsources.PU(1.2)
    else
        Vsources.PU(1.0)
    end
    # Since we have set the mode to snapshot, let's solve a snapshot here
    Solution.Solve()

    # Some properties are updated only after a Solve. It's not the case for Vsources.PU(),
    # but it's a good habit to avoid checking before solving.
    V_source = Vsources.PU()
    @show i
    @show V_source
    # Check the load status
    @show CktElement.Name(), real(CktElement.TotalPowers())
    @show _ODSS.Circuit.AllNodeVmagPUByPhase(1)
end


w=[1, 2, 3, 4, 5, 6 , 6, 7, 8, 9, 10]
