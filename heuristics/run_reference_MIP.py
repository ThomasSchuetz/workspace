import design_computation_5R1C_multi_obj as opti_model

# Compute limits (min costs, min emissions)
emissions_max = 1000 # ton CO2 per year
# Minimize costs
filename_min_costs = "res_ga_fixed_devices/result_mip.pkl"
options={"filename_results" : filename_min_costs,
         "enev_restrictions": False,
         "pv_scenario": False,
         "opt_costs": True,
         "store_start_vals": False,
         "load_start_vals": False,
         "filename_start_vals": "start_values_without_enev.csv"}
(min_costs, max_emissions) = opti_model.optimize(emissions_max, options)