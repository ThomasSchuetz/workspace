import design_optimization_fixed_devices as opti_model


def recalc(individual, number):
    tag = "no_restrictions"
    filename_start_values = "start_values_without_enev.csv"

    # Compute limits (min costs, min emissions)
    emissions_max = 1000 # ton CO2 per year
    # Minimize costs
    filename_min_costs = "res_ga_fixed_devices/" + tag + str(number) + ".pkl"
    options={"filename_results" : filename_min_costs,
             "enev_restrictions": False,
             "pv_scenario": False,
             "opt_costs": True,
             "store_start_vals": False,
             "load_start_vals": False,
             "filename_start_vals": filename_start_values}

    components = {"roofs":   {i: individual[i]   for i in range(4)},
                  "walls":   {i: individual[i+4] for i in range(4)},
                  "windows": {i: individual[i+8] for i in range(4)}}

    (min_costs, max_emissions) = opti_model.optimize(emissions_max, options, components)

    print min_costs


if __name__ == "__main__":
    individual_ga = [0,1,0,0,1,0,0,0,1,0,0,0]
    individual_guess = [1,0,0,0,1,0,0,0,1,0,0,0]
    
    recalc(individual_ga, 1)
    recalc(individual_guess, 2)