#==========

Builds model costs from physical and financial assumptions, so they stay consistent when a fuel
price or the discount rate changes. Units: capex $/kW, fixed O&M $/kW-yr, heat rate MMBtu/MWh,
fuel price $/MMBtu, variable O&M $/MWh. Results are $/MW-yr and $/MWh.

==========#

# Share of capex recovered each year by a level annuity over `life_years`
function capital_recovery_factor(discount_rate::Real, life_years::Real)
    life_years > 0 || throw(ArgumentError("Life must be positive, got $life_years."))
    discount_rate >= 0 || throw(ArgumentError("Discount rate can't be negative, got $discount_rate."))
    discount_rate == 0 && return 1 / life_years
    return discount_rate / (1 - (1 + discount_rate)^(-life_years))
end

# $/MW-yr. The factor of 1000 converts $/kW to $/MW.
function annual_investment_cost(capex_per_kw::Real, fixed_om_per_kw_yr::Real, life_years::Real, discount_rate::Real; multiplier::Real = 1.0)
    return 1000 * multiplier * (capex_per_kw * capital_recovery_factor(discount_rate, life_years) + fixed_om_per_kw_yr)
end

# $/MWh
variable_cost(heat_rate::Real, fuel_price::Real, variable_om::Real) = heat_rate * fuel_price + variable_om

# ton CO2/MWh
emission_intensity(heat_rate::Real, carbon_per_mmbtu::Real) = heat_rate * carbon_per_mmbtu
