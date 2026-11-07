# Converted from GenX 0.3.7

Source: `.claude/files for claude/Benders_MGA_ScalingFY16z`. Built by `convert_genx_case` (`src/inputs/Convert_genx.jl`).

## Representation
- 12 representative weeks per weather year, including peak demand, peak residual load and lowest renewable weeks. Each zone's annual energy is restored after clustering.
- 5 weather years: scenario labels in `Zone_profiles` show which is the source and which are synthetic. Constraint scores, least to most: [902342.0, 906300.0, 912936.0, 913511.0, 921763.0].
- Synthetic years are block-resampled from the one source year with annual wind and solar anomalies. They aren't observed weather.

## Assumptions added
- Fuel: repo gas distribution, coal 1.47 and uranium 0.85 replace the case's prices (about $14.50 gas, $21 coal). Oil ($38.82) and biomass ($5) keep the case's values.
- Thermal outages: the case's availability is 1, so rates come from `GENX_OUTAGE_RATES`.
- New batteries are fixed at 4 hours (the case allowed 1 to 10), with energy costs included. CRM penalty is the repo's 5000.

## Kept from the case
Annual investment and fixed O&M by zone, heat rates and variable O&M (capacity-weighted across zones per technology), existing fleets, build limits, transmission, demand, availability profiles, NSE costs, CRM margins (flag off, as in the case).

## Availability repaired
The source has wind, solar and distributed resources at exactly 1.0 in every hour, most of them absent from the top-level Generators_variability.csv. These 108 resources (16 Distributed, 1 Offshore existing, 48 Offshore fixed, 12 Offshore floating, 3 Solar, 16 Solar existing, 12 Wind existing) now use the mean of real profiles in their zone, and 24 of them, with none in the zone, use a mean across zones. Fixed offshore mostly borrows floating offshore profiles, since few fixed clusters have real ones. `Availability_repairs.csv` lists each one. They're estimates, and wind and solar tiers are ranked on the repaired values.

Distributed generation is assumed to be solar and uses its zone's new-build solar profile.

## Hydro limits
The source has hydro flat at 1.0. It's now capped at 80% of capacity in any hour (fleet availability, ORNL/DOE U.S. Hydropower Market Report, 2023 edition) and at 38% of full-capacity energy in each representative week (2016-2025 mean capacity factor, EIA Electric Power Monthly Table 6.07.B). The same values apply to every zone and season, and to conventional and small hydro alike, so this has no seasonal or regional shape. The energy limit uses a dispatch result as a stand-in for water and operating limits.

## Not represented
- Dropped resources: Coal_CCS90AvgCF_Moderate, Nuclear_Nuclear_Moderate, comm_space_heat_cool, comm_water_heat, res_space_heat_cool, res_water_heat, trans_light_duty.
- Unit commitment, ramping and minimum power, line losses, zone-specific fuel prices, storage variable O&M, and moving hydro energy between weeks.
- Build-limit policy: the case's inactive minimum capacity and CO2 inputs aren't carried over.

