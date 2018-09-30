using CSV
using DataFrames
using SeedMix

tablefile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Plant_data_updated_8_3_18.csv"
phylofile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Distance_ALL_7_31_18.csv"

sd = SpeciesData(tablefile, phylofile)
