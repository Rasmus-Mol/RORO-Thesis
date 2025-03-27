# Script to do something with Olivers variance of weight data


# Uses data from variance of weight from Oliver

using DataFrames, XLSX, CSV
file_path = "/home/rasmusmolsted/Documents/DTU/Onedrive-filer/Speciale/RORO-Thesis/OliverModelStochastic/data/CargoWeights.csv"
data = CSV.read(file_path, DataFrame; delim=",")
names(data)
data[!,"ID"] # <- access column by name
# Removing rows with missing values in specific column
df = dropmissing(data, :"Count.Booked Weight")
dropmissing!(df, :"Corrected Weight")
dropmissing!(df, :"Variance")