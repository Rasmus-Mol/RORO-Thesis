# Script to do something with Olivers variance of weight data
# Uses data from variance of weight from Oliver
#using DataFrames, XLSX, CSV, Plots, Distributions, CategoricalArrays, StatsBase

function load_Weight_Variance_data(file_path::String)
    # load data
    data = CSV.read(file_path, DataFrame; delim=",")
    # Removing rows with missing values in specific column
    df = dropmissing(data, :"Count.Booked Weight")
    dropmissing!(df, :"Corrected Weight")
    dropmissing!(df, :"Variance")
    # Rename columns
    rename!(df, :"Count.StartTime" => :CountStartTime)
    rename!(df, :"Count.UnitTypeID" => :CountUnitTypeID)
    rename!(df, :"Count.Booked Weight" => :CountBookedWeight)
    rename!(df, :"Count.VGM Weight" => :CountVGMWeight)
    rename!(df, :"Corrected Weight" => :CorrectedWeight)
    return df
end

function weight_difference_info(df,print::Bool=true)
    df_secu = df[df.CountUnitTypeID .== "SWAP45", :]
    df_trailer = df[df.CountUnitTypeID .== "TRA", :]

    if print
        names(df)
        unit_types = unique(df.CountUnitTypeID)
        println(unit_types)
        plot(df_secu.Variance, ylabel = "Weight variance (kg)", xlabel = "Count", title = "Variance of weight for Secu-box", label = "Secu-box")
        savefig("Plots/Data/secu_weight_diff.png")
        plot(df_trailer.Variance, ylabel = "Weight variance (kg)", xlabel = "Count", title = "Variance of weight for trailer", label = "trailer")
        savefig("Plots/Data/trailer_weight_diff.png")
        mean_secu = mean(df_secu.Variance)
        mean_trailer = mean(df_trailer.Variance)
        println("Mean variance for Secu-box: $mean_secu")
        println("Mean variance for trailer: $mean_trailer")
        count_secu = count(x -> x>0 , df_secu.Variance)
        count_trailer = count(x -> x>0 , df_trailer.Variance)
        println("Number of Trailers: $(length(df_trailer.Variance))")
        println("Number of times Trailers weighs less than informed: $count_trailer")
        println("Number of times Trailers weighs more than informed: $(count(x -> x<0 , df_trailer.Variance))")
        println("Number of times Trailers weighs were accurately informed: $(count(x -> x==0 , df_trailer.Variance))")
        println("Number of Secu-boxes: $(length(df_secu.Variance))")
        println("Number of times Secu-boxes weighs less than informed: $count_secu")
        println("Number of times Secu-boxes weighs more than informed: $(count(x -> x<0 , df_secu.Variance))")
        println("Number of times Secu-boxes weighs were accurately informed: $(count(x -> x==0 , df_secu.Variance))")    
    end
    return df_secu, df_trailer
end

# Seperates data into equally spaced bins
function seperate_data_into_bins(df,bins,remove_outlier::Bool=true)
    df_new = copy(df)
    if remove_outlier
        df_new = filter(row -> abs(row.CountBookedWeight - mean(df.CountBookedWeight)) / std(df.CountBookedWeight) ≤ 3, df)
    end
    df_sort = sort(df_new, :CountBookedWeight)
    min_val = df_sort[1,:CountBookedWeight]
    max_val = df_sort[end,:CountBookedWeight]
    # Calculating bin-sizes
    bin_size = (max_val - min_val)/bins
    bin_arr = []
    bin_start = min_val
    bin_end = min_val
    for i in 1:bins
        bin_end = bin_end + bin_size
        push!(bin_arr, [bin_start, bin_end]) 
        bin_start = bin_end
    end
    # adding bin numbers to each data point
    bin_number = []
    for i in 1:length(df_sort.CountBookedWeight)
        for j in 1:bins
            if df_sort.CountBookedWeight[i] <= bin_arr[j][2]
                push!(bin_number, j)
                break
            end
        end
    end
    println(length(bin_number))
    println(length(df_sort.CountBookedWeight))
    df_sort.BinNumber = bin_number
    return df_sort, bin_arr
end
# Seperate data into quantiles based on booked weight
function seperate_data_into_quantiles(df,quantiles,remove_outlier::Bool=true)
    df_new = copy(df)
    if remove_outlier
        df_new = filter(row -> abs(row.CountBookedWeight - mean(df.CountBookedWeight)) / std(df.CountBookedWeight) ≤ 3, df)
    end
    df_sort = sort(df_new, :CountBookedWeight)
    # Calculate quantile values (e.g., 0.25, 0.5, 0.75 for quartiles)
    quantile_arr = [quantile(df_sort.CountBookedWeight, i / quantiles) for i in 1:quantiles]
    
    # adding quantile numbers to each data point
    quantile_number = []
    for i in 1:length(df_sort.CountBookedWeight)
        for j in 1:quantiles
            if df_sort.CountBookedWeight[i] <= quantile_arr[j]
                push!(quantile_number, j)
                break
            end
        end
    end
    df_sort.QuantileNumber = quantile_number
    return df_sort, quantile_arr
end
