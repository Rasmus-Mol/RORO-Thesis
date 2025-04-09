# Script to do something with Olivers variance of weight data
# Uses data from variance of weight from Oliver

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

using DataFrames, XLSX, CSV, Plots, Distributions, CategoricalArrays, StatsBase

# Function to plot distribution and compare with normal distribution - Weight
# df_secu and df_trailer should have the same number of quantiles as input
function scenario_distribution_with_normal(cargoC::CargoCollection, df_secu, df_trailer, n_quantiles::Int)
    if maximum(df_secu.QuantileNumber) != n_quantiles || maximum(df_trailer.QuantileNumber) != n_quantiles
        error("Number of quantiles in df does not match n_quantiles")
        return
    end
    # Get the weights of the cargo
    # Seperate scenario weight according to type
    scenario_secu_weight = collect(filter(items -> items.cargo_type_id == 4, cargoC.items).weight)
    scenario_trailer_weight = collect(filter(items -> items.cargo_type_id == 1, cargoC.items).weight)
    # Sort scenario weight into quantiles
    quantile_secu = [quantile(df_secu.CorrectedWeight, i / n_quantiles) for i in 1:n_quantiles]./1000
    quantile_trailer = [quantile(df_trailer.CorrectedWeight, i / n_quantiles) for i in 1:n_quantiles]./1000

    if n_quantiles >= 3
        secu_weight = [sort(filter(w -> w <= quantile_secu[i] && w > quantile_secu[i-1] , scenario_secu_weight)) for i in 2:(n_quantiles-1)]
        push!(secu_weight, sort(filter(w -> w > quantile_secu[n_quantiles-1], scenario_secu_weight)))
        pushfirst!(secu_weight, sort(filter(w -> w <= quantile_secu[1], scenario_secu_weight)))
        trailer_weight = [sort(filter(w -> w <= quantile_trailer[i] && w>quantile_trailer[i-1], scenario_trailer_weight)) for i in 2:(n_quantiles-1)]
        push!(trailer_weight, sort(filter(w -> w > quantile_trailer[n_quantiles-1], scenario_trailer_weight)))
        pushfirst!(trailer_weight, sort(filter(w -> w <= quantile_trailer[1], scenario_trailer_weight)))
        # Sort scenario weight into quantiles
        historic_secu_weight = [collect(filter(x -> x.QuantileNumber == i, df_secu).CorrectedWeight) for i in 1:n_quantiles]./1000
        historic_trailer_weight = [collect(filter(x -> x.QuantileNumber == i, df_trailer).CorrectedWeight) for i in 1:n_quantiles]./1000
    elseif n_quantiles == 2 # 2 quantiles
        secu_weight = [sort(filter(w -> w <= quantile_secu[1], scenario_secu_weight)),
                       sort(filter(w -> w > quantile_secu[1], scenario_secu_weight))]
        trailer_weight = [sort(filter(w -> w <= quantile_trailer[1], scenario_trailer_weight)),
                          sort(filter(w -> w > quantile_trailer[1], scenario_trailer_weight))]
        # Sort scenario weight into quantiles
        historic_secu_weight = [collect(filter(x -> x.QuantileNumber == i, df_secu).CorrectedWeight) for i in 1:n_quantiles]./1000
        historic_trailer_weight = [collect(filter(x -> x.QuantileNumber == i, df_trailer).CorrectedWeight) for i in 1:n_quantiles]./1000
    else # 1 quantile, ie. all data together
        secu_weight = [sort(scenario_secu_weight)]
        trailer_weight = [sort(scenario_trailer_weight)]
        # Sort scenario weight into quantiles
        historic_secu_weight = [collect(filter(x -> x.QuantileNumber == i, df_secu).CorrectedWeight) for i in 1:n_quantiles]./1000
        historic_trailer_weight = [collect(filter(x -> x.QuantileNumber == i, df_trailer).CorrectedWeight) for i in 1:n_quantiles]./1000
    end
    p = []
    for i in 1:n_quantiles
        if length(secu_weight[i])>0
            pl = histogram(secu_weight[i], normalize=true, label="Secu weight",
            ylabel = "Frequency", xlabel = "Weight", title = "Secu weight quantile $(i)")
            fit_secu = fit(Normal, historic_secu_weight[i])
            x = range(minimum(historic_secu_weight[i]), stop=maximum(historic_secu_weight[i]), length=100)
            plot!(x, pdf.(fit_secu, x), label="Fitted Normal", lw=2)
            push!(p,pl)
        end
    end
    xplots = Int(ceil(n_quantiles/2))
    yplots = Int(ceil(n_quantiles/xplots))
    plo = [plot(p..., layout=(xplots,yplots))]
    p = []
    for i in 1:n_quantiles
        if length(trailer_weight[i])>0
            pl = histogram(trailer_weight[i], normalize=true, label="Secu weight",
            ylabel = "Frequency", xlabel = "Weight", title = "Trailer weight quantile $(i)")
            fit_trailer = fit(Normal, historic_trailer_weight[i])
            x = range(minimum(trailer_weight[i]), stop=maximum(trailer_weight[i]), length=100)
            plot!(x, pdf.(fit_trailer, x), label="Fitted Normal", lw=2)
            push!(p,pl)
        end
    end
    xplots = Int(ceil(n_quantiles/2))
    yplots = Int(ceil(n_quantiles/xplots))
    push!(plo,plot(p..., layout=(xplots,yplots)))
    return plo
end


# Function to plot distribution and compare with normal distribution - Variance
# df_secu and df_trailer should have the same number of quantiles as input
function scenario_distribution_with_normal(cargoC::CargoCollection, df_secu, df_trailer, n_quantiles::Int)
    if maximum(df_secu.QuantileNumber) != n_quantiles || maximum(df_trailer.QuantileNumber) != n_quantiles
        error("Number of quantiles in df does not match n_quantiles")
        return
    end
    # Get the weights of the cargo
    # Seperate scenario weight according to type
    scenario_secu_weight = collect(filter(items -> items.cargo_type_id == 4, cargoC.items).weight)
    scenario_trailer_weight = collect(filter(items -> items.cargo_type_id == 1, cargoC.items).weight)
    # Sort scenario weight into quantiles
    quantile_secu = [quantile(df_secu.CorrectedWeight, i / n_quantiles) for i in 1:n_quantiles]./1000
    quantile_trailer = [quantile(df_trailer.CorrectedWeight, i / n_quantiles) for i in 1:n_quantiles]./1000

    if n_quantiles >= 3
        secu_weight = [sort(filter(w -> w <= quantile_secu[i] && w > quantile_secu[i-1] , scenario_secu_weight)) for i in 2:(n_quantiles-1)]
        push!(secu_weight, sort(filter(w -> w > quantile_secu[n_quantiles-1], scenario_secu_weight)))
        pushfirst!(secu_weight, sort(filter(w -> w <= quantile_secu[1], scenario_secu_weight)))
        trailer_weight = [sort(filter(w -> w <= quantile_trailer[i] && w>quantile_trailer[i-1], scenario_trailer_weight)) for i in 2:(n_quantiles-1)]
        push!(trailer_weight, sort(filter(w -> w > quantile_trailer[n_quantiles-1], scenario_trailer_weight)))
        pushfirst!(trailer_weight, sort(filter(w -> w <= quantile_trailer[1], scenario_trailer_weight)))
        # Sort scenario weight into quantiles
        historic_secu_weight = [collect(filter(x -> x.QuantileNumber == i, df_secu).CorrectedWeight) for i in 1:n_quantiles]./1000
        historic_trailer_weight = [collect(filter(x -> x.QuantileNumber == i, df_trailer).CorrectedWeight) for i in 1:n_quantiles]./1000
    elseif n_quantiles == 2 # 2 quantiles
        secu_weight = [sort(filter(w -> w <= quantile_secu[1], scenario_secu_weight)),
                       sort(filter(w -> w > quantile_secu[1], scenario_secu_weight))]
        trailer_weight = [sort(filter(w -> w <= quantile_trailer[1], scenario_trailer_weight)),
                          sort(filter(w -> w > quantile_trailer[1], scenario_trailer_weight))]
        # Sort scenario weight into quantiles
        historic_secu_weight = [collect(filter(x -> x.QuantileNumber == i, df_secu).CorrectedWeight) for i in 1:n_quantiles]./1000
        historic_trailer_weight = [collect(filter(x -> x.QuantileNumber == i, df_trailer).CorrectedWeight) for i in 1:n_quantiles]./1000
    else # 1 quantile, ie. all data together
        secu_weight = [sort(scenario_secu_weight)]
        trailer_weight = [sort(scenario_trailer_weight)]
        # Sort scenario weight into quantiles
        historic_secu_weight = [collect(filter(x -> x.QuantileNumber == i, df_secu).CorrectedWeight) for i in 1:n_quantiles]./1000
        historic_trailer_weight = [collect(filter(x -> x.QuantileNumber == i, df_trailer).CorrectedWeight) for i in 1:n_quantiles]./1000
    end
    p = []
    for i in 1:n_quantiles
        if length(secu_weight[i])>0
            pl = histogram(secu_weight[i], normalize=true, label="Secu weight",
            ylabel = "Frequency", xlabel = "Weight", title = "Secu weight quantile $(i)")
            fit_secu = fit(Normal, historic_secu_weight[i])
            x = range(minimum(historic_secu_weight[i]), stop=maximum(historic_secu_weight[i]), length=100)
            plot!(x, pdf.(fit_secu, x), label="Fitted Normal", lw=2)
            push!(p,pl)
        end
    end
    xplots = Int(ceil(n_quantiles/2))
    yplots = Int(ceil(n_quantiles/xplots))
    plo = [plot(p..., layout=(xplots,yplots))]
    p = []
    for i in 1:n_quantiles
        if length(trailer_weight[i])>0
            pl = histogram(trailer_weight[i], normalize=true, label="Secu weight",
            ylabel = "Frequency", xlabel = "Weight", title = "Trailer weight quantile $(i)")
            fit_trailer = fit(Normal, historic_trailer_weight[i])
            x = range(minimum(trailer_weight[i]), stop=maximum(trailer_weight[i]), length=100)
            plot!(x, pdf.(fit_trailer, x), label="Fitted Normal", lw=2)
            push!(p,pl)
        end
    end
    xplots = Int(ceil(n_quantiles/2))
    yplots = Int(ceil(n_quantiles/xplots))
    push!(plo,plot(p..., layout=(xplots,yplots)))
    return plo
end