using Random, DataFrames, DataFramesMeta, CSV

# Create jitter data of horizontal axis in Figure 4 and Figure 7
function jitter_csv(;
                    # Nalgos = 11, Ncombinations = 30,
                    Nalgos = 9, Ncombinations = 42,
                    savepath = homedir(),
                    # filestr="gauss"
                    filestr="multi"
                    )
   df = DataFrame()
   from = 0.
   for i in 1:Nalgos
       df[Symbol("x"*string(i))] = rand(Ncombinations) .+ from
       from += 1.5
    end
    CSV.write(joinpath(savepath, "jitter_"* filestr *".csv"),
        df, delim = " ")
end
