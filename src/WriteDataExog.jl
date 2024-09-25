module WriteDataExog
import DelimitedFiles: writedlm
using Printf
using CSV
using DataFrames

function writedata_exog(TE::NamedTuple, exogindex::Int, R::String)
    # Constructing a label string
    labeller = "exog" * @sprintf("%02d", exogindex)

    # Generating year indices for different categories
    yearindex_share = Vector{Int64}(undef, 30)
    yearindex_share .= collect(1:30) .+ 2020

    # Creating the share path matrix
    sharepath = hcat(
        yearindex_share,
        100 .* TE.renewshare_path_region[:, 1:30]',
        100 .* TE.renewshareUS[1:30],
        100 .* TE.renewshare_path_world[:, 1:30]'
    )

    writedlm("$R/Renewable_share$(labeller).csv", sharepath, ",")

    
end

end