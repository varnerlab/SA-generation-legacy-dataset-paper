import Base.indexin

function indexin(model::BSTModel, species_symbol::String; 
    key::String="total_species_list")::Union{Nothing,Int}

    # get the total_species_list from the model
    total_species_list = model.total_species_list;

    # get the total species list -
    if (in(species_symbol,total_species_list) == false)
        return nothing
    end

    # check -
    return findfirst(x->x==species_symbol,total_species_list)
end

function findfiles(path::String, pattern::String)::Array{String,1}

    # initialize -
    file_list = Array{String,1}()

    # get the list of files -
    for file in readdir(path)
        if occursin(pattern, file)
            push!(file_list, file)
        end
    end

    # return -
    return file_list
end
