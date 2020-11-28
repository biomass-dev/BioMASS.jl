const observables = [
    "Phosphorylated_MEKc"
    "Phosphorylated_ERKc"
    "Phosphorylated_RSKw"
    "Phosphorylated_CREBw"
    "dusp_mRNA"
    "cfos_mRNA"
    "cFos_Protein"
    "Phosphorylated_cFos"
]

function observables_index(observable_name::String)::Int
    if !(observable_name in observables)
        error("$observable_name is not defined in observables.")
    end
    return findfirst(isequal(observable_name),observables)
end