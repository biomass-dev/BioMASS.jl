module V

const NAMES = [
    
]

for (idx,name) in enumerate(NAMES)
    eval(Meta.parse("const $name = $idx"))
end

const NUM = length(NAMES)

end  # module