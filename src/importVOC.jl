# Helper functions to render *Vision of Chaos* LSystem rule files
function read2Lstrs(lines::Vector{String})
    Lstrs = Vector{String}[]
    params = Vector{String}(undef,4); params[4] = "" #name, angle, axiom, rules
    isinbracket = false
    for line in lines
        if !isinbracket && occursin('{',line) #beginning line of a Lstrs
            isinbracket = true
            params[1] = match(r"\s*(\S+)\s*\{",line).captures[1]
        elseif isinbracket && occursin(r"\s*angle\s\d+"i,line) #Angle
            params[2] = match(r"\s*angle\s(\d+)"i,line).captures[1]
        elseif isinbracket && occursin(r"\s*axiom\s\S+"i,line) #Axiom
            params[3] = match(r"\s*axiom\s(\S+)"i,line).captures[1]
        elseif isinbracket && occursin(r"\s*=\S*",line) #rule
            params[4] *= strip(line) * "\n"
        elseif isinbracket && occursin('}',line)
            isinbracket = false
            push!(Lstrs,params)
            params = Vector{String}(undef,4); params[4] = ""
        else
            # println(line)
        end
    end
    Lstrs
end
isbad_rulestr(rulestr::String) = occursin(r"m|M|G|g|<|>|@|\/|[\\]",rulestr) || !occursin(r"f|F|d|D",rulestr)
function Lstrs_filterbad!(Lstrs)
    for i in length(Lstrs):-1:1
        isbad_rulestr(Lstrs[i][4]) && deleteat!(Lstrs,i)
    end
    Lstrs
end
function Lstrs_process!(Lstrs)
    for params in Lstrs
        params[3] = uppercase(params[3])
        params[4] = uppercase(replace(params[4],r"c\d\d" => ""))
        if occursin('|',params[4])
            anglediv = parse(Int,params[2])
            if iseven(anglediv)
                params[4] = replace(params[4],"|" => "+"^div(anglediv,2))
            else
                println("Odd turnaround at ", params[1])
            end
        end
    end
    Lstrs
end

function importVoC(DVoC::Dict{String,Vector{Vector{String}}})
    Dout = Dict{String,Dict{String,Lsystem}}()
    steps = Dict(U8('F') => one(Float), U8('D') => one(Float))
    for (modulename,Lstrs) in DVoC
        isempty(Lstrs) && continue
        D = Dict{String,Lsystem}()
        for params in Lstrs
            name,angdivstr,axiom,rulestr = params
            L = Lsystem(Seq(axiom),importrules(rulestr),steps,parse(Int,angdivstr))
            if pointpile(L,10)[end] != 0
                D[name] = L
            end
        end
        Dout[modulename] = D
    end
    return Dout
end

# rendering workflow
# cd("Visions of Chaos/Examples/Data")
# DVOC = Dict{String,Vector{Vector{String}}}()
# for filename in filter(name->occursin(r"\.l$",name),readdir())
#     lines = open(filename) do f
#         readlines(f)
#     end
#     Lstrs = read2Lstrs(lines)
#     Lstrs_filterbad!(Lstrs)
#     Lstrs_process!(Lstrs)
#     global DVOC[filename[1:end-2]] = Lstrs
# end
# LVOC = importVoC(DVOC)
# maxiter(L::Lsystem,point_count_threshold::Int) = (for i = 1:30; pointpile(L,i)[end] > point_count_threshold && return i end; return 1)
# cd("renders")
# for modulename in union(setdiff(keys(LVOC),readdir()),["MCWORTER1"])
#     println(modulename,"========================================")
#     modulename in readdir() || mkdir(modulename)
#     cd(modulename)
#     for (name,L) in LVOC[modulename]
#         print(name)
#         PLIout = execute_lines_outmemo(L,maxiter(L,100000))
#         begin
#             BA = plotlines(PLIout,2000,2000)
#             filename = string(name,"-",hash(BA),".pbm")
#             if filename in readdir()
#                 @warn "Already exists"
#             else
#                 savePBM(filename,BA)
#             end
#         end
#         println()
#     end
#     cd("..")
# end
