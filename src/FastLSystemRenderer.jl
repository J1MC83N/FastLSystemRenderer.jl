
# Naming conventions:
#     seq - array of U8 representing an ASCII string
#     tally - count of number of a char in a seq; indexed to charset
#     pile - cumulation of a quantity through a seq; indexed to 1:length(seq)+1; always has first element zero
#     P - point matrix: 2xN Float
#     LI - line index offset array: index offset of point to draw a line to; similarly indexed to P
#         1 for no line; otherwise non-positive

module FastLSystemRenderer
export Lrender, Lcalculate, Lrender_bitmap, Lsystem


using DataStructures, StaticArrays, Octavian, Base.Threads, LinearAlgebra
using DataStructures: Forward
using LinearAlgebra: I

const Seq = Vector{UInt8}
const U8 = UInt8
const Float = Float32
const Charset = Vector{U8}

# Seq(str::String) = Array(codeunits(str))
charset(seq::Seq) = seq |> sort |> unique
charset(str::String) = charset(Seq(str))
charset(rules::Dict{UInt8,Seq}) = reduce(vcat,[vcat(Charset(seq),c) for (c,seq) in rules]) |> sort |> unique
charset(axiom::Union{U8,Seq},rules::Dict{UInt8,Seq}) = reduce(vcat,[vcat(Charset(seq),c) for (c,seq) in rules], init = axiom)  |> sort |> unique
charset(steplens::Dict{U8,Float}) = keys(steplens) |> collect |> sort |> unique

chartally(char::U8,charset::Charset=[char])::Vector{Int} = (out = zeros(Int,length(charset)); out[searchsortedfirst(charset,char)] = 1; out)
function chartally(seq::Seq, charset::Charset=charset(seq))::Vector{Int} #tally up the number of each char in a seq, return in sorted ordering 
    tallydict = Dict{U8,Int}()
    for c in seq
        tallydict[c] = 1 + get!(tallydict,c,0)
    end
    tally = zeros(Int,length(charset))
    for i in eachindex(tally)
        tally[i] += get!(tallydict,charset[i],0)
    end
    return tally
end

function tallymatrix(rules::Dict{U8,Seq},charset::Charset=charset(rules))::Matrix{Int} 
    #generate the tally iteration matrix; the tally on any iteration is M^niter*tally
    l = length(charset)
    M = Matrix{Int}(I,l,l)
    
    for i = 1:l
        c = charset[i]
        if haskey(rules,c)
            M[:,i] = chartally(rules[c],charset)
        end
    end
    return M
end

# given a Dict of stepping characters and their step lengths, cache their stepping vectors for all possible angles
function gensteps(steplens::Dict{U8,Float},anglediv::Int)
    steps = Dict{U8,LSC}()
    sins,coses = zeros(Float,anglediv), zeros(Float,anglediv)
    for i = 1:anglediv
        sins[i],coses[i] = round(sin(2π/anglediv*i),digits=15), round(cos(2π/anglediv*i),digits=15)
    end
    for (c,len) in steplens
        steps[c] = (len,sins.*len,coses.*len)
    end
    return steps
end

#Acronym for: step length, sines and cosines
const LSC = Tuple{Float,Vector{Float},Vector{Float}}
mutable struct Lsystem
    axiom::Seq
    rules::Dict{UInt8,Seq}
    anglediv::Int
    charset::Charset #TODO: use SortedSet
    tallymatrix::Matrix{Int}
    steps::Dict{U8,LSC}
    memo::Vector{Dict{UInt8,Seq}}
end
function Lsystem(axiom::Seq, rules::Dict{U8,Seq}, steplens::Dict{U8,Float}, anglediv::Int)::Lsystem #constructor
    cset = vcat(charset(axiom,rules),charset(steplens),U8('+'),U8('-')) |> sort |> unique # always insert the characters +-
    tallyM = tallymatrix(rules,cset)
    steps = gensteps(steplens,anglediv)
    memo_rule = Dict{U8,Seq}()
    for c in cset #completing the memoized replacement
        memo_rule[c] = haskey(rules,c) ? copy(rules[c]) : [c]
    end
    return Lsystem(axiom,rules,anglediv,cset,tallyM,steps,[memo_rule])
end

chartally(axiom::Union{U8,Seq}, rules::Dict{U8,Seq}, niter::Int, cset::Charset=charset(axiom,rules), tallymatrix::Matrix{Int}=tallymatrix(rules,cset))::Vector{Int} = tallymatrix^niter*chartally(axiom,cset)
chartally(axiom::Union{U8,Seq}, niter::Int, tallymatrix::Matrix{Int}, charset::Charset=charset(axiom,rules))::Vector{Int} = tallymatrix^niter*chartally(axiom,charset)
chartally(lsystem::Lsystem,niter::Int) = chartally(lsystem.axiom,lsystem.rules,niter,lsystem.charset,lsystem.tallymatrix)
chartally(lsystem::Lsystem,axiom::Union{U8,Seq},niter::Int) = chartally(axiom,lsystem.rules,niter,lsystem.charset,lsystem.tallymatrix)


###################################################
# Preperations & pre-processing: pile calculation #
###################################################

function debranch!(seq::Seq) #remove all branching subsequences
    pointer = 1
    depth = 0
    pointer_depth0 = 0
    while pointer <= length(seq)
        c = seq[pointer]
        if c == 0x5B
            depth == 0 && (pointer_depth0 = pointer)
            depth += 1
        elseif c == 0x5D
            depth -= 1
            if depth == 0
                deleteat!(seq,pointer_depth0:pointer)
                pointer = pointer_depth0-1
            end
        end
        pointer += 1
    end
    return seq
end
function debranch(seq::Seq)
    out = copy(seq)
    debranch!(out)
end
function debranch(rules::Dict{U8,Seq})
    out = deepcopy(rules)
    for (c,seq) in out
        debranch!(seq)
    end
    return out
end
hasbranch(charset::Charset) = 0x5B in charset && 0x5D in charset

norm_angle_anglediv(anglediv::Int) = ang -> mod(ang-1,anglediv)+1
function anglepile(axiom::Seq, niter::Int, rules::Dict{U8,Seq}, anglediv::Int, charset::Charset=charset(axiom,rules))::Vector{Int}
    # @assert 0x2B in charset && 0x2D in charset
    norm_angle = norm_angle_anglediv(anglediv)
    tallyM = tallymatrix(hasbranch(charset) ? debranch(rules) : rules, charset)
    index_plus,index_minus = searchsortedfirst(charset,0x2B),searchsortedfirst(charset,0x2D)
    l = length(axiom)
    pile = Vector{Int}(undef,l+1); pile[1] = anglediv #pile[i] is the angle to begin evaluating axiom[i]; pile[1]=0, pile[end=l+1] is the net angle
    stack = Stack{Int}()
    for i = 1:l
        c = axiom[i]
        if c == 0x5B # [
            push!(stack,pile[i])
            pile[i+1] = pile[i]
        elseif c == 0x5D # ]
            pile[i+1] = pop!(stack)
        else
            Δ = if haskey(rules,axiom[i])
                tally = chartally(c,rules,niter,charset,tallyM)
                tally[index_plus]-tally[index_minus]
            else
                (c == 0x2B) - (c == 0x2D)
            end
            pile[i+1] = norm_angle(Δ+pile[i])
        end
    end
    return pile
end
anglepile(lsystem::Lsystem,niter::Int) = anglepile(lsystem.axiom,niter,lsystem.rules,lsystem.anglediv,lsystem.charset)
anglepile(lsystem::Lsystem,axiom::Seq,niter::Int) = anglepile(axiom,niter,lsystem.rules,lsystem.anglediv,lsystem.charset)

const Coord = SArray{Tuple{2},Float,1,2} # 2D coordinate
@inline R(coord::Coord,θ::Float) = (@inbounds Coord((cos(θ)*coord[1]-sin(θ)*coord[2], sin(θ)*coord[1]+cos(θ)*coord[2])))
function displacementpile(axiom::Seq, niter::Int, rules::Dict{U8,Seq}, disps::Dict{U8,Coord}, anglediv::Int, charset::Charset=charset(axiom,rules))::Vector{Coord}
    l = length(axiom)
    baseangle = Float(2π)/anglediv
    angpile = anglepile(axiom,niter,rules,anglediv,charset)
    pile = [Coord(0,0)]
    stack = Stack{Coord}()
    for i = 1:l
        c = axiom[i]
        coord = if c == 0x5B # [
            push!(stack, copy(pile[i]))
            copy(pile[i])
        elseif c == 0x5D # ]
            pop!(stack)
        elseif haskey(disps,c) # has step length properties
            #TODO?: cache sin&cos here
            pile[i] .+ R(disps[c],baseangle*angpile[i])
        else
            copy(pile[i])
        end
        push!(pile,coord)
    end
    return pile
end
displacementpile(char::U8, niter::Int, rules::Dict{U8,Seq}, disps::Dict{U8,Coord}, anglediv::Int, charset::Charset=charset(axiom,rules)) = displacementpile(rules[char],niter-1,rules,disps,anglediv,charset)
function displacementpile(axiom::Seq, niter::Int, rules::Dict{U8,Seq}, steps::Dict{U8,LSC}, anglediv::Int, charset::Charset=charset(axiom,rules))
    disps = Dict{U8,Coord}() #dict containing all (potentially) coordinate displacing chars, and their respective displacement.
    for (c,seq) in rules; disps[c] = Coord(0,0) end #adding variable entires
    for (c,lsc) in steps; disps[c] = Coord(lsc[1],0) end #adding stepping chars, overriding the default 0-disp variables
    for i in 1:niter
        disps_next = Dict{U8,Coord}()
        for (c,lsc) in steps; disps_next[c] = Coord(lsc[1],0); end #adding stepping chars as default (for non-variables)
        for (c,seq) in rules
            disps_next[c] = displacementpile(c,i,rules,disps,anglediv,charset)[end]
        end
        disps = disps_next
    end
    return displacementpile(axiom,niter,rules,disps,anglediv,charset)
end
displacementpile(L::Lsystem,niter::Int) = displacementpile(L.axiom,niter,L.rules,L.steps,L.anglediv,L.charset)
displacementpile(L::Lsystem,axiom::Seq,niter::Int) = displacementpile(axiom,niter,L.rules,L.steps,L.anglediv,L.charset)

function pointpile(axiom::Seq, niter::Int, rules::Dict{U8,Seq}, steps::Dict{U8,LSC}, charset::Charset=charset(axiom,rules))
    tallyM = tallymatrix(rules,charset)
    l = length(axiom)
    pile = zeros(Int,l+1)
    for i in 1:l
        s = 0
        c = axiom[i]
        tally = chartally(c,niter,tallyM,charset)
        
        for j in eachindex(tally)
            s += haskey(steps,charset[j])*tally[j]
        end
        pile[i+1] = pile[i] + s
    end
    return pile
end
pointpile(L::Lsystem,niter::Int) = pointpile(L.axiom,niter,L.rules,L.steps,L.charset)
pointpile(L::Lsystem,axiom::Seq,niter::Int) = pointpile(axiom,niter,L.rules,L.steps,L.charset)

###############################################################################
# Derivation & execution #
###############################################################################

function Lderive_naive(axiom::Seq, rules::Dict{U8,Seq}, niter::Int)::Seq
    # @assert niter >= 0
    for i in 1:niter
        buff = IOBuffer(sizehint=sum(chartally(axiom,rules,1)))
        for c in axiom
            write(buff,get(rules,c,c))
        end
        axiom = buff.data
        #close(buff)
    end
    return axiom
end
Lderive_naive(lsystem::Lsystem,niter::Int) = Lderive_naive(lsystem.axiom,lsystem.rules,niter)

function Lmemoize!(L::Lsystem,niter::Int,MAX_MEMO_LENGTH::Int)
    memo,rules = L.memo,L.rules
    max_reached, breaking_iter = false, niter
    for i = length(memo)+1:niter
        D = Dict{U8,Seq}()
        for c in L.charset
            buff_size = sum(chartally(L,c,i))
            buff_size > MAX_MEMO_LENGTH && ((max_reached,breaking_iter) = (true,i))
            buff = IOBuffer(sizehint=buff_size)
            for c_rule in memo[1][c]
                write(buff,memo[i-1][c_rule])
            end
            D[c] = buff.data
        end
        push!(memo,D) #; close(buff)
        max_reached && break
    end
    return breaking_iter
end




const PLI = NamedTuple{(:P, :LI), Tuple{Matrix{Float32}, Vector{Int32}}}
const PLIi = NamedTuple{(:P, :LI, :lineind), Tuple{Matrix{Float32}, Vector{Int32}, Int32}}
PLIi(P::Matrix{Float32}, LI::Vector{Int32}, lineind::Int32) = PLIi((P,LI,lineind))
function initPLIi(L::Lsystem,axiom::Seq,iter::Int,pad::Bool)
    if !pad
        l = pointpile(L,axiom,iter)[end]
        P = Matrix{Float}(undef,2,l)
        LI = Vector{Int32}(undef,l)
        lineind = zero(Int32)
    else
        l = pointpile(L,axiom,iter)[end]+1
        P = Matrix{Float}(undef,2,l)
        P[1,1] = P[2,1] = zero(Float)
        LI = Vector{Int32}(undef,l)
        LI[1] = zero(Int32)
        lineind = one(Int32)
    end
    return PLIi(P,LI,lineind)
end

function execute_lines_0!(P_LIi::PLIi,pointer::Int,axiom::Seq,ang::Int,x::Float,y::Float,anglediv::Int,steps::Dict{U8,LSC})
    P,LI,lineind = P_LIi
    norm_angle = norm_angle_anglediv(anglediv) # bound ang to within [1,anglediv]
    stack_ang, stack_dist, stack_lineind = Stack{Int}(),Stack{Coord}(),Stack{Int32}()
    for c in axiom #loop through chars in the memoized sequence for this axiom char
        if c == 0x2B # +
            ang = norm_angle(ang+1)
        elseif c == 0x2D # -
            ang = norm_angle(ang-1)
        elseif c == 0x5B # [
            push!(stack_ang,ang)
            push!(stack_dist,Coord((x,y)))
            push!(stack_lineind,lineind)
        elseif c == 0x5D # ]
            ang = pop!(stack_ang)
            x,y = pop!(stack_dist)
            lineind = pop!(stack_lineind)
        elseif haskey(steps,c) # stepping char
            x += steps[c][3][ang] #  third of LSC -> cos cache
            y += steps[c][2][ang] # second of LSC -> sin cache
            P[1,pointer],P[2,pointer] = x,y
            LI[pointer] = lineind-pointer
            pointer += 1
            lineind = Int32(pointer-1)
        end
    end
    # @assert isempty(stack_ang)
    # @assert isempty(stack_dist)
    return PLIi(P,LI,lineind)
end
execute_lines_0!(P_LIi::PLIi,pointer::Int,axiom::Seq,anglediv::Int,steps::Dict{U8,LSC}) = execute_lines_0!(P_LIi,pointer,axiom,anglediv,zero(Float),zero(Float),anglediv,steps)
execute_lines_0(L::Lsystem,axiom::Seq;pad::Bool) = (P_LIi = initPLIi(L,axiom,0,pad); execute_lines_0!(P_LIi,pad+1,axiom,L.anglediv,L.steps))

function execute_lines_inmemo(L::Lsystem,axiom::Seq,niter::Int)
    memo, rules, anglediv, steps = L.memo, L.rules, L.anglediv, L.steps
    niter_memo = length(memo)
    @assert 0 < niter <= niter_memo #already memoized
    angpile, disppile, ptpile = anglepile(L,niter), displacementpile(L,niter), pointpile(L,niter)
    P = Matrix{Float}(undef,2,ptpile[end])
    LI = Vector{Int32}(undef,ptpile[end])
    #last/nearest point of the current scope
    stack_SLI = Stack{Int32}()
    lineind = Int32(0)
    for i_axiom in eachindex(axiom) #loop through chars in the axiom
        c_axiom = axiom[i_axiom]
        if c_axiom == 0x5B #[
            push!(stack_SLI,lineind)
        elseif c_axiom == 0x5D #]
            lineind = pop!(stack_SLI)
        elseif ptpile[i_axiom+1] - ptpile[i_axiom] > 0 # only when points are added
            pointer = ptpile[i_axiom]+1 #1 due to ptpile starting at zero and first point at pointer 1 not 0
            seq_memo = memo[niter][c_axiom]
            ang = angpile[i_axiom]; x,y = disppile[i_axiom]
            lineind = execute_lines_0!(PLIi(P,LI,lineind),pointer,seq_memo,ang,x,y,anglediv,steps).lineind
        end
    end
    return PLIi(P,LI,lineind)
end
execute_lines_inmemo(L::Lsystem,niter::Int) = execute_lines_inmemo(L,L.axiom,niter)

mycopyto!(P::Matrix{Float}, col::Int, Pc::Matrix{Float}) = (@boundscheck (@assert length(P) >= length(Pc)+2*(col-1)); unsafe_copyto!(P,col*2-1,Pc,1,length(Pc)))
mycopyto!(LI::Vector{Int32}, i::Int, LIc::Vector{Int32}) = (@boundscheck (@assert length(LI) >= length(LIc)+i-1); unsafe_copyto!(LI,i,LIc,1,length(LIc)))
function execute_lines_outmemo!(P_LIi::PLIi,pointer_offset::Int,seq::Seq,iter::Int,L::Lsystem,rotations::Vector{Matrix{Float}},exememo::Dict{U8,PLIi})
    rules, steps = L.rules, L.steps
    angpile, disppile, ptpile = anglepile(L,seq,iter), displacementpile(L,seq,iter), pointpile(L,seq,iter)
    P,LI,lineind = P_LIi
    @boundscheck (@assert size(P,2) >= pointer_offset+ptpile[end]-1)
    stack_SLI = Stack{Int32}()
    for i_seq in eachindex(seq)
        c = seq[i_seq]
        pointer = pointer_offset+ptpile[i_seq]+1
        if c == 0x5B #[
            push!(stack_SLI,lineind)
        elseif c == 0x5D #]
            lineind = pop!(stack_SLI)
        elseif haskey(rules,c) #move by memoized replacing char
            Pc = Octavian.matmul(rotations[angpile[i_seq]],exememo[c].P) .+ disppile[i_seq]
            mycopyto!(P,pointer,Pc)
            mycopyto!(LI,pointer,exememo[c].LI)
            #since exememo is generated from execute_line_0, which assumes the present of the local origin point at ind 0 which LI can point to, we need to reassign the LI elements pointing to ind0 to instead point to the current lineind, i.e. the last point of the scope
            for i_exememoLI in eachindex(exememo[c].LI)
                i_LI = i_exememoLI + pointer-1 #absolute index to LI; -1 due to 1-based array indexing
                if LI[i_LI] == -i_exememoLI # i points Δ=-i, i.e. ind0
                    LI[i_LI] = lineind-i_LI
                end
            end
            if exememo[c].lineind != 0 #the elsecase is when c doesn't alter the previous lineind, e.g. "F=", or "F=+[-F]-"
                lineind = Int32(pointer-1+exememo[c].lineind)
            end
        elseif haskey(steps,c) #move by non-replacing stepping char
            ang = angpile[i_seq]
            P[1,pointer] = steps[c][3][ang] + disppile[i_seq][1]
            P[2,pointer] = steps[c][2][ang] + disppile[i_seq][2]
            LI[pointer] = lineind-pointer
            lineind = Int32(pointer)
        end
    end
    return PLIi(P,LI,lineind)
end
execute_lines_outmemo(seq::Seq,iter::Int,L::Lsystem,rotations::Vector{Matrix{Float}},exememo::Dict{U8,PLIi}) = (P_LIi = initPLIi(L,seq,iter,false); execute_lines_outmemo!(P_LIi,0,seq,iter,L,rotations,exememo))
function execute_lines_outmemo(L::Lsystem,niter::Int)
    memo, axiom, anglediv, steps, rules = L.memo, L.axiom, L.anglediv, L.steps, L.rules
    niter_memo = length(memo)
    @assert niter_memo < niter
    rotations = [(ang->[cos(ang) -sin(ang); sin(ang) cos(ang)])(Float(angdiv*2π/anglediv)) for angdiv in 1:anglediv]
    #executing the furtherest memoized strings as a base start
    exememo = Dict{U8,PLIi}()
    for (c_rule,seq) in rules; exememo[c_rule] = execute_lines_0(L,memo[end][c_rule],pad=false) end
    #iteration start
    exememo_next = Dict{U8,PLIi}()
    for iter in niter_memo:niter-1
        for (c_rule,seq_rule) in rules
            exememo_next[c_rule] = execute_lines_outmemo(seq_rule,iter,L,rotations,exememo)
        end
        exememo = exememo_next
        exememo_next = Dict{U8,PLIi}()
    end
    P_LIi = initPLIi(L,axiom,niter,true) #padding a zero
    return execute_lines_outmemo!(P_LIi,1,axiom,niter,L,rotations,exememo)
end


###############################################################################
# Plotting & saving 
###############################################################################


function divlabor(NP::Int,NT::Int)
    n = div(NP,NT)
    out = [1+n*(i-1):n*i for i in 1:NT]
    out[end] = out[end].start:NP
    return out
end
function boundary(P::Matrix{Float})
    xmin = xmax = ymin = ymax = zero(Float)
    for i in 1:size(P,2)
        x,y = P[1,i], P[2,i]
        x < xmin && (xmin = x)
        x > xmax && (xmax = x)
        y < ymin && (ymin = y)
        y > ymax && (ymax = y)
    end
    return (xmin,xmax),(ymin,ymax)
end
expand(a::Float,b::Float,f::Float) = a-f*(b-a),b+f*(b-a)
expand(interval::NTuple{2,Float},f::Float) = interval[1]-f*(interval[2]-interval[1]),interval[2]+f*(interval[2]-interval[1])
topixel(min::Float,max::Float,n::Int,pad::Int) = function(pos::Float) round(Int,n*(pos-min)/(max-min))+pad end
topixel(interval::NTuple{2,Float},n::Int,pad::Int) = topixel(interval[1],interval[2],n,pad)
topixelrev(min::Float,max::Float,n::Int,pad::Int) = function(pos::Float) n-round(Int,n*(pos-min)/(max-min))+pad end
topixelrev(interval::NTuple{2,Float},n::Int,pad::Int) = topixelrev(interval[1],interval[2],n,pad)
function ORBuffs!(Buffs::Vector{BitMatrix})
    while length(Buffs) > 1
        rB = 1:2:length(Buffs)-1
        newBuffs = Vector{BitMatrix}(undef,length(rB)+isodd(length(Buffs)))
        isodd(length(Buffs)) && (newBuffs[end] = Buffs[end])
        @threads for iB in eachindex(rB)
            base = Buffs[rB[iB]]
            base .|= Buffs[rB[iB]+1]
            newBuffs[iB] = base
        end
        Buffs = newBuffs
        # print(" len", length(Buffs))
    end
    return Buffs[1]
end

#TODO: find a better to define plotlines and plotlinesT together, instead of using @eval
plotblock_LI = :(il0 = il+LI[il];
    if il0 == il-1; plotline!(Buff,x0,y0,x1,y1);
    elseif il0 != -1; plotline!(Buff,x2pix(P[1,il0]),y2pix(P[2,il0]),x1,y1) end)
plotblock_LIT = :(il0 = il+LI[il]; 
    if il0 == il-1; plotline!(Buff,y0,x0,y1,x1);
    elseif il0 != -1; plotline!(Buff,y2pix(P[2,il0]),x2pix(P[1,il0]),y1,x1) end)
plotblock_all = :(plotline!(Buff,x0,y0,x1,y1))
plotblock_allT = :(plotline!(Buff,y0,x0,y1,x1))
pix0_LI = :(LI[laborstart] == 1 ? (0,0) : (x2pix(P[1,laborstart+LI[laborstart]]), y2pix(P[2,laborstart+LI[laborstart]])))
pix0_all = :(laborstart == 1 ? (x2pix(0f0),y2pix(0f0)) : (x2pix(P[1,laborstart-1]),x2pix(P[2,laborstart-1])))
for (funcname,maindim,secddim,pix0,plotblock) in 
   ((:plotlines,     :h,:w,pix0_LI, plotblock_LI),
    (:plotlines_all, :h,:w,pix0_all,plotblock_all))
    # (:plotlinesT,    :w,:h,pix0_LI, plotblock_LIT)
    # (:plotlinesT_all,:w,:h,pix0_all,plotblock_allT)
    @eval function $funcname(P::Matrix{Float},LI::Vector{Int32},w::Int,h::Int,pad::Int=-1)
        # @assert size(P,2) == length(LI) #similarly indexed
        if pad == -1 #automatic padding
            pad = div(128 - h%64,2)
        end
        NT = nthreads()
        NL = size(P,2)
        Buffs = [zeros(Bool,h+2pad+isodd(h),w+2pad)|>BitArray for i in 1:NT]
        labordivs = divlabor(NL,NT)
        xintv, yintv = boundary(P)
        x2pix, y2pix = topixel(xintv,$maindim,pad), topixelrev(yintv,$secddim,pad)
        @threads for it in 1:NT
            # print(" ",it)
            id = threadid()
            isempty(labordivs[id]) && continue
            Buff = Buffs[id]
            labordiv = labordivs[id]
            laborstart = labordiv[1]
            x0,y0 = $pix0
            for il in labordiv
                x1,y1 = x2pix(P[1,il]), y2pix(P[2,il])
                $plotblock
                x0,y0 = x1,y1
            end
            # print("_", it)
        end
        Buff = ORBuffs!(Buffs)
        # plotline!(Buff,x2pix(0f0),y2pix(-1f1),x2pix(0f0),y2pix(1f1))
        # plotline!(Buff,x2pix(-1f1),y2pix(0f0),x2pix(1f1),y2pix(0f0))
        return Buff
    end
end
plotlines(P::Matrix{Float},w::Int,h::Int,pad::Int=-1) = plotlines_all(P,Int32[],w,h,pad)
plotlines(P_LIi::PLIi,w::Int,h::Int,pad::Int=-1) = plotlines(P_LIi.P,P_LIi.LI,w,h,pad)

birange(i::Int,j::Int) = range(i,j,step=(j>=i)*2-1)
function plotline!(buff::BitMatrix, x0::Int, y0::Int, x1::Int, y1::Int)
    # this is Bresenham's line algorithm that I can't recall for the life of me where I adapted from
    (x0==x1 && y0==y1) && (buff[x0,y0]=true; return buff)
    δx = x1 - x0
    δy = y1 - y0
    if abs(δx) >= abs(δy)
        yincr = δy >= 0 ? 1 : -1
        erdecr = δy >= 0 ? 1f0 : -1f0
        δe = δy / abs(δx)
        er = 0f0
        y = y0
        for x in birange(x0,x1)
            buff[x, y] = true
            er += δe
            if abs(er) > 0.5f0
                y  += yincr
                er -= erdecr
            end
        end
    else
        xincr = δx >= 0 ? 1 : -1
        erdecr = δx >= 0 ? 1f0 : -1f0
        δe = δx / abs(δy)
        er = 0f0
        x = x0
        for y in birange(y0,y1)
            buff[x, y] = true
            er += δe
            if abs(er) > 0.5f0
                x  += xincr
                er -= erdecr
            end
        end
    end
    return buff
end

function revU8(b::U8)::U8
    b = (b & 0xF0) >> 4 | (b & 0x0F) << 4
    b = (b & 0xCC) >> 2 | (b & 0x33) << 2
    b = (b & 0xAA) >> 1 | (b & 0x55) << 1
    return b
end

function savePBM(filename::String, B::BitMatrix)
    @assert filename[end-3:end] == ".pbm"
    size1,size2 = size(B)
    @assert size1 % 64 == 0 #this way chunks for columns of B stay clear cut
    chunks = B.chunks
    open(filename, "w") do s
        REVU8 = [revU8(i) for i = 0x00:0xFF]
        write(s, "P4\n")
        write(s, "# pbm file written by Julia\n")
        write(s, "$(size1) $(size2)\n")
        for chunk in chunks
            for i in 1:8
                write(s, REVU8[chunk%UInt8 + 1])
                chunk = chunk >> 8
            end
        end
    end
end

function importrules(str::String)
    rules = Dict{U8,Seq}()
    equate_str = "\\s*(?:=|[=|-]>)\\s*"
    for line in split(str,'\n')
        isempty(strip(line)) && continue
        line[1] == ';' && continue
        char, sub = strip(match(Regex("^([A-Z])$equate_str"),line).captures[1]), strip(match(Regex("$equate_str(.*)\\s*\$"),line).captures[1])
        @assert length(char) == 1
        rules[U8(char[1])] = Seq(sub)
    end
    return rules
end

function Lsystem(axiom::String, rulestr::String, anglediv::Integer, steplens::Dict{Char,<:Real}=Dict('F'=>1))
    _steplens = Dict{U8,Float}()
    for (c,len) in steplens
        _steplens[U8(c)] = Float(len)
    end
    Lsystem(Seq(axiom),importrules(rulestr),_steplens,Int(anglediv))
end

Lcalculate(L::Lsystem, niter::Int) = execute_lines_outmemo(L,niter)
Lrender_bitmap(res_calc::PLIi, width::Int, height::Int) = plotlines(res_calc, height, width)
Lsave(filename::String, bitmap::BitMatrix) = savePBM(filename, bitmap)
function Lrender(filename::String, L::Lsystem, niter::Int, width::Int, height::Int)
    timed = @timed Lsave(filename, Lrender_bitmap(Lcalculate(L,niter), width, height))
    nchar = sum(chartally(L,niter))
    return "L-system w/ $nchar characters successfully rendered in $(timed.time) seconds"
end

end