function execute_0!(P::Matrix{Float},pointer::Int,axiom::Seq,ang::Int,x::Float,y::Float,anglediv::Int,steps::Dict{U8,LSC})
    norm_angle = norm_angle_anglediv(anglediv) # bound ang to within [1,anglediv]
    stack_ang, stack_dist = Stack{Int}(),Stack{Coord}()
    for c in axiom #loop through chars in the memoized sequence for this axiom char
        if c == 0x2B # +
            ang = norm_angle(ang+1)
        elseif c == 0x2D # -
            ang = norm_angle(ang-1)
        elseif c == 0x5B # [
            push!(stack_ang,ang)
            push!(stack_dist,Coord((x,y)))
        elseif c == 0x5D # ]
            ang = pop!(stack_ang)
            x,y = pop!(stack_dist)
        elseif haskey(steps,c) # stepping char
            x += steps[c][3][ang] #  third of LSC -> cos cache
            y += steps[c][2][ang] # second of LSC -> sin cache
            P[1,pointer],P[2,pointer] = x,y
            pointer += 1
        end
    end
    # @assert isempty(stack_ang)
    # @assert isempty(stack_dist)
    return P
end
initP(L::Lsystem,axiom::Seq) = Matrix{Float}(undef,2,pointpile(L,axiom,0)[end])
execute_0!(P::Matrix{Float},axiom::Seq,anglediv::Int,steps::Dict{U8,LSC}) = execute_0!(P,1,axiom,anglediv,zero(Float),zero(Float),anglediv,steps)
execute_0(L::Lsystem,axiom::Seq) = (P = initP(L,axiom); execute_0!(P,axiom,L.anglediv,L.steps); P)
execute_0(L::Lsystem) = execute_0(L,L.axiom)

function execute_inmemo(L::Lsystem,axiom::Seq,niter::Int)
    memo, rules, anglediv, steps = L.memo, L.rules, L.anglediv, L.steps
    niter_memo = length(memo)
    @assert niter <= niter_memo #already memoized
    @assert niter > 0
    angpile, disppile, ptpile = anglepile(L,niter), displacementpile(L,niter), pointpile(L,niter)
    P = Matrix{Float}(undef,2,ptpile[end])
    for i_axiom in eachindex(axiom) #loop through chars in the axiom
        c_axiom = axiom[i_axiom]
        if ptpile[i_axiom+1] - ptpile[i_axiom] > 0 # only when points are added
            pointer = ptpile[i_axiom]+1 #1 due to ptpile starting at zero and first point at pointer 1 not 0
            seq_memo = memo[niter][c_axiom]
            ang = angpile[i_axiom]; x,y = disppile[i_axiom]
            execute_0!(P,pointer,seq_memo,ang,x,y,anglediv,steps)
            # @assert pointer == ptpile[i_axiom+1]+1
            # @assert ang == angpile[i_axiom+1]
            # @assert all(isapprox.((x,y), disppile[i_axiom+1].data, atol = 1E-4))
        end
    end
    return P
end
execute_inmemo(L::Lsystem,niter::Int) = execute_inmemo(L,L.axiom,niter)

function execute_outmemo(seq::Seq,iter::Int,L::Lsystem,rotations::Vector{Matrix{Float}},exememo::Dict{U8,Matrix{Float}})
    rules, steps = L.rules, L.steps
    angpile, disppile, ptpile = anglepile(L,seq,iter), displacementpile(L,seq,iter), pointpile(L,seq,iter)
    P = Matrix{Float}(undef,2,ptpile[end])
    for i_seq in eachindex(seq)
        c = seq[i_seq]
        if haskey(rules,c) #move by memoized replacing char
            Pc = Octavian.matmul(rotations[angpile[i_seq]],exememo[c]) .+ disppile[i_seq]
            mycopyto!(P,ptpile[i_seq]+1,Pc)
        elseif haskey(steps,c) #move by non-replacing stepping char
            ang = angpile[i_seq]
            P[1,ptpile[i_seq]+1] = steps[c][3][ang] + disppile[i_seq][1]
            P[2,ptpile[i_seq]+1] = steps[c][2][ang] + disppile[i_seq][2]
        end
    end
    return P
end
mycopyto!(P::Matrix{Float}, col::Int, Pc::Matrix{Float}) = (@boundscheck (@assert length(P) >= length(Pc)+2*(col-1)); unsafe_copyto!(P,col*2-1,Pc,1,length(Pc)))
function execute_outmemo(L::Lsystem,niter::Int)
    memo, axiom, anglediv, steps, rules = L.memo, L.axiom, L.anglediv, L.steps, L.rules
    niter_memo = length(memo)
    @assert niter_memo < niter
    rotations = [(ang->[cos(ang) -sin(ang); sin(ang) cos(ang)])(Float(angdiv*2π/anglediv)) for angdiv in 1:anglediv]
    #executing the furtherest memoized strings as a base start
    exememo = Dict{U8,Matrix{Float}}()
    for (c_rule,seq) in rules; exememo[c_rule] = execute_0(L,memo[end][c_rule]) end
    #iteration start
    exememo_next = Dict{U8,Matrix{Float}}()
    for iter in niter_memo:niter-1
        for (c_rule,seq_rule) in rules
            exememo_next[c_rule] = execute_outmemo(seq_rule,iter,L,rotations,exememo)
        end
        exememo = exememo_next
        exememo_next = Dict{U8,Matrix{Float}}()
    end
    P = execute_outmemo(axiom,niter,L,rotations,exememo)
    return P
end



for (funcname,maindim,img_maindim,buffinit,ipix,jpix) in ((:plotpoints, :h,:imgh,:(zeros(Bool,w+2pad,imgh)),:xpix,:ypix),
                                                          (:plotpointsT,:w,:imgw,:(zeros(Bool,imgw,h+2pad)),:ypix,:xpix))
    @eval function $funcname(P::Matrix{Float},w::Int,h::Int,pad::Int=-1)
        imgw,imgh = w,h
        if pad == -1 #automatic padding
            pad = div(64 - $maindim%64,2)
            $img_maindim = $maindim + 2pad + isodd($maindim)
        end
        NT = nthreads()
        NP = size(P,2)
        Buffs = [$buffinit|>BitArray for i in 1:NT]
        labordivs = divlabor(NP,NT)
        xintv, yintv = boundary(P)
        x2pix, y2pix = topixel(xintv,w,pad), topixel(yintv,h,pad)
        
        @threads for it in 1:NT
            print(" ",it)
            id = threadid()
            Buff = Buffs[id]
            xpix_prev,ypix_prev = 0,0
            for ip in labordivs[id]
                x, y = P[1,ip], P[2,ip]
                xpix,ypix = x2pix(x),y2pix(y)
                if xpix!=xpix_prev || ypix!=ypix_prev
                    Buff[$ipix,$jpix] = true
                end
                xpix_prev,ypix_prev = xpix,ypix
            end
            print("_", it)
        end
        Buff = ORBuffs!(Buffs)
        Buff[1][x2pix(zero(Float)),y2pix(zero(Float))] = true
        return Buff
    end
end
