# FastLSystemRenderer

[![Build Status](https://github.com/J1MC83N/FastLSystemRenderer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/J1MC83N/FastLSystemRenderer.jl/actions/workflows/CI.yml?query=branch%3Amain)

Multi-threaded Lindenmayer Systems (L-systems) rendering in Julia. 


## Examples
```julia
using FastLSystemRenderer

# simple Koch snowflake for 10 iterations
rulestr = "F -> F-F++F-F"
axiom = "F++F++F"
turn_angle = 6 # 360∘/6 = 60∘
niter = 10
koch_snowflake = Lsystem(axiom,rulestr,turn_angle)
Lrender("koch_snowflake.pbm",koch_snowflake,niter,2000,2300)
```

![koch_snowflake](https://user-images.githubusercontent.com/40559557/215311192-31731801-d9d6-4897-bdec-085a6490ac89.png)


```julia
#Plant02: Adrian Mariano, from the Algorithmic Beauty of Plants
plant02 = begin
    rulestr = "F -> F[+F]F[-F][F]"
    axiom = "F"
    turn_angle = 18 #360∘/18 = 20∘
    Lsystem(axiom,rulestr,turn_angle)
end
niter = 9
FastLSystemRenderer.Lrender_time("plant02-n$niter.pbm",plant02,niter,8000,5000)
```
![plant02-n9](https://user-images.githubusercontent.com/40559557/215312164-ea794a30-0711-419a-bf99-a5ba5cd61c3a.png)

```julia
# two symmetric binomial trees
binotree = begin
    rulestr = "X -> F[+X][FX]
    Y -> F[-Y][FY]
    F -> H
    H -> FF"
    stepping_chars = Dict('F' => 1, 'H' => √2)
    turn_angle = 13
    axiom = "[X]Y"
    Lsystem(axiom,rulestr,turn_angle,stepping_chars)
end
niter = 18
Lrender("binotree-n$niter.pbm",binotree,niter,6000,5000)
```

![binotree-n18](https://user-images.githubusercontent.com/40559557/215314040-6c3eef21-356c-4a42-a599-b3ce953a57f0.png)



## Features & limitations

Currently supported rule characters: 
* `+`: turn anticlockwise by turning angle
* `-`: turn clockwise by turning angle
* `[`: push current drawing state onto stack
* `]`: pop current drawing state from the stack
* Customizable set of stepping (moving forward while drawing a line) characters with arbitrary stepping distance (see third example – binomial tree – for a use case. Default stepping character is `F` with unit step length. 

Features TODOs:
- [ ] Render to point set (without drawing lines) via old methods
- [ ] Render with equalized axes (currently the axes scale ratio changes with given image dimension)
- [ ] Switchable backend for matrix multiplication (currently uses Octavian.jl)


Limitations (that probably will stay for a good while):
* Does not support moving without drawing and other common character-level functionalites found in existing L-system implementations
* Does not support arbitrary integer divisions

Internal TODOs:
- [ ] Cleanup multiple dispatch mess with many of the functions
- [ ] Avoid type pirating with certain type alias
- [ ] Represent `charset` with `SortedSet` from DataStructures.jl instead of arrays
- [ ] Define struct `Lsystem` with Parameters.jl
- [ ] Comments!






