# Quick Start
This tutorial is for those already familiar with Julia, Git, Docker and / or Amazon AWS EC2.

#### Table of Contents
```@contents
Pages = ["quickstart.md"]
Depth = 2:3
```


## Quick Start using Git
Download data from [RWCodeData](https://github.com/vrg2121/RWCodeData) to RWCodeLocal directory.

!!! tip "Recommended File Structure"
    ```
    RWCodeLocal/
    ├── Data/                   # Data downloaded from RWCodeData
        └── ...                
    ├── Guesses/                # Data downloaded from RWCodeData
        └── ...
    ├── Results/                # Empty folder for saving output results
    └── RWCode/                 # Code from RWCode.jl
        └── src/
            └── functions/
            └── ...     
    ```

Start your Julia session from the CLI with 12+ threads:

```bash
julia +1.10 -t 12
```

Pull the package into the RWCodeLocal directory:
    
```julia
julia> cd("path/to/RWCodeLocal")

julia> VERSION
v"1.10.5"

julia> ] 

(@v1.10) pkg> add https://github.com/vrg2121/RWCode.jl
```


Activate the directory environment:
    
```julia
julia> cd("./RWCode")

julia> ]

pkg> activate .

pkg> instantiate

julia> using RWCode
```


Set up model inputs:
    
- See Model Configuration for more information about `config` options.

```julia
julia> D = "FULL/PATH/TO/DATA"
"FULL/PATH/TO/DATA"

julia> G = "FULL/PATH/TO/Guesses"
"FULL/PATH/TO/Guesses"

julia> R = "FULL/PATH/TO/Results"
"FULL/PATH/TO?Results"

julia> config = ModelConfig()
```


Run the model:
    
```julia
julia> run_rwcode(config, D, G, R)
```


## Quick Start using Docker
Retrieve the docker image and run interactively:
    
```bash
docker run -it vrg2121/rwcode /bin/bash
```
Start up julia with multithreading:

* From the terminal / CLI...
```bash
julia -t 12
```

```julia
julia> ]

pkg> activate .

julia> using RWCode

```

Set up model inputs:
    
```julia
julia> D = "/home/jl/Data"

julia> G = "/home/jl/Guesses"

julia> R = "/home/jl/Results"

julia> config = ModelConfig()

```

!!! note
    The location of the Data, Guesses and Results files are pre-set when using the Docker image.

Run the model:
   
```julia
julia> run_rwcode(config, D, G, R)
```

## Quick Start using Cloud: AWS EC2
Initialize an EC2 instance with at least 8 CPUs.

Load Docker onto the EC2 instance.

Follow the Quick Start using Docker guide.

Retrieve the outputs from EC2.

