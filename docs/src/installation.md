# Installation
RWCode is a model of the global energy grid that has been developed to run on personal computers and the cloud. The code is designed to efficiently parallelize using hyperthreading. It is recommended that you run the model on the cloud or a cluster. Most of the models take less than a day to run on a laptop, but take up significant storage/CPU.

#### Table of Contents
```@contents
Pages = ["installation.md"]
Depth = 2:3
```

## Installation with PC and Git
This installation process assumes that [Julia](https://julialang.org/downloads/) has been downloaded and can be accessed from the command line interface (CLI) or an IDE ([VS Code](https://code.visualstudio.com/docs/languages/julia), [Jupyter Notebooks](https://github.com/JuliaLang/IJulia.jl)). An IDE is not necessary, as everything can be run from the CLI.

### Downloading the data
Due to size constraints, the data for the model is contained in a separate github account.

1. Navigate to the data at [RWCodeData](https://github.com/vrg2121/RWCodeData).
2. Above the list of files, click <> **Code**.
3. Click "Download Zip."

If you are running the code on a laptop or desktop, it is recommended that you unzip these files into an overarching folder for the project. See the recommended file structure below.

### Initialize a Results Folder
In the same directory as the RWCodeData download, create an empty `Results` folder. All data outputs will be saved to this folder.

### Downloading the package
RWCode.jl is not registered and must be pulled from github. It is recommended that you pull it into the same folder as the RWCodeData.

Install Julia from the official binaries [here](https://julialang.org/downloads/#install_julia).


Check that the appropriate version of Julia is installed

```julia
julia> VERSION
v"1.10.5"
```

!!! tip
	If the version does not start with `1.10`, run the following commands in the CLI:
	```bash
	$ juliaup add 1.10.5
	$ julia +1.10
	```


Navigate to the path where RWCode data is stored

```bash
cd PATH/TO/RWCODEDATA
julia +1.10
```

!!! tip
	To exit the Julia REPL in the CLI, just use command `exit()` in the Julia REPL.


Pull the RWCode.jl package from GIT

```julia
julia> ] 
(@v1.10) pkg> add https://github.com/vrg2121/RWCode.jl
```

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

### Set Up Multithreading
Check the number of threads available.

```julia
julia> Threads.nthreads()
1
```

- This is too few threads to run the model efficiently. Need at least 12 threads for the model to run efficiently. More threads will speed up the transition solves.

Set up number of threads available on personal computer.
- Most consumer CPUs support two threads per core. So if your computer has 8 CPUs, you can call up to 16 threads. It is recommended to use at least 12 threads.

- To initiate more threads, execute the following in the terminal:

```bash
julia +1.10 -t 16
```

- You can double check the number of threads available by using the `Threads.nthreads()` command again.

!!! note "Threads in VS Code"
	Refer to [this discourse](https://discourse.julialang.org/t/julia-num-threads-in-vs-code-windows-10-wsl/28794/6) for setting up the number of threads in Julia using the VS Code IDE.

#### Follow the [Tutorial: Run the Model](@ref).

## Installation with Docker
If you do not want to install Julia, download a compatible IDE, download the data and package, then you can use the docker image. The docker image can be used on a personal computer or on the Cloud. The docker image contains the requisite version of Julia

### Running with Docker: Personal Computer
To run RWCode via docker, first set up docker for your computer. See the [dockerdocs](https://docs.docker.com/desktop/?_gl=1*1y6i8pd*_gcl_au*MjExMzU5MjMxMi4xNzI2NzYwMjEy*_ga*MTk3MTkwNDMwOC4xNzI2NzYwMjEy*_ga_XJWPQMJYHQ*MTcyOTA5MjkzMC4xNi4xLjE3MjkwOTI5MzAuNjAuMC4w) for installation intructions for each OS. The docker image contains all of the data, package code and Julia v1.10 necessary for running the code.

Once, docker is started on your computer, execute the following command in your terminal:

```bash
docker run -it vrg2121/rwcode /bin/bash
```

#### Initiate the model run:

Navigate to the Julia REPL with the desired number of threads

```bash
jl@f416a6d0e29c:~/RWCode$ julia -t 12
```
!!! info "Number of Threads"
	There are 2 threads on each CPU available. It is recommended that you run the program with at least 12 threads.

Initiate the package

```julia
julia> ]

(@v1.10) pkg> activate .
Activating project at `~/RWCode`

(RWCode) pkg>
julia> using RWCode
Precompiling RWCode
	Progress [============>                 ] x/120
```

#### Follow the [Tutorial: Run the Model](@ref)

## Installation with Docker on the Cloud
There are various cloud interfaces that can be used to run this model. In fact, Julia has it's own supercomputing provider, [JuliaHub](https://juliahub.com/). 

AWS EC2 is the focus of this tutorial. This tutorial assumes you have already set up an [AWS account](https://aws.amazon.com/resources/create-account/) and launched an [EC2 instance](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html).

### AWS: EC2
1. Install Docker on Your EC2 Instance
2. Pull the docker image from DockerHub
3. Run the model.