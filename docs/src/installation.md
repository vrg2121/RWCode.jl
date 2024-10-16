# Installation
RWCode is a model build in Julia that has also been developed to run on the cloud. The code is designed to efficiently parallelize using hyperthreading. It is recommended that you run the model on the cloud or a cluster. Most of the models take less than a day to run on a laptop, but take up significant storage/CPU.

# Standard installation

## Downloading the data
Due to size constraints, the data for the model is contained in a separate github account.

1. Navigate to the data. The data can be found at [RWCodeData](https://github.com/vrg2121/RWCodeData).
2. Above the list of files, click <> **Code**.
3. Click "Download Zip."

If you are running the code on a laptop or desktop, it is recommneded that you unzip these files into an overarching folder for the project.


## Downloading the package
RWCode.jl is not registered and must be pulled from github.

```julia
] add https://github.com/vrg2121/RWCode.jl
```

## Running with Docker
If you do not want to install Julia, start an IDE, download the data and package, then you can use the docker image. The docker image can be used on a personal computer or on the Cloud. 

### Running with Docker: Personal Computer
To run RWCode via docker, first set up docker for your computer. See the [dockerdocs](https://docs.docker.com/desktop/?_gl=1*1y6i8pd*_gcl_au*MjExMzU5MjMxMi4xNzI2NzYwMjEy*_ga*MTk3MTkwNDMwOC4xNzI2NzYwMjEy*_ga_XJWPQMJYHQ*MTcyOTA5MjkzMC4xNi4xLjE3MjkwOTI5MzAuNjAuMC4w) for details for each type of computer. The docker image contains all of the data, package code and Julia v.1.10 necessary for running the code.

Once, docker is started on your computer, execute the following command in your terminal:
```bash
docker run -it vrg2121/rwcode:dev /bin/bash
```

Once the docker image is running, follow these commmands to initiate the model run:

1. Navigate to the Julia REPL
```bash
jl@f416a6d0e29c:~/RWCode$ julia
```

2. Initiate the package
```julia
julia> ]
(@v1.10) pkg> activate .
  Activating project at `~/RWCode`

(RWCode) pkg>
julia> using RWCode
Precompiling RWCode
    Progress [============>                 ] x/118
```

3. Set up vaariable paths to data, guesses and results
```julia
julia> D = "/home/jl/Data"
julia> G = "/home/jl/Guesses"
julia> R = "/home/jl/Results"
```
Note that the location of the Data, Guesses and Results files are pre-set when using the Docker image.

4. Configure the model
```julia
julia> config = ModelConfig()
Enter RunTransition (0 or 1, default = 1):

Enter RunBatteries (0 or 1, default=0):

Enter RunExog (0 or 1, default=0):

Enter RunCurtailment (0 or 1, default=0):

Enter the Number of Transition Iterations (recommend 0-100 iterations, default=2):

Enter Initial Production (default = 100):

Enter hoursofstorage (default=0):

Enter hoursvec (comma-separated, default = 2,4,6):
```
In this example, all of the values are automatically set to the default. For more information about the model configurations see this page.

5. Run the model
```julia
julia> run_rwcode(config, D, G, R)
```
All of the outputs will be labeled and stored as .csv files in the Results file at /home/jl/Results.

6. Copy the results from container to local
Return to the local Terminal.

```bash
docker ps -a
docker cp <containerID>:/home/jl/Results </host/path/target>
```
**Finding the containerID**

Remember to change `</host/path/target>` to the local path you want to copy the data into.

### Running with Docker: Cloud
There are various cloud interfaces that can be used to run this model. In fact, Julia has it's own supercomputing provider, [JuliaHub](https://juliahub.com/). 

While JuliaHub does provide a quick start up and smooth IDE, other services like AWS EC2 are the focus of this tutorial due to lower cost. This tutorial assumes you have already set up an AWS EC2 account. 

#### AWS: EC2