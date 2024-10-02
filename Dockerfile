#./Dockerfile
FROM julia:latest

# install git and dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*
 
# authorize SSH host

RUN useradd --create-home --shell /bin/bash jl


USER jl

RUN mkdir /home/jl/RWCode

WORKDIR /home/jl/RWCode

RUN git clone https://github.com/vrg2121/RWCodeData.git

COPY . .

RUN julia --project -e "using Pkg; Pkg.instantiate();"

EXPOSE 8080

ENV JULIA_DEPOT_PATH="/home/jl/.julia"
ENV JULIA_REVISE="off"
ENV EARLYBIND="true"

CMD ["julia", "--project", "./src/RWCode.jl"]