#./Dockerfile
FROM julia:latest

RUN useradd --create-home --shell /bin/bash jl 

RUN mkdir /home/jl/RWCode

WORKDIR /home/jl/RWCode

RUN chown -R jl:jl /home/jl/

USER jl

COPY . .

RUN julia --project -e "using Pkg; Pkg.instantiate();"

EXPOSE 8080

ENV JULIA_DEPOT_PATH="/home/jl/.julia"
ENV JULIA_REVISE="off"
ENV EARLYBIND="true"

CMD ["julia", "--project", "./src/RWCode.jl"]