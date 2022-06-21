FROM julia:1.6.3

WORKDIR /home/

COPY . .

RUN julia -E ' \
            using Pkg; \
            Pkg.activate(".");\
            pkg"instantiate" ;\
            pkg"precompile" ;\
            '

RUN echo 'using Pkg; Pkg.activate(".")' > /usr/local/julia/etc/julia/startup.jl 
CMD ["julia"]