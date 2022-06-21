FROM julia:1.6.3

WORKDIR /home/

COPY . .

RUN apt-get -qq update; apt-get -y --no-install-recommends install git qt5-default at-spi2-core libgtk-3-dev xauth xvfb gconf-gsettings-backend
RUN export GKSwstype="100"

RUN julia -E ' \
            using Pkg; \
            Pkg.activate(".");\
            pkg"instantiate" ;\
            pkg"precompile" ;\
            '

RUN echo 'using Pkg; Pkg.activate(".")' > /usr/local/julia/etc/julia/startup.jl 
CMD ["julia"]