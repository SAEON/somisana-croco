ARG MATLAB_RELEASE=r2022a
FROM mathworks/matlab-deps:${MATLAB_RELEASE}

# Declare the global argument to use at the current build stage
ARG MATLAB_RELEASE

# Install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && \
    apt-get install --no-install-recommends --yes \
        wget \
        unzip \
        ca-certificates && \
    apt-get clean && apt-get autoremove

# Run mpm to install MATLAB in the target location and delete the mpm installation afterwards
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \ 
    chmod +x mpm && \
    ./mpm install \
      --release=${MATLAB_RELEASE} \
      --destination=/opt/matlab \
      --products MATLAB && \
    rm -f mpm /tmp/mathworks_root.log && \
    ln -s /opt/matlab/bin/matlab /usr/local/bin/matlab

# bake in the official matlab croco-tools 
RUN mkdir /croco_tools-v1.3.1 && \
    wget -q https://gitlab.inria.fr/croco-ocean/croco_tools/-/archive/v1.3.1/croco_tools-v1.3.1.tar.gz && \
    tar -xvzf croco_tools-v1.3.1.tar.gz -C /croco_tools-v1.3.1 && \
    rm croco_tools-v1.3.1.tar.gz

# bake in the code from our repo
RUN mkdir /somisana-croco
ADD . /somisana-croco

# Add "matlab" user and grant sudo permission.
RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab \
    && echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab \
    && chmod 0440 /etc/sudoers.d/matlab

ARG LICENSE_SERVER
ENV MLM_LICENSE_FILE=$LICENSE_SERVER
ENV MW_DDUX_FORCE_ENABLE=true MW_CONTEXT_TAGS=MATLAB:DOCKERFILE:V1

# Set user and work directory
USER matlab
WORKDIR /home/matlab
# ENTRYPOINT ["matlab"]
