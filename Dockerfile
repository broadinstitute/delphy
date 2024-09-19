# global-scope variables to persist across image build stages
ARG REPO_DIR="delphy"
ARG REPO_TAG=0.999
ARG BUILD_TYPE="Release" # 'Release' or 'Debug'
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python+deps will not be found)


# ========== [ stage 1: create build environment image ] =======================
FROM --platform=$BUILDPLATFORM mambaorg/micromamba:1.5.9 AS builder

# pull in arguments to the scope of this image build stage
ARG REPO_DIR REPO_TAG BUILD_TYPE
# ensure conda environment is on the PATH
ARG MAMBA_DOCKERFILE_ACTIVATE

# target arguments passed in to the image build stage 
# via the --platform flag (ARG for build-local scope)
#   see: https://www.docker.com/blog/faster-multi-platform-builds-dockerfile-cross-compilation-guide/
ARG TARGETPLATFORM TARGETOS TARGETARCH TARGETVARIANT

COPY --chown=$MAMBA_USER:$MAMBA_USER ./packaging/build-conda-env.yaml ./
# install dependencies via conda (micromamba, in this case)
RUN micromamba install --yes --file ./build-conda-env.yaml && \
    micromamba clean --all --yes


# ========== [ stage 2: compilation ] ==========================================
# For additional information on compilation, see:
#   https://github.com/broadinstitute/delphy/blob/main/INSTALL.md
FROM --platform=$BUILDPLATFORM builder AS compilation

# pull in arguments to the scope of this image build stage
ARG REPO_DIR REPO_TAG BUILD_TYPE
# ensure conda environment is on the PATH
ARG MAMBA_DOCKERFILE_ACTIVATE

ARG TARGETPLATFORM TARGETOS TARGETARCH TARGETVARIANT

# write Python stderr without buffering (does not impact stdout)
ENV PYTHONUNBUFFERED=TRUE

RUN printf "I'm building for TARGETPLATFORM=${TARGETPLATFORM}" \
    && printf ", TARGETARCH=${TARGETARCH}" \
    && printf ", TARGETVARIANT=${TARGETVARIANT} \n" \
    && printf "With uname -s : " && uname -s \
    && printf "and  uname -m : " && uname -m

USER root
RUN apt-get update && apt-get -y install git && apt-get clean
USER $MAMBA_USER

# copy in from full repo instance containing this Dockerfile
COPY --chown=$MAMBA_USER:$MAMBA_USER ./ ${REPO_DIR}
# or...
# make a shallow clone the delphy source tree
# --depth=2 because CheckGit.cmake looks back a commit
#RUN git config --global advice.detachedHead false && \
#    git clone --depth 2 --recursive --shallow-submodules --branch $REPO_TAG \
#    https://github.com/broadinstitute/delphy.git ./${REPO_DIR}

WORKDIR "$REPO_DIR"

RUN mkdir -p build/${BUILD_TYPE}
WORKDIR "build/${BUILD_TYPE}"

# configure conan
RUN conan profile new default --detect && \
    conan profile update settings.compiler.version=13.2 default && \
    conan profile update settings.compiler.libcxx=libstdc++11 default && \
    conan profile update settings.os=$(uname -s) default && \
    conan profile update settings.arch=$(uname -m) default

# check out the linked dependencies that are included as git sub modules
RUN git submodule update --init --recursive

# install build dependencies not
RUN conan install ../..

RUN cmake ../.. -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
                # flags to build as a static binary
                # since the lightweight Alpine Linux image
                # used for distribution favors musl over glibc.
                # Even with static linking, the resulting Alpine-based image
                # is quite reasonable in size
                -DCMAKE_FIND_LIBRARY_SUFFIXES="*.a" \
                -DBUILD_SHARED_LIBS=OFF \
                -DCMAKE_EXE_LINKER_FLAGS="-static"

# compile; use all cores available, reserving one for system housekeeping
RUN make VERBOSE=1 -j $(nproc)

# quick checdk that the binary will run
RUN ./delphy --version | grep "Delphy Version"


# ========== [ stage 3: packing into image for distribution ] ==================
# After building, copy the binary to a minimal image for execution
# (without the bulky compilation toolchain, or the source for linked dependencies)
FROM alpine AS packaging

# It may be an option to base the build on the 'scratch' pseudo-image.
#   see: https://hub.docker.com/_/scratch
# An image built "from scratch" would be extremely minimal, 
# but the result would likely be less useful to users than an Alpine Linux-based
# image, which is slightly larger but also offers the various
# GNU/Linux built-ins helpful for pre- or post- processing of data
# FROM scratch 

# pull in arguments to the scope of this image build stage
ARG REPO_DIR BUILD_TYPE

LABEL software="delphy" \
      about.summary="Delphy is a fast, scalable, accurate and accessible tool \
        for Bayesian phylogenetics based on Explicit Mutation-Annotated Trees" \
      about.home="https://github.com/broadinstitute/delphy"

# copy the binary from the compilation stage to 
COPY --from=compilation --chown=$USER:$USER --chmod=755 \
      /tmp/${REPO_DIR}/build/${BUILD_TYPE}/delphy /usr/local/bin/delphy

# quick checdk that the binary will run
# under the base image used for distribution
RUN /usr/local/bin/delphy --version

CMD ["/usr/local/bin/delphy]