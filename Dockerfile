FROM debian:bullseye as prism

RUN apt-get update && apt-get install -y autoconf automake g++ gcc libtool make patch swig wget openjdk-17-jdk-headless \
    && rm -rf /var/lib/apt/lists/* \
    && cd /opt \
    && wget -O- -q https://github.com/prismmodelchecker/prism/archive/refs/tags/v4.7.tar.gz | tar -xzf- \
    && mv /opt/prism-4.7 /opt/prism \
    && cd /opt/prism/prism \
    && make


FROM gradle:jdk17 as gradle

COPY --chown=gradle:gradle sources.tar.gz /opt/sources.tar.gz
COPY --from=prism /opt/prism /opt/lib/models/lib/prism/
RUN cd /opt \
    && tar -xf sources.tar.gz \
    && rm -f /opt/sources.tar.gz \
    && gradle distTar --no-daemon -Pskip-prism-make


FROM debian:bullseye

WORKDIR /opt
RUN apt-get update && apt-get install -y openjdk-17-jre-headless python3 python3-tabulate time \
    && rm -rf /var/lib/apt/lists/*

COPY --from=prism /opt/prism /opt/prism
COPY --from=gradle /opt/build/distributions/stationary-distribution-sampling-0.1.tar /opt/sds/
COPY --from=gradle /opt/models /opt/models
COPY --from=gradle /opt/run.py /opt/eval.py /opt/

RUN cd /opt/sds && tar -xf stationary-distribution-sampling-0.1.tar --strip-components=1 \
    && rm stationary-distribution-sampling-0.1.tar \
    && ln -s /opt/prism/prism/bin/prism /usr/bin/prism \
    && ln -s /opt/sds/bin/stationary-distribution-sampling /usr/bin/sds
