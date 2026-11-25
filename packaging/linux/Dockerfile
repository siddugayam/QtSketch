# Build environment for the Linux AppImage: Ubuntu 24.04 (glibc 2.39) with Qt 6.8 LTS from qt.io.
# Used by packaging/linux/build_in_docker.sh.
FROM ubuntu:24.04
ARG QT_VER=6.8.3
ENV DEBIAN_FRONTEND=noninteractive

RUN echo 'Acquire::http::User-Agent "qtsketch-build";' > /etc/apt/apt.conf.d/90useragent

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build git curl ca-certificates file patchelf \
        desktop-file-utils squashfs-tools zsync libfuse2t64 \
        python3 python3-pip python3-venv \
        libgl1-mesa-dev libopengl0 libglx0 libegl1 \
        libxkbcommon0 libxkbcommon-x11-0 libfontconfig1 libfreetype6 libfreetype-dev \
        libdbus-1-3 libglib2.0-0 \
        libxcb-cursor0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 \
        libxcb-randr0 libxcb-render-util0 libxcb-shape0 libxcb-shm0 \
        libxcb-sync1 libxcb-xfixes0 libxcb-xinerama0 libxcb-xkb1 \
        libx11-xcb1 libxrender1 libxi6 libsm6 libice6 \
        libwayland-client0 libwayland-cursor0 libwayland-egl1 \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/aqt && /opt/aqt/bin/pip install --no-cache-dir aqtinstall \
    && if [ "$(uname -m)" = aarch64 ]; then H=linux_arm64; A=linux_gcc_arm64; else H=linux; A=linux_gcc_64; fi \
    && /opt/aqt/bin/aqt install-qt -O /opt/Qt "$H" desktop "$QT_VER" "$A" -m qtshadertools \
    && ln -s /opt/Qt/$QT_VER/gcc_* /opt/qt

ENV PATH=/opt/qt/bin:$PATH CMAKE_PREFIX_PATH=/opt/qt
