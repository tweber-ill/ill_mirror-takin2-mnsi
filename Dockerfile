FROM ubuntu:22.04
WORKDIR /root

#All build tools should be here (everyone usually already has these)
RUN apt-get update && apt-get -y install wget unzip git build-essential cmake

#Make sure that we have ALL dependencies
##The following are often not available
RUN apt-get --no-install-recommends -y install libboost-dev libboost-filesystem-dev libboost-program-options-dev libboost-iostreams-dev \
	qtbase5-dev liblapacke-dev libpng-dev
##Minuit2 is pretty much never available
RUN git clone https://github.com/GooFit/Minuit2.git && \
	mkdir Minuit2/build && cd Minuit2/build && \
	cmake .. && make && make install

#Compilation
RUN git clone https://github.com/tweber-ill/ill_mirror-takin2-mnsi.git && \
	cd /root/ill_mirror-takin2-mnsi/ext && ./setup_externals.sh && \
	cd /root/ill_mirror-takin2-mnsi && make

WORKDIR /root/ill_mirror-takin2-mnsi/bin
