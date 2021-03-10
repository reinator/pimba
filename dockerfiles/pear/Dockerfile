FROM itvdsbioinfo/pimba_ubuntu:20.04

MAINTAINER OLIVEIRA Renato <renato.renison@gmail.com>

#INSTALL DEPENDENCIES
RUN apt-get update
RUN apt-get -y install git gcc g++ python2.7 make automake perl
RUN apt-get -y install automake1.11
RUN apt-get -y install autoconf libtool build-essential 
RUN apt-get -y install wget

#RUN git clone https://github.com/tseemann/PEAR.git 
RUN wget https://m1q00q.bl.files.1drv.com/y4m31S8Im8x5u3UACBSPI728HFVaxMsbcWwhSC0gz6Br0nHlFFnIuxxai7Uw5HNjjk7J-r6nDSh5-yfSj7HSIixz8CSKNudJDQpJORNMRCh4g0jVfCvYrdF_rqHzAeXQ8iHyMbtvTju8UkSUHnEC1pMXB6T1Ex2N-2MMmtBzNYgKCg7Kx4UofZDd-2y0Hc_CfkG72viept-oN1wUwIeB6Z1mg
RUN mv y4m31S8Im8x5u3UACBSPI728HFVaxMsbcWwhSC0gz6Br0nHlFFnIuxxai7Uw5HNjjk7J-r6nDSh5-yfSj7HSIixz8CSKNudJDQpJORNMRCh4g0jVfCvYrdF_rqHzAeXQ8iHyMbtvTju8UkSUHnEC1pMXB6T1Ex2N-2MMmtBzNYgKCg7Kx4UofZDd-2y0Hc_CfkG72viept-oN1wUwIeB6Z1mg pear-master.tar.gz
RUN chmod +x pear-master.tar.gz && tar -xzvf pear-master.tar.gz
#RUN tar -xzfv pear-master.tar.gz
RUN cd PEAR-master && ./configure && make && make install
#RUN cd PEAR
#RUN ./configure
#RUN make 
#RUN make install

RUN ln -s /usr/bin/python2.7 /usr/bin/python

WORKDIR /pear/