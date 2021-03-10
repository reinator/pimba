FROM alpine:3.4

MAINTAINER OLIVEIRA Renato <renato.renison@gmail.com>

RUN apk add --update --no-cache libstdc++ bash perl\
    && apk add --update --no-cache --virtual=deps perl-dev\
    && apk del deps

WORKDIR /data/
COPY bmp-otuName.pl /data/

#CMD ["perl fastx_trimmer", "-h"]
#ENTRYPOINT ["/usr/bin/perl", "/data/bmp-otuName.pl"]