FROM ubuntu 
MAINTAINER <enter your name here> <enter your email here>

RUN apt-get update -y
RUN apt-get upgrade -y

COPY requirements.txt .

RUN DEBIAN_FRONTEND=noninteractive xargs apt-get install -y < requirements.txt 

WORKDIR "/root"

