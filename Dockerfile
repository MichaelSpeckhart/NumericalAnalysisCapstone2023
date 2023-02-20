FROM ubuntu 
MAINTAINER Alex Clevenger (abc324@lehigh.edu)

RUN apt-get update -y
RUN apt-get upgrade -y

COPY requirements.txt .
RUN DEBIAN_FRONTEND=interactive xargs apt-get install -y < requirements.txt 

COPY py_modules.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r py_modules.txt

WORKDIR "/root"

