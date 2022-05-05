from python:3.9
run pip install numpy
run apt-get update
run apt-get -y install zlib1g 
copy combine_splice.py /bin
#run wget https://storage.googleapis.com/pub/gsutil.tar.gz
#run tar xfz gsutil.tar.gz -C /bin
#ENV PATH="/bin/gsutil:${PATH}"
ENTRYPOINT ["python", "/bin/combine_splice.py"]